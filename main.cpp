// main.cpp
// Coordina l'intero workflow del tool FLoRe, gestendo:
// - Parsing degli argomenti da linea di comando
// - Avvio del logging asincrono (e salvataggio su file se verbose è attivo)
// - Lettura del file FASTA
// - Costruzione della struttura ReadData (calcolo dei fingerprint, reverse complement, ecc.)
// - Indicizzazione e pre-filtraggio delle coppie candidate
// - Elaborazione parallelo per il calcolo degli overlap tramite FGOE, adaptive_overlap_extension, PSH e core
// - Profiling e benchmarking (tempo e memoria) con tabella salvata su file resume_{modalità}.txt
// - Scrittura dei risultati in formato JSON

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <tuple>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <thread>
#include <mutex>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <iterator>
#include <cassert>
#include <sys/resource.h>
#include <queue>
#include <condition_variable>
#include <atomic>
#include <iomanip>   
#include <climits>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

#include "logging.hpp"
#include "read.hpp"
#include "util.hpp"
#include "index.hpp"
#include "overlap.hpp"
#include "jsonoutput.hpp"
#include "profiling.hpp"
#include "filter.hpp"
#include "cfl.hpp"
#include "icfl.hpp"
#include "cfl_icfl.hpp"
using namespace std;

/* ------------------------------------------------------------------------- */
/* Strumentazione tempi/memoria per funzione/step */

struct StatRec {
    long long ns = 0;
    size_t     calls = 0;
    long long  mem_kb_delta = 0;
};

static std::mutex g_stats_mtx;
static std::unordered_map<std::string, StatRec> g_stats;

static inline void record_stat(const std::string& name, long long ns, long long mem_kb_delta)
{
    std::lock_guard<std::mutex> lk(g_stats_mtx);
    auto &s = g_stats[name];
    s.ns += ns;
    s.calls += 1;
    s.mem_kb_delta += mem_kb_delta;
}

//costruisce la tabella come stringa (riuso per console + file)
static std::string build_stats_table_string()
{
    std::vector<std::pair<std::string, StatRec>> items;
    items.reserve(g_stats.size());
    for (auto &kv : g_stats) items.emplace_back(kv.first, kv.second);

    std::sort(items.begin(), items.end(),
              [](const auto& a, const auto& b){ return a.second.ns > b.second.ns; });

    std::ostringstream out;
    out << "\n================= PROFILING SUMMARY =================\n";
    out << std::left  << std::setw(34) << "Function"
        << std::right << std::setw(10) << "Calls"
        << std::setw(14) << "Total s"
        << std::setw(14) << "Avg ms"
        << std::setw(16) << "MemΔ total KB"
        << std::setw(16) << "MemΔ/call KB" << "\n";
    out << std::string(34+10+14+14+16+16, '-') << "\n";

    for (auto &it : items)
    {
        const auto &name = it.first;
        const auto &s = it.second;
        double total_s = s.ns / 1e9;
        double avg_ms  = s.calls ? (s.ns / 1e6 / (double)s.calls) : 0.0;
        long long mem_total = s.mem_kb_delta;
        long long mem_per_call = s.calls ? (s.mem_kb_delta / (long long)s.calls) : 0;

        out << std::left  << std::setw(34) << name.substr(0,33)
            << std::right << std::setw(10) << s.calls
            << std::setw(14) << std::fixed << std::setprecision(3) << total_s
            << std::setw(14) << std::fixed << std::setprecision(3) << avg_ms
            << std::setw(16) << mem_total
            << std::setw(16) << mem_per_call
            << "\n";
    }
    out << std::string(34+10+14+14+16+16, '-') << "\n";
    return out.str();
}
/* ------------------------------------------------------------------------- */

static void print_usage(const char *prog_name)
{
    cerr << "Usage: " << prog_name << " -f <file_fasta> [options]\n"
         << "  -f, --fasta <file>         File FASTA da processare\n"
         << "  -m, --min_overlap <int>    Lunghezza minima overlap [13]\n"
         << "  -r, --max_repeat_threshold <int> Soglia max ripetizioni consecutive [10]\n"
         << "  -k, --kmer <int>           Lunghezza k-mer [15]\n"
         << "  -j, --json <file>          File JSON di output [results.json]\n"
         << "  -t, --threads <int>        Numero thread [hw_concurrency]\n"
         << "  -v, --verbose              Modalità verbosa (log anche su file)\n"
         << "  --solid_fingerprint        Abilita fingerprint solido\n"
         << "  --solid_min_freq <int>     Frequenza minima [1]\n"
         << "  --solid_max_freq <int>     Frequenza massima [100000]\n"
         << "  --cfl_threshold <int>      Soglia fattori Lyndon lunghi [30]\n"
         << "  --no_cfl                   Disattiva decomposizione CFL\n"
         << "  --icfl_threshold <int>     Soglia fattori Lyndon lunghi [30]\n"
         << "  --icfl                     Abilita decomposizione ICFL\n"
         << "  --cfl_icfl_threshold <int> Soglia fattori CFL per subdecomposizione ICFL [30]\n"
         << "  --cfl_icfl                 Abilita decomposizione CFL+ICFL\n"
         << "  -h, --help                 Mostra questo messaggio\n";
}

static void print_config(const std::string &fasta_file,
    int  min_overlap,
    int  max_repeat_threshold,
    int  k,
    const std::string &json_filename,
    unsigned int num_threads,
    bool verbose,
    bool use_solid_fingerprint,
    int  solid_min_freq,
    int  solid_max_freq,
    bool use_cfl,
    int  cfl_long_threshold,
    bool use_icfl,
    int  icfl_long_threshold,
    bool use_cflicfl,
    int  cflicfl_long_threshold)
{
    auto onoff = [](bool b) { return b ? "ON" : "OFF"; };
    std::ostringstream cfg;

    cfg << "\n================= CONFIG SUMMARY =================\n"
       << "FASTA file              : " << fasta_file                << '\n'
       << "JSON output             : " << json_filename             << '\n'
       << "Threads                 : " << num_threads               << '\n'
       << "Verbose                 : " << onoff(verbose)            << '\n'
       << "--------------------------------------------------\n"
       << "k-mer length            : " << k                         << '\n'
       << "min_overlap             : " << min_overlap               << '\n'
       << "max_repeat_threshold    : " << max_repeat_threshold      << '\n'
       << "--------------------------------------------------\n"
       << "Solid fingerprint       : " << onoff(use_solid_fingerprint) << '\n'
       << "  ├─ min_freq           : " << solid_min_freq            << '\n'
       << "  └─ max_freq           : " << solid_max_freq            << '\n'
       << "CFL                     : " << onoff(use_cfl)            << '\n'
       << "  └─ long_threshold     : " << cfl_long_threshold        << '\n'
       << "ICFL                    : " << onoff(use_icfl)           << '\n'
       << "  └─ long_threshold     : " << icfl_long_threshold       << '\n'
       << "CFL→ICFL                : " << onoff(use_cflicfl)        << '\n'
       << "  └─ long_threshold     : " << cflicfl_long_threshold    << '\n'
       << "==================================================\n";

    async_log(cfg.str());
}

/*
 * getMemoryUsageKB: memoria (KB) del processo.
 */
static size_t getMemoryUsageKB()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc, sizeof(pmc)))
        return pmc.PeakWorkingSetSize / 1024;
    return 0;
#elif defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS)
        return 0;
    return info.resident_size / 1024;
#else
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
#endif
}

int main(int argc, char *argv[])
{
    thread loggerThread(loggingThreadFunction);

    string fasta_file;
    int min_overlap = 13;
    int max_repeat_threshold = 10;
    int k = 15;
    string json_filename = "results.json";
    unsigned int num_threads = 0;
    bool verbose = false;
    bool use_cfl = false;
    int cfl_long_threshold = 30;
    bool use_icfl = false;
    bool use_cflicfl = false;
    int cflicfl_long_threshold = 30;
    int icfl_long_threshold = 30;

    bool use_solid_fingerprint = false;
    int solid_min_freq = 1;
    int solid_max_freq = 100000;

    struct option longopts[] = {
        {"fasta", required_argument, 0, 'f'},
        {"min_overlap", required_argument, 0, 'm'},
        {"max_repeat_threshold", required_argument, 0, 'r'},
        {"kmer", required_argument, 0, 'k'},
        {"json", required_argument, 0, 'j'},
        {"threads", required_argument, 0, 't'},
        {"verbose", no_argument, 0, 'v'},
        {"solid_fingerprint", no_argument, 0, 1000},
        {"solid_min_freq", required_argument, 0, 1001},
        {"solid_max_freq", required_argument, 0, 1002},
        {"cfl_threshold", required_argument, 0, 1003},
        {"no_cfl", no_argument, 0, 1004},
        {"icfl_threshold", required_argument, 0, 1005},
        {"icfl", no_argument, 0, 1006},
        {"cfl_icfl", no_argument, 0, 1007},
        {"cfl_icfl_threshold", required_argument, 0, 1008},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};

    int opt, longindex = 0;
    while ((opt = getopt_long(argc, argv, "f:m:r:k:j:t:vh", longopts, &longindex)) != -1)
    {
        switch (opt)
        {
        case 'f': fasta_file = optarg; break;
        case 'm': min_overlap = atoi(optarg); break;
        case 'r': max_repeat_threshold = atoi(optarg); break;
        case 'k': k = atoi(optarg); break;
        case 'j': json_filename = optarg; break;
        case 't': num_threads = (unsigned int)atoi(optarg); break;
        case 'v': verbose = true; break;
        case 1000: use_solid_fingerprint = true; use_cfl = false; use_icfl = false; break;
        case 1001: solid_min_freq = atoi(optarg); break;
        case 1002: solid_max_freq = atoi(optarg); break;
        case 1003: cfl_long_threshold = atoi(optarg); use_cfl = true; break;
        case 1004: use_cfl = false; use_cflicfl = false; use_icfl = false; break;
        case 1005: icfl_long_threshold = atoi(optarg); break;
        case 1006: use_icfl = true; use_cfl = false; use_cflicfl = false; break;
        case 1007: use_cflicfl = true; use_cfl = false; use_icfl = false; break;
        case 1008: cflicfl_long_threshold = atoi(optarg); break;
        case 'h':
            print_usage(argv[0]);
            loggingDone.store(true);
            logQueueCV.notify_one();
            loggerThread.join();
            return 0;
        default:
            print_usage(argv[0]);
            loggingDone.store(true);
            logQueueCV.notify_one();
            loggerThread.join();
            return 1;
        }
    }
    if (fasta_file.empty())
    {
        cerr << "Errore: specificare il file FASTA di input\n";
        print_usage(argv[0]);
        loggingDone.store(true);
        logQueueCV.notify_one();
        loggerThread.join();
        return 1;
    }
    if (num_threads == 0)
    {
        num_threads = thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;
    }

    print_config(fasta_file, min_overlap, max_repeat_threshold, k,
                 json_filename, num_threads, verbose,
                 use_solid_fingerprint, solid_min_freq, solid_max_freq,
                 use_cfl, cfl_long_threshold, use_icfl, icfl_long_threshold,
                 use_cflicfl, cflicfl_long_threshold);

    if (verbose)
    {
        logToFile = true;
        logFile.open("verbose.log", std::ios::out | std::ios::trunc);
        if (!logFile.is_open())
            cerr << "Errore nell'apertura del file di log verbose.log" << endl;
    }

    Profiling profiler;
    profiler.start("Inizio workflow FLoRe");
    if (use_cfl)            async_log("CFL attiva – threshold = " + to_string(cfl_long_threshold) + "\n");
    else if (use_icfl)      async_log("ICFL attiva – threshold = " + to_string(icfl_long_threshold) + "\n");
    else if (use_cflicfl)   async_log("CFL→ICFL attiva – threshold = " + to_string(cflicfl_long_threshold) + "\n");
    else                    async_log("CFL/ICFL disattivata – fingerprint classico\n");

    // --- 1) Lettura FASTA ---
    async_log("Lettura del file FASTA: " + fasta_file + "\n");
    {
        size_t mem_b = getMemoryUsageKB();
        auto t0 = chrono::steady_clock::now();
        vector<string> raw_reads = read_fasta_buffered(fasta_file);
        auto t1 = chrono::steady_clock::now();
        record_stat("read_fasta_buffered",
                    chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count(),
                    (long long)getMemoryUsageKB() - (long long)mem_b);

        async_log("MEM_USATA dopo lettura FASTA: " + to_string(getMemoryUsageKB() / 1024.0) + " MB\n");
        if (raw_reads.empty())
        {
            cerr << "Nessuna read trovata nel file FASTA.\n";
            profiler.stop();
            loggingDone.store(true);
            logQueueCV.notify_one();
            loggerThread.join();
            return 1;
        }
        async_log("Numero di read lette: " + to_string(raw_reads.size()) + "\n");

        // --- 2) Preprocess reads ---
        size_t mem_b2 = getMemoryUsageKB();
        auto t2 = chrono::steady_clock::now();
        vector<string> processed_reads = preprocess_reads(raw_reads);
        auto t3 = chrono::steady_clock::now();
        record_stat("preprocess_reads",
                    chrono::duration_cast<chrono::nanoseconds>(t3 - t2).count(),
                    (long long)getMemoryUsageKB() - (long long)mem_b2);

        raw_reads.clear();
        raw_reads.shrink_to_fit();
        async_log("MEM_USATA dopo preprocess_reads: " + to_string(getMemoryUsageKB() / 1024.0) + " MB\n");

        // --- 2.1) Fattori Lyndon (CFL/ICFL/CFL→ICFL) ---
        size_t max_factor = 0;
        if (use_cfl)
        {
            size_t mem_b3 = getMemoryUsageKB();
            auto t4 = chrono::steady_clock::now();
            for (auto &r : processed_reads)
            {
                auto fac = cfl_factors(r);
                for (auto &f : fac) max_factor = max(max_factor, f.size());
            }
            auto t5 = chrono::steady_clock::now();
            record_stat("cfl_factors(all)",
                        chrono::duration_cast<chrono::nanoseconds>(t5 - t4).count(),
                        (long long)getMemoryUsageKB() - (long long)mem_b3);
            async_log("Massimo fattore CFL      : " + to_string(max_factor) + "\n");
        }
        else if (use_icfl)
        {
            size_t mem_b3 = getMemoryUsageKB();
            auto t4 = chrono::steady_clock::now();
            for (auto &r : processed_reads)
            {
                auto fac = icfl_factors(r);
                for (auto &f : fac) max_factor = max(max_factor, f.size());
            }
            auto t5 = chrono::steady_clock::now();
            record_stat("icfl_factors(all)",
                        chrono::duration_cast<chrono::nanoseconds>(t5 - t4).count(),
                        (long long)getMemoryUsageKB() - (long long)mem_b3);
            async_log("Massimo fattore ICFL     : " + to_string(max_factor) + "\n");
        }
        else if (use_cflicfl)
        {
            size_t mem_b3 = getMemoryUsageKB();
            auto t4 = chrono::steady_clock::now();
            for (auto &r : processed_reads)
            {
                auto fac = cfl_icfl(r, cflicfl_long_threshold);
                for (auto &f : fac) max_factor = max(max_factor, f.size());
            }
            auto t5 = chrono::steady_clock::now();
            record_stat("cfl_icfl(all)",
                        chrono::duration_cast<chrono::nanoseconds>(t5 - t4).count(),
                        (long long)getMemoryUsageKB() - (long long)mem_b3);
            async_log("Massimo fattore CFL→ICFL : " + to_string(max_factor) + "\n");
        }

        // --- 3) Build ReadData ---
        vector<ReadData> all_reads(processed_reads.size());

        unordered_set<unsigned int> solid_fingerprint_set;
        if (use_solid_fingerprint)
        {
            size_t mem_b4 = getMemoryUsageKB();
            auto t6 = chrono::steady_clock::now();
            solid_fingerprint_set = buildSolidFingerprintSet(processed_reads, k, solid_min_freq, solid_max_freq);
            auto t7 = chrono::steady_clock::now();
            record_stat("buildSolidFingerprintSet",
                        chrono::duration_cast<chrono::nanoseconds>(t7 - t6).count(),
                        (long long)getMemoryUsageKB() - (long long)mem_b4);
        }

        size_t mem_b5 = getMemoryUsageKB();
        auto t8 = chrono::steady_clock::now();
        buildAllReadsData(all_reads, processed_reads, k,
                          use_solid_fingerprint, solid_fingerprint_set,
                          use_cfl, cfl_long_threshold, use_icfl, icfl_long_threshold, use_cflicfl, cflicfl_long_threshold);
        auto t9 = chrono::steady_clock::now();
        record_stat("buildAllReadsData",
                    chrono::duration_cast<chrono::nanoseconds>(t9 - t8).count(),
                    (long long)getMemoryUsageKB() - (long long)mem_b5);

        processed_reads.clear();
        processed_reads.shrink_to_fit();
        async_log("MEM_USATA dopo buildAllReadsData: " + to_string(getMemoryUsageKB() / 1024.0) + " MB\n");

        // --- 4) Fill fingerprint sets ---
        size_t mem_b6 = getMemoryUsageKB();
        auto t10 = chrono::steady_clock::now();
        fillFingerprintSets(all_reads);
        auto t11 = chrono::steady_clock::now();
        record_stat("fillFingerprintSets",
                    chrono::duration_cast<chrono::nanoseconds>(t11 - t10).count(),
                    (long long)getMemoryUsageKB() - (long long)mem_b6);

        async_log("MEM_USATA dopo fillFingerprintSets: " + to_string(getMemoryUsageKB() / 1024.0) + " MB\n");

        // --- 5) Indici + candidate pairs ---
        async_log("Costruzione degli indici invertiti...\n");
        unordered_map<long long, vector<int>> index_fwd, index_rev;
        index_fwd.reserve(all_reads.size() * 4);
        index_rev.reserve(all_reads.size() * 4);

        size_t mem_b7 = getMemoryUsageKB();
        auto t12 = chrono::steady_clock::now();
        buildInvertedIndex(all_reads, index_fwd, index_rev);
        auto t13 = chrono::steady_clock::now();
        record_stat("buildInvertedIndex",
                    chrono::duration_cast<chrono::nanoseconds>(t13 - t12).count(),
                    (long long)getMemoryUsageKB() - (long long)mem_b7);

        async_log("Generazione delle coppie candidate...\n");
        size_t mem_b8 = getMemoryUsageKB();
        auto t14 = chrono::steady_clock::now();
        vector<Pair> candidate_pairs = generateCandidatePairs(index_fwd, index_rev);
        auto t15 = chrono::steady_clock::now();
        record_stat("generateCandidatePairs",
                    chrono::duration_cast<chrono::nanoseconds>(t15 - t14).count(),
                    (long long)getMemoryUsageKB() - (long long)mem_b8);

        async_log("MEM_USATA dopo generateCandidatePairs: " + to_string(getMemoryUsageKB() / 1024.0) + " MB\n");
        index_fwd.clear(); index_rev.clear();
        index_fwd.rehash(0); index_rev.rehash(0);

        // --- 6) Pre‐filtraggio coppie ---
        async_log("Pre-filtraggio delle coppie candidate...\n");
        size_t mem_b9 = getMemoryUsageKB();
        auto t16 = chrono::steady_clock::now();
        vector<Pair> filtered;
        filtered.reserve(candidate_pairs.size());
        for (auto &c : candidate_pairs)
        {
            auto &rd1 = all_reads[c.i];
            auto &rd2 = all_reads[c.j];
            int ff = hybrid_sorted_intersection_size(rd1.sorted_fwd, rd2.sorted_fwd);
            int fr = hybrid_sorted_intersection_size(rd1.sorted_fwd, rd2.sorted_rev);
            int rf = hybrid_sorted_intersection_size(rd1.sorted_rev, rd2.sorted_fwd);
            int rr = hybrid_sorted_intersection_size(rd1.sorted_rev, rd2.sorted_rev);
            int best_val = max({ff, fr, rf, rr});
            if (best_val >= min_overlap) filtered.push_back(c);
        }
        auto t17 = chrono::steady_clock::now();
        record_stat("prefilter_pairs",
                    chrono::duration_cast<chrono::nanoseconds>(t17 - t16).count(),
                    (long long)getMemoryUsageKB() - (long long)mem_b9);

        candidate_pairs.clear(); candidate_pairs.shrink_to_fit();
        async_log("MEM_USATA dopo filtro coppie: " + to_string(getMemoryUsageKB() / 1024.0) + " MB\n");
        async_log("Coppie candidate dopo filtro: " + to_string(filtered.size()) + "\n");

        // --- 7) Overlap parallelo ---
        profiler.mark("Inizio elaborazione overlap");
        mutex json_mutex;
        vector<JsonResult> json_results;
        json_results.reserve(filtered.size());
        atomic<size_t> candidateIndex(0);

        auto worker = [&]()
        {
            while (true)
            {
                size_t idx = candidateIndex.fetch_add(1);
                if (idx >= filtered.size()) break;
                auto &p = filtered[idx];

                size_t mem_bw = getMemoryUsageKB();
                auto tw0 = chrono::steady_clock::now();
                OverlapResult best_ov = compare_candidate_pair(
                    all_reads[p.i], all_reads[p.j],
                    k, min_overlap, verbose, max_repeat_threshold);
                auto tw1 = chrono::steady_clock::now();
                record_stat("compare_candidate_pair",
                            chrono::duration_cast<chrono::nanoseconds>(tw1 - tw0).count(),
                            (long long)getMemoryUsageKB() - (long long)mem_bw);

                if (best_ov.overlap_len <= 0) continue;

                string region_r1 = safe_substr(best_ov.r1, best_ov.start1, best_ov.end1 - best_ov.start1);
                string region_r2 = safe_substr(best_ov.r2, best_ov.start2, best_ov.end2 - best_ov.start2);

                string annotation = get_overlap_annotation(region_r1, best_ov.overlap_len, min_overlap, max_repeat_threshold);

                // FILTRO A: Low complexity
                {
                    size_t mem_bA = getMemoryUsageKB();
                    auto tA0 = chrono::steady_clock::now();
                    bool lc = (pseudo_overlap::low_complexity(region_r1) || pseudo_overlap::low_complexity(region_r2));
                    auto tA1 = chrono::steady_clock::now();
                    record_stat("low_complexity",
                                chrono::duration_cast<chrono::nanoseconds>(tA1 - tA0).count(),
                                (long long)getMemoryUsageKB() - (long long)mem_bA);
                    if (lc) { if (verbose) { async_log("-----------------------------------\n");
                            async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                            async_log("Pseudo-overlap scartato: bassa complessità\n"); }
                            continue; }
                }

                // FILTRO B: FCLA
                {
                    size_t mem_bB = getMemoryUsageKB();
                    auto tB0 = chrono::steady_clock::now();
                    bool ok = pseudo_overlap::fingerprint_chained_local_align(region_r1, region_r2, k, 0.80, k * 2);
                    auto tB1 = chrono::steady_clock::now();
                    record_stat("FCLA",
                                chrono::duration_cast<chrono::nanoseconds>(tB1 - tB0).count(),
                                (long long)getMemoryUsageKB() - (long long)mem_bB);
                    if (!ok) { if (verbose) { async_log("-----------------------------------\n");
                            async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                            async_log("Pseudo-overlap scartato: identità insufficiente (FCLA)\n"); }
                            continue; }
                }

                // FILTRO C: Spectrum similarity
                {
                    size_t mem_bC = getMemoryUsageKB();
                    auto tC0 = chrono::steady_clock::now();
                    bool ok = pseudo_overlap::spectrum_similarity(region_r1, region_r2);
                    auto tC1 = chrono::steady_clock::now();
                    record_stat("spectrum_similarity",
                                chrono::duration_cast<chrono::nanoseconds>(tC1 - tC0).count(),
                                (long long)getMemoryUsageKB() - (long long)mem_bC);
                    if (!ok) { if (verbose) { async_log("-----------------------------------\n");
                            async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                            async_log("Pseudo-overlap scartato: spectrum divergence alta\n"); }
                            continue; }
                }

                // FILTRO D: Block entropy consistency
                {
                    size_t mem_bD = getMemoryUsageKB();
                    auto tD0 = chrono::steady_clock::now();
                    bool ok = (pseudo_overlap::block_entropy_consistency(region_r1) &&
                               pseudo_overlap::block_entropy_consistency(region_r2));
                    auto tD1 = chrono::steady_clock::now();
                    record_stat("block_entropy_consistency",
                                chrono::duration_cast<chrono::nanoseconds>(tD1 - tD0).count(),
                                (long long)getMemoryUsageKB() - (long long)mem_bD);
                    if (!ok) { if (verbose) { async_log("-----------------------------------\n");
                            async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                            async_log("Pseudo-overlap scartato: incongruenza entropica a blocchi\n"); }
                            continue; }
                }

                if (verbose)
                {
                    async_log("-----------------------------------\n");
                    async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                    async_log("Overlap = " + to_string(best_ov.overlap_len) + " " + annotation + "\n");
                    async_log("Algoritmo: " + best_ov.used_algorithm + "\n");
                }

                if (annotation.find("SCARTATA") == string::npos)
                {
                    ostringstream oss;
                    oss << "{"
                        << "\"read1\":" << (p.i + 1) << ","
                        << "\"read2\":" << (p.j + 1) << ","
                        << "\"orientation1\":\"" << best_ov.orientation1 << "\","
                        << "\"orientation2\":\"" << best_ov.orientation2 << "\","
                        << "\"start1\":" << best_ov.start1 << ","
                        << "\"end1\":" << best_ov.end1 << ","
                        << "\"len_read1\":" << best_ov.r1.size() << ","
                        << "\"start2\":" << best_ov.start2 << ","
                        << "\"end2\":" << best_ov.end2 << ","
                        << "\"len_read2\":" << best_ov.r2.size() << ","
                        << "\"overlap_length\":" << best_ov.overlap_len << ","
                        << "\"overlap_region_read1\":\"" << safe_substr(best_ov.r1, best_ov.start1, best_ov.end1 - best_ov.start1) << "\","
                        << "\"fingerprint_read1\":\"" << best_ov.fingerprint_r1 << "\","
                        << "\"overlap_region_read2\":\"" << safe_substr(best_ov.r2, best_ov.start2, best_ov.end2 - best_ov.start2) << "\","
                        << "\"fingerprint_read2\":\"" << best_ov.fingerprint_r2 << "\","
                        << "\"used_algorithm\":\"" << best_ov.used_algorithm << "\""
                        << ",\"fp_type\":\""
                        << (use_icfl                ? "ICFL"
                            : use_cfl               ? "CFL"
                            : use_solid_fingerprint ? "SOLID"
                            : use_cflicfl           ? "CFL_ICFL"
                                                    : "CLASSIC")
                        << "\"}";
                    lock_guard<mutex> lk(json_mutex);
                    json_results.emplace_back(p.i + 1, p.j + 1, oss.str());
                }
            }
        };

        vector<thread> pool;
        pool.reserve(num_threads);
        for (unsigned int th = 0; th < num_threads; ++th) pool.emplace_back(worker);
        for (auto &t : pool) t.join();
        filtered.clear(); filtered.shrink_to_fit();
        profiler.mark("Fine elaborazione overlap");
        size_t total_overlaps = json_results.size();
        async_log("Totale overlap trovati: " + std::to_string(total_overlaps) + "\n");

        // ── FASTA di fingerprint delle read coinvolte ──────────────────────
        {
            std::unordered_set<int> overlap_reads;
            for (auto &jr : json_results) { overlap_reads.insert(jr.read1 - 1); overlap_reads.insert(jr.read2 - 1); }
            auto dump_fp = [&](const std::vector<unsigned int> &fp){
                std::ostringstream os; for (size_t j = 0; j < fp.size(); ++j) { os << fp[j]; if (j + 1 < fp.size()) os << '-'; }
                return os.str();
            };
            std::ofstream fa_out("overlap_reads.fa");
            if (fa_out)
            {
                for (int i : overlap_reads)
                {
                    fa_out << ">read" << (i + 1) << "_fp_forward\n"
                           << dump_fp(all_reads[i].pr_fwd.comp.comp_fp) << "\n"
                           << ">read" << (i + 1) << "_fp_reverse\n"
                           << dump_fp(all_reads[i].pr_rev.comp.comp_fp) << "\n";
                }
                fa_out.close();
                async_log("Salvate " + std::to_string(overlap_reads.size()) + " fingerprint (forward+reverse) in overlap_reads.fa\n");
            }
        }

        // --- 8) Salvataggio e fine ---
        auto end_time = chrono::steady_clock::now();
        double total_seconds = chrono::duration<double>(end_time - profiler.start_time()).count();
        size_t mem_kb = getMemoryUsageKB();

        async_log("Salvataggio risultati in " + json_filename + "...\n");
        {
            size_t mem_b_json = getMemoryUsageKB();
            auto tj0 = chrono::steady_clock::now();
            write_sorted_json(json_results, json_filename);
            auto tj1 = chrono::steady_clock::now();
            record_stat("write_sorted_json",
                        chrono::duration_cast<chrono::nanoseconds>(tj1 - tj0).count(),
                        (long long)getMemoryUsageKB() - (long long)mem_b_json);
        }

        async_log("\nTempo totale di esecuzione: " + to_string(total_seconds) + " s\n");
        async_log("Memoria usata (attuale): " + to_string(mem_kb) + " KB\n");
        profiler.stop();

        // Tabella riassuntiva -> console + file resume_{modalità}.txt
        std::string table = build_stats_table_string();
        std::cout << table;

        // stessa priorità di fp_type: ICFL > CFL > SOLID > CFL_ICFL > CLASSIC (lowercase per filename)
        std::string mode =
            (use_icfl                ? "icfl" :
             use_cfl                ? "cfl"  :
             use_solid_fingerprint  ? "solid":
             use_cflicfl            ? "cfl_icfl" :
                                       "classic");

        std::string fname = "resume_" + mode + ".txt";
        std::ofstream resume_out(fname);
        if (resume_out) {
            resume_out << table;
            resume_out.close();
            async_log("Tabella di profiling salvata in " + fname + "\n");
        } else {
            async_log("Impossibile scrivere il file di riepilogo: " + fname + "\n");
        }

        // cleanup logging
        loggingDone.store(true);
        logQueueCV.notify_one();
        loggerThread.join();
        if (logFile.is_open()) logFile.close();

        return 0;
    }
}
