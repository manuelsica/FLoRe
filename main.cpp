// main.cpp
// Coordina l'intero workflow del tool FLoRe, gestendo:
// - Parsing degli argomenti da linea di comando
// - Avvio del logging asincrono (e salvataggio su file se verbose è attivo)
// - Lettura del file FASTA
// - Costruzione della struttura ReadData (calcolo dei fingerprint, reverse complement, ecc.)
// - Indicizzazione e pre-filtraggio delle coppie candidate
// - Elaborazione parallela per il calcolo degli overlap tramite FGOE, adaptive_overlap_extension, PSH e core (Combined Overlap Refinement Engine)
// - Profiling e benchmarking (tempo e memoria)
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
#ifdef __APPLE__
#include <mach/mach.h>
#endif

#include "logging.hpp"    // Logging asincrono
#include "read.hpp"       // Funzioni per la lettura dei file FASTA
#include "util.hpp"       // Funzioni utilitarie (encoding, reverse complement, ecc.)
#include "index.hpp"      // Indicizzazione e pre-filtraggio delle coppie candidate
#include "overlap.hpp"    // Calcolo degli overlap (FGOE, adaptive_overlap_extension, PSH, core, fallback)
#include "jsonoutput.hpp" // Scrittura dei risultati in formato JSON
#include "profiling.hpp"  // Profiling e benchmarking
#include "filter.hpp"
#include "cfl.hpp"
#include "icfl.hpp"
#include "cfl_icfl.hpp"
using namespace std;

/*
 * print_usage:
 * Stampa le istruzioni d'uso del tool.
 */
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

/* ------------------------------------------------------------------------- */
/* print_config:
 * Riassume tutte le opzioni effettivamente in uso.
 * Va richiamata dopo il parsing degli argomenti.
 */
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
<< "k‑mer length            : " << k                         << '\n'
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
/* ------------------------------------------------------------------------- */


/*
 * getMemoryUsageKB:
 * Restituisce la memoria utilizzata (in KB) dal processo corrente.
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

/*
 * main:
 * Esegue il workflow completo del tool FLoRe.
 */
int main(int argc, char *argv[])
{
    // Avvio del thread di logging asincrono
    thread loggerThread(loggingThreadFunction);

    // Parametri di default
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

    // Parametri per il fingerprint solido
    bool use_solid_fingerprint = false;
    int solid_min_freq = 1;
    int solid_max_freq = 100000;

    // Parsing degli argomenti da linea di comando
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
        case 'f':
            fasta_file = optarg;
            break;
        case 'm':
            min_overlap = atoi(optarg);
            break;
        case 'r':
            max_repeat_threshold = atoi(optarg);
            break;
        case 'k':
            k = atoi(optarg);
            break;
        case 'j':
            json_filename = optarg;
            break;
        case 't':
            num_threads = (unsigned int)atoi(optarg);
            break;
        case 'v':
            verbose = true;
            break;
        case 1000:
            use_solid_fingerprint = true;
            use_cfl = false;
            use_icfl = false;
            break;
        case 1001:
            solid_min_freq = atoi(optarg);
            break;
        case 1002:
            solid_max_freq = atoi(optarg);
            break;
        case 1003:
            cfl_long_threshold = atoi(optarg);
            use_cfl = true;
            break;
        case 1004:
            use_cfl = false;
            use_cflicfl = false;
            use_icfl = false;
            break;
        case 1005: // --icfl_threshold
            icfl_long_threshold = atoi(optarg);
            break;
        case 1006: // --icfl
            use_icfl = true;
            use_cfl = false; // disattiva automaticamente la CFL
            use_cflicfl = false;
            break;
        case 1007: // --cflicfl
            use_cflicfl = true;
            use_cfl = false;
            use_icfl = false;
            break;
        case 1008: // --cflicfl_threshold
            cflicfl_long_threshold = atoi(optarg);
            break;
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
        if (num_threads == 0)
            num_threads = 2;
    }

    print_config(fasta_file, min_overlap, max_repeat_threshold, k,
        json_filename, num_threads, verbose,
        use_solid_fingerprint, solid_min_freq, solid_max_freq,
        use_cfl, cfl_long_threshold, use_icfl, icfl_long_threshold,
        use_cflicfl, cflicfl_long_threshold);
    // Se l'opzione verbose è attiva, abilita il logging su file
    if (verbose)
    {
        logToFile = true;
        logFile.open("verbose.log", std::ios::out | std::ios::trunc);
        if (!logFile.is_open())
        {
            cerr << "Errore nell'apertura del file di log verbose.log" << endl;
        }
    }

    // Avvio del profiling per misurare le performance
    Profiling profiler;
    profiler.start("Inizio workflow FLoRe");
    if (use_cfl)
    {
        async_log("CFL attiva – threshold = " + to_string(cfl_long_threshold) + "\n");
    }
    else if (use_icfl)
    {
        async_log("ICFL attiva – threshold = " + to_string(icfl_long_threshold) + "\n");
    }
    else if (use_cflicfl)
    {
        async_log("CFL→ICFL attiva – threshold = " + to_string(cflicfl_long_threshold) + "\n");
    }
    else
    {
        async_log("CFL/ICFL disattivata – fingerprint classico\n");
    }
    // --- 1) Lettura FASTA ---
    async_log("Lettura del file FASTA: " + fasta_file + "\n");
    vector<string> raw_reads = read_fasta_buffered(fasta_file);
    async_log("MEM_USATA dopo lettura FASTA: " +
              to_string(getMemoryUsageKB() / 1024.0) + " MB\n");
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
    vector<string> processed_reads = preprocess_reads(raw_reads);
    raw_reads.clear();
    raw_reads.shrink_to_fit();
    async_log("MEM_USATA dopo preprocess_reads: " +
              to_string(getMemoryUsageKB() / 1024.0) + " MB\n");
    // --- 2.1) Calcolo lunghezza massima fattori Lyndon di tutti gli algoritmi ---
    // --- 2.1) Calcolo lunghezza massima fattori Lyndon ---
    size_t max_factor = 0;
    if (use_cfl)
    {
        for (auto &r : processed_reads)
        {
            auto fac = cfl_factors(r);
            for (auto &f : fac)
                max_factor = max(max_factor, f.size());
        }
        async_log("Massimo fattore CFL      : " + to_string(max_factor) + "\n");
    }
    else if (use_icfl)
    {
        for (auto &r : processed_reads)
        {
            auto fac = icfl_factors(r);
            for (auto &f : fac)
                max_factor = max(max_factor, f.size());
        }
        async_log("Massimo fattore ICFL     : " + to_string(max_factor) + "\n");
    }
    else if (use_cflicfl)
    {
        for (auto &r : processed_reads)
        {
            auto fac = cfl_icfl(r, cflicfl_long_threshold);
            for (auto &f : fac)
                max_factor = max(max_factor, f.size());
        }
        async_log("Massimo fattore CFL→ICFL : " + to_string(max_factor) + "\n");
    }

    // --- 3) Build ReadData ---
    vector<ReadData> all_reads(processed_reads.size());
    unordered_set<unsigned int> solid_fingerprint_set = use_solid_fingerprint
                                                            ? buildSolidFingerprintSet(processed_reads, k, solid_min_freq, solid_max_freq)
                                                            : unordered_set<unsigned int>();
    buildAllReadsData(all_reads, processed_reads, k,
                      use_solid_fingerprint, solid_fingerprint_set,
                      use_cfl, cfl_long_threshold, use_icfl, icfl_long_threshold, use_cflicfl, cflicfl_long_threshold);
    processed_reads.clear();
    processed_reads.shrink_to_fit();
    async_log("MEM_USATA dopo buildAllReadsData: " +
              to_string(getMemoryUsageKB() / 1024.0) + " MB\n");

    // --- 4) Fill fingerprint sets ---
    fillFingerprintSets(all_reads);
    async_log("MEM_USATA dopo fillFingerprintSets: " +
              to_string(getMemoryUsageKB() / 1024.0) + " MB\n");

    // --- 5) Costruzione indici + candidate pairs ---
    async_log("Costruzione degli indici invertiti...\n");
    unordered_map<long long, vector<int>> index_fwd, index_rev;
    index_fwd.reserve(all_reads.size() * 4);
    index_rev.reserve(all_reads.size() * 4);
    buildInvertedIndex(all_reads, index_fwd, index_rev);

    async_log("Generazione delle coppie candidate...\n");
    vector<Pair> candidate_pairs = generateCandidatePairs(index_fwd, index_rev);
    async_log("MEM_USATA dopo generateCandidatePairs: " +
              to_string(getMemoryUsageKB() / 1024.0) + " MB\n");
    index_fwd.clear();
    index_rev.clear();
    index_fwd.rehash(0);
    index_rev.rehash(0);

    // --- 6) Pre‐filtraggio coppie ---
    async_log("Pre-filtraggio delle coppie candidate...\n");
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
        if (best_val >= min_overlap)
            filtered.push_back(c);
    }
    candidate_pairs.clear();
    candidate_pairs.shrink_to_fit();
    async_log("MEM_USATA dopo filtro coppie: " +
              to_string(getMemoryUsageKB() / 1024.0) + " MB\n");
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
            if (idx >= filtered.size())
                break;
            auto &p = filtered[idx];

            // Overlap classico
            OverlapResult best_ov = compare_candidate_pair(
                all_reads[p.i], all_reads[p.j],
                k, min_overlap, verbose, max_repeat_threshold);
            if (best_ov.overlap_len <= 0)
                continue;

            // Estrazione regioni
            string region_r1 = safe_substr(
                best_ov.r1, best_ov.start1, best_ov.end1 - best_ov.start1);
            string region_r2 = safe_substr(
                best_ov.r2, best_ov.start2, best_ov.end2 - best_ov.start2);

            // Annotazione minima
            string annotation = get_overlap_annotation(
                region_r1, best_ov.overlap_len,
                min_overlap, max_repeat_threshold);

            // FILTRO A: Low complexity
            if (pseudo_overlap::low_complexity(region_r1) ||
                pseudo_overlap::low_complexity(region_r2))
            {

                if (verbose)
                {
                    async_log("-----------------------------------\n");
                    async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                    async_log("Pseudo-overlap scartato: bassa complessità\n");
                    async_log("Region1: " + region_r1 + "\n");
                    async_log("Region2: " + region_r2 + "\n");
                }
                continue;
            }

            // FILTRO B: FCLA
            if (!pseudo_overlap::fingerprint_chained_local_align(
                    region_r1, region_r2, k, 0.80, k * 2))
            {
                if (verbose)
                {
                    async_log("-----------------------------------\n");
                    async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                    async_log("Pseudo-overlap scartato: identità insufficiente (FCLA)\n");
                }
                continue;
            }

            // FILTRO C: Spectrum similarity
            if (!pseudo_overlap::spectrum_similarity(region_r1, region_r2))
            {
                if (verbose)
                {
                    async_log("-----------------------------------\n");
                    async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                    async_log("Pseudo-overlap scartato: spectrum divergence alta\n");
                }
                continue;
            }

            // FILTRO D: Block entropy consistency
            if (!pseudo_overlap::block_entropy_consistency(region_r1) ||
                !pseudo_overlap::block_entropy_consistency(region_r2))
            {
                if (verbose)
                {
                    async_log("-----------------------------------\n");
                    async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                    async_log("Pseudo-overlap scartato: incongruenza entropica a blocchi\n");
                }
                continue;
            }
            if (verbose)
            {
                async_log("-----------------------------------\n");
                async_log("Coppia read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                async_log("Overlap = " + to_string(best_ov.overlap_len) + " " + annotation + "\n");
                async_log("Algoritmo: " + best_ov.used_algorithm + "\n");
                async_log("Orientazioni: " + best_ov.orientation1 + " - " + best_ov.orientation2 + "\n");
                async_log("Region r1: " + region_r1 + "\n");
                async_log("Fingerprint r1: " + best_ov.fingerprint_r1 + "\n");
                async_log("Region r2: " + region_r2 + "\n");
                async_log("Fingerprint r2: " + best_ov.fingerprint_r2 + "\n");
            }
            // Emissione JSON
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
                    << "\"overlap_region_read1\":\"" << region_r1 << "\","
                    << "\"fingerprint_read1\":\"" << best_ov.fingerprint_r1 << "\","
                    << "\"overlap_region_read2\":\"" << region_r2 << "\","
                    << "\"fingerprint_read2\":\"" << best_ov.fingerprint_r2 << "\","
                    << "\"used_algorithm\":\"" << best_ov.used_algorithm << "\""
                    // debug fields for fingerprint mode and marker
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
    for (unsigned int th = 0; th < num_threads; ++th)
        pool.emplace_back(worker);
    for (auto &t : pool)
        t.join();
    filtered.clear();
    filtered.shrink_to_fit();
    profiler.mark("Fine elaborazione overlap");
    size_t total_overlaps = json_results.size();
    async_log("Totale overlap trovati: " + std::to_string(total_overlaps) + "\n");
    // ── INIZIO BLOCCO FASTA ──────────────────────────────────────────────
    {
        // Raccogli tutti gli indici di read coinvolte in un overlap
        std::unordered_set<int> overlap_reads;
        for (auto &jr : json_results)
        {
            overlap_reads.insert(jr.read1 - 1);
            overlap_reads.insert(jr.read2 - 1);
        }
        // Funzione di utilità per serializzare la fingerprint
        auto dump_fp = [&](const std::vector<unsigned int> &fp)
        {
            std::ostringstream os;
            for (size_t j = 0; j < fp.size(); ++j)
            {
                os << fp[j];
                if (j + 1 < fp.size())
                    os << '-';
            }
            return os.str();
        };
        // Apri il “FASTA” di output con le fingerprint
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
    // ── FINE BLOCCO FASTA ────────────────────────────────────────────────

    // --- 8) Salvataggio e fine ---
    auto end_time = chrono::steady_clock::now();
    double total_seconds = chrono::duration<double>(end_time - profiler.start_time()).count();
    size_t mem_kb = getMemoryUsageKB();

    async_log("Salvataggio risultati in " + json_filename + "...\n");
    write_sorted_json(json_results, json_filename);
    async_log("\nTempo totale di esecuzione: " + to_string(total_seconds) + " s\n");
    async_log("Memoria usata (attuale): " + to_string(mem_kb) + " KB\n");
    profiler.stop();

    // Termina logging
    loggingDone.store(true);
    logQueueCV.notify_one();
    loggerThread.join();
    if (logFile.is_open())
        logFile.close();

    return 0;
}
