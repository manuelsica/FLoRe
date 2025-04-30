// main.cpp
// Coordina l'intero workflow del tool FLoRe, gestendo:
// - Parsing degli argomenti da linea di comando
// - Avvio del logging asincrono (e salvataggio su file se verbose è attivo)
// - Lettura del file FASTA
// - Costruzione della struttura ReadData (calcolo dei fingerprint, reverse complement, ecc.)
// - Indicizzazione e pre-filtraggio delle coppie candidate
// - Elaborazione parallela per il calcolo degli overlap tramite FGOE, adaptive_overlap_extension, KHS e core (Combined Overlap Refinement Engine)
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
#include "overlap.hpp"    // Calcolo degli overlap (FGOE, adaptive_overlap_extension, KHS, core, fallback)
#include "jsonoutput.hpp" // Scrittura dei risultati in formato JSON
#include "profiling.hpp"  // Profiling e benchmarking
#include "filter.hpp"

using namespace std;

/*
 * print_usage:
 * Stampa le istruzioni d'uso del tool.
 */
static void print_usage(const char *prog_name)
{
    cerr << "Usage: " << prog_name << " -f <file_fasta> [options]\n"
         << "  -f, --fasta <file>         File FASTA da processare\n"
         << "  -m, --min_overlap <int>    Lunghezza minima dell'overlap [default: 13]\n"
         << "  -r, --max_repeat_threshold <int> Soglia massima per ripetizioni consecutive [default: 10]\n"
         << "  -k, --kmer <int>           Lunghezza del k-mer [default: 15]\n"
         << "  -j, --json <file>          File JSON di output [default: results.json]\n"
         << "  -t, --threads <int>        Numero di thread [default: hardware_concurrency]\n"
         << "  -v, --verbose              Modalità verbosa (salva log anche in 'verbose.log')\n"
         << "  --solid_fingerprint        Abilita logica di fingerprint solido\n"
         << "  --solid_min_freq <int>     Frequenza minima per fingerprint solido [default: 1]\n"
         << "  --solid_max_freq <int>     Frequenza massima per fingerprint solido [default: 100000]\n"
         << "  -h, --help                 Mostra questo messaggio\n";
}

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
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

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
            break;
        case 1001:
            solid_min_freq = atoi(optarg);
            break;
        case 1002:
            solid_max_freq = atoi(optarg);
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

    // --- 3) Build ReadData ---
    vector<ReadData> all_reads(processed_reads.size());
    unordered_set<unsigned int> solid_fingerprint_set = use_solid_fingerprint
        ? buildSolidFingerprintSet(processed_reads, k, solid_min_freq, solid_max_freq)
        : unordered_set<unsigned int>();
    buildAllReadsData(all_reads, processed_reads, k, use_solid_fingerprint, solid_fingerprint_set);
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
            if (idx >= filtered.size()) break;
            auto &p = filtered[idx];

            // Overlap classico
            OverlapResult best_ov = compare_candidate_pair(
                all_reads[p.i], all_reads[p.j],
                k, min_overlap, verbose, max_repeat_threshold);
            if (best_ov.overlap_len <= 0) continue;

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
                pseudo_overlap::low_complexity(region_r2)) continue;

            // FILTRO B: FCLA
            if (!pseudo_overlap::fingerprint_chained_local_align(
                region_r1, region_r2, k, 0.80, k * 2)) continue;

            // FILTRO C: Spectrum similarity
            if (!pseudo_overlap::spectrum_similarity(region_r1, region_r2)) continue;

            // FILTRO D: Block entropy consistency
            if (!pseudo_overlap::block_entropy_consistency(region_r1) ||
                !pseudo_overlap::block_entropy_consistency(region_r2)) continue;

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
                    << "}";
                lock_guard<mutex> lk(json_mutex);
                json_results.emplace_back(p.i + 1, p.j + 1, oss.str());
            }
        }
    };

    vector<thread> pool;
    pool.reserve(num_threads);
    for (unsigned int th = 0; th < num_threads; ++th)
        pool.emplace_back(worker);
    for (auto &t : pool) t.join();
    filtered.clear();
    filtered.shrink_to_fit();
    profiler.mark("Fine elaborazione overlap");

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
    if (logFile.is_open()) logFile.close();

    return 0;
}
