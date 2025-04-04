// main.cpp
// Questo file contiene la funzione main() che coordina l'esecuzione dell'intero tool FLoRe.
// Include la gestione degli argomenti da linea di comando, il lancio del logging asincrono,
// la lettura del file FASTA, la costruzione degli indici, il calcolo degli overlap in parallelo
// e infine la scrittura dei risultati in formato JSON.

#include <iostream>             // per std::cout, std::cerr
#include <string>               // per std::string
#include <vector>               // per std::vector
#include <getopt.h>             // per la gestione degli argomenti da linea di comando
#include <tuple>                // per std::tuple
#include <sstream>              // per std::ostringstream (costruzione stringhe)
#include <cstdlib>              // per std::atoi
#include <algorithm>            // per std::sort, std::min, std::max
#include <thread>               // per std::thread
#include <mutex>                // per std::mutex, std::lock_guard
#include <fstream>              // per std::ofstream, std::ifstream
#include <unordered_set>        // per std::unordered_set
#include <unordered_map>        // per std::unordered_map
#include <chrono>               // per std::chrono::steady_clock
#include <iterator>             // per std::istreambuf_iterator
#include <cassert>              // per assert()
#include <sys/resource.h>       // per getrusage (Linux)
#include <queue>                // per std::queue
#include <condition_variable>   // per std::condition_variable
#include <atomic>               // per std::atomic
#ifdef __APPLE__
#include <mach/mach.h>          // per task_info su macOS
#endif

// Includiamo i moduli sviluppati
#include "logging.hpp"       // Funzioni e variabili per il logging asincrono
#include "read.hpp"          // Funzioni per leggere file FASTA
#include "util.hpp"          // Funzioni utilitarie (reverse complement, k-mer encoding, ecc.)
#include "index.hpp"         // Funzioni per la costruzione degli indici e per il pre-filtraggio
#include "overlap.hpp"       // Funzioni per il calcolo degli overlap tra le read
#include "jsonoutput.hpp"    // Funzioni per scrivere i risultati in formato JSON

using namespace std;  // per abbreviare i nomi standard (es. std::string -> string)

/*
 * Funzione print_usage:
 * Stampa a stderr il messaggio d'uso del programma.
 */
static void print_usage(const char *prog_name)
{
    cerr << "Usage: " << prog_name << " -f <file_fasta> [options]\n"
         << "  -f, --fasta <file>                  File FASTA da processare\n"
         << "  -m, --min_overlap <int>             Lunghezza minima dell'overlap [default: 13]\n"
         << "  -r, --max_repeat_threshold <int>    Soglia massima per ripetizioni consecutive [default: 10]\n"
         << "  -k, --kmer <int>                    Lunghezza del k-mer [default: 15]\n"
         << "  -j, --json <file>                   File JSON di output [default: results.json]\n"
         << "  -t, --threads <int>                 Numero di thread da utilizzare [default: hardware_concurrency]\n"
         << "  -v, --verbose                       Modalità verbosa (stampa maggiori dettagli)\n"
         << "  --solid_fingerprint                 Abilita la logica di fingerprint solido\n"
         << "  --solid_min_freq <int>              Frequenza minima per fingerprint solido [default: 1]\n"
         << "  --solid_max_freq <int>              Frequenza massima per fingerprint solido [default: 100000]\n"
         << "  -h, --help                          Mostra questo messaggio\n";
}

/*
 * Funzione getMemoryUsageKB:
 * Calcola la memoria usata (in KB) dall'applicazione.
 * Varia a seconda del sistema operativo.
 */
static size_t getMemoryUsageKB()
{
#ifdef _WIN32
    // Su Windows si potrebbe usare GetProcessMemoryInfo
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc, sizeof(pmc)))
    {
        return pmc.PeakWorkingSetSize / 1024;
    }
    return 0;
#elif defined(__APPLE__)
    // Su macOS usiamo le API Mach per ottenere la memoria residente attuale
    mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    kern_return_t kr = task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                                 (task_info_t)&info, &infoCount);
    if (kr != KERN_SUCCESS)
    {
        return 0;
    }
    return info.resident_size / 1024;
#else
    // Su Linux, ru_maxrss (resource usage) indica la memoria usata in KB
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
#endif
}

/*
 * Funzione main:
 * Coordina l'esecuzione del tool FLoRe.
 */
int main(int argc, char *argv[])
{
    // Avvio del thread per il logging asincrono
    thread loggerThread(loggingThreadFunction);

    // Dichiarazione e inizializzazione dei parametri di default
    string fasta_file;
    int min_overlap = 13;
    int max_repeat_threshold = 10;
    int k = 15;
    string json_filename = "results.json";
    unsigned int num_threads = 0;
    bool verbose = false;

    // Parametri per il "solid fingerprint"
    bool use_solid_fingerprint = false;
    int solid_min_freq = 1;
    int solid_max_freq = 100000;

    // Definizione delle opzioni da linea di comando con getopt_long
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
        {0,0,0,0}
    };

    int opt;
    int longindex = 0;
    // Parsing degli argomenti della riga di comando
    while ((opt = getopt_long(argc, argv, "f:m:r:k:j:t:vh", longopts, &longindex)) != -1)
    {
        switch(opt)
        {
        case 'f':
            fasta_file = optarg;  // file FASTA
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
            // Segnala al logging la terminazione
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

    // Controlla che il file FASTA sia stato specificato
    if (fasta_file.empty())
    {
        cerr << "Errore: specificare il file FASTA\n";
        print_usage(argv[0]);
        loggingDone.store(true);
        logQueueCV.notify_one();
        loggerThread.join();
        return 1;
    }
    // Se il numero di thread non è specificato, usa il numero di core hardware
    if (num_threads == 0)
    {
        num_threads = thread::hardware_concurrency();
        if (num_threads == 0)
            num_threads = 2;  // fallback a 2 thread
    }

    // Log iniziale: lettura del file FASTA
    async_log("Lettura del file FASTA: " + fasta_file + "\n");
    vector<string> reads = read_fasta_buffered(fasta_file);
    if (reads.empty())
    {
        cerr << "Nessuna read trovata.\n";
        loggingDone.store(true);
        logQueueCV.notify_one();
        loggerThread.join();
        return 1;
    }
    async_log("Numero di read: " + to_string(reads.size()) + "\n");
    async_log("Fingerprint/k (k-mer): " + to_string(k) + "\n");

    // Se abilitata, costruisce la tabella dei fingerprint solidi
    unordered_set<unsigned int> solid_fingerprint_set;
    if (use_solid_fingerprint)
    {
        async_log("Costruzione tabella di frequenza globale (fingerprint)...\n");
        auto freq_map = buildGlobalFingerprintFrequency(reads, k);
        async_log("Filtraggio fingerprint (min_freq=" + to_string(solid_min_freq) +
                  ", max_freq=" + to_string(solid_max_freq) + ")...\n");
        // Filtra i fingerprint in base alla frequenza
        for (auto &kv : freq_map)
        {
            if (kv.second >= solid_min_freq && kv.second <= solid_max_freq)
            {
                solid_fingerprint_set.insert(kv.first);
            }
        }
        freq_map.clear();
        freq_map.rehash(0);
        async_log("Totale fingerprint solidi: " + to_string(solid_fingerprint_set.size()) + "\n");
    }

    // Costruisce la struttura delle read con i relativi fingerprint
    vector<ReadData> all_reads(reads.size());
    buildAllReadsData(all_reads, reads, k, use_solid_fingerprint, solid_fingerprint_set);

    // Libera la memoria occupata dalle read originali
    reads.clear();
    reads.shrink_to_fit();

    // Prepara le strutture per il pre-filtraggio (HSIF)
    fillFingerprintSets(all_reads);

    // Costruisce gli indici invertiti (OIIB)
    async_log("Costruzione indici invertiti (OIIB)...\n");
    unordered_map<long long, vector<int>> index_fwd, index_rev;
    index_fwd.reserve(all_reads.size() * 4);
    index_rev.reserve(all_reads.size() * 4);
    buildInvertedIndex(all_reads, index_fwd, index_rev);

    // Genera le coppie candidate (HCPG)
    async_log("Generazione coppie candidate (HCPG)...\n");
    vector<Pair> candidate_pairs = generateCandidatePairs(index_fwd, index_rev);

    // Libera la memoria degli indici
    index_fwd.clear();
    index_rev.clear();
    index_fwd.rehash(0);
    index_rev.rehash(0);

    // Pre-filtraggio delle coppie candidate (HSIF)
    async_log("Pre-filtraggio delle coppie candidate (HSIF)...\n");
    vector<Pair> filtered;
    filtered.reserve(candidate_pairs.size());
    // Per ogni coppia, si calcola l'intersezione tra i fingerprint ordinati per decidere se tenere la coppia
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
        {
            filtered.push_back(c);
        }
    }
    candidate_pairs.clear();
    candidate_pairs.shrink_to_fit();
    async_log("Coppie candidate dopo filtro: " + to_string(filtered.size()) + "\n");

    // Avvio del timer per misurare il tempo di elaborazione degli overlap
    using Clock = chrono::steady_clock;
    auto start_time = Clock::now();

    async_log("Elaborazione candidate in multi-thread (" + to_string(num_threads) + " thread)...\n");

    // Prepariamo un vettore per salvare i risultati JSON e una mutex per la sincronizzazione
    mutex json_mutex;
    vector<JsonResult> json_results;
    json_results.reserve(filtered.size());
    // Variabile atomica per indicizzare in maniera concorrente le coppie candidate
    atomic<size_t> candidateIndex(0);

    // Funzione lambda worker per il processing in thread
    auto worker = [&]()
    {
        while(true)
        {
            size_t idx = candidateIndex.fetch_add(1);
            if (idx >= filtered.size())
                break;
            auto &p = filtered[idx];
            // Calcola l'overlap migliore per la coppia
            OverlapResult best_ov = compare_candidate_pair(
                all_reads[p.i],
                all_reads[p.j],
                k,
                min_overlap,
                verbose,
                max_repeat_threshold
            );
            if (best_ov.overlap_len > 0)
            {
                // Estrae le regioni di overlap dalle read
                string region_r1 = safe_substr(best_ov.r1, best_ov.start1,
                                               best_ov.end1 - best_ov.start1);
                string region_r2 = safe_substr(best_ov.r2, best_ov.start2,
                                               best_ov.end2 - best_ov.start2);
                // Calcola un'annotazione basata sui match e ripetizioni
                string annotation = get_overlap_annotation(
                    region_r1,
                    best_ov.overlap_len,
                    min_overlap,
                    max_repeat_threshold
                );
                if (verbose)
                {
                    // Se in modalità verbosa, stampa dettagli
                    async_log("-----------------------------------\n");
                    async_log("Coppia di read " + to_string(p.i+1) + " - " + to_string(p.j+1) + "\n");
                    async_log("Overlap = " + to_string(best_ov.overlap_len) + " " + annotation + "\n");
                    async_log("Algoritmo usato: " + best_ov.used_algorithm + "\n");
                    async_log("Orientation: " + best_ov.orientation1 + " - " + best_ov.orientation2 + "\n");
                    async_log("Region r1: " + region_r1 + "\n");
                    async_log("Fingerprint r1: " + best_ov.fingerprint_r1 + "\n");
                    async_log("Region r2: " + region_r2 + "\n");
                    async_log("Fingerprint r2: " + best_ov.fingerprint_r2 + "\n");
                }
                // Se l'annotazione non indica un match scartato, costruisce l'output JSON
                if (annotation.find("SCARTATA") == string::npos)
                {
                    ostringstream oss;
                    oss << "{";
                    oss << "\"read1\":" << (p.i+1) << ","
                        << "\"read2\":" << (p.j+1) << ","
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
                        << "\"used_algorithm\":\"" << best_ov.used_algorithm << "\"";
                    oss << "}";
                    // Aggiunge il risultato al vettore JSON in modo thread-safe
                    lock_guard<mutex> lk(json_mutex);
                    json_results.push_back(JsonResult{(int)(p.i+1), (int)(p.j+1), oss.str()});
                }
            }
        }
    };

    // Avvio dei thread worker
    vector<thread> pool;
    pool.reserve(num_threads);
    for (unsigned int th = 0; th < num_threads; th++)
    {
        pool.emplace_back(worker);
    }
    for (auto &t : pool)
    {
        t.join();
    }
    filtered.clear();
    filtered.shrink_to_fit();

    // Calcola il tempo totale di elaborazione
    auto end_time = Clock::now();
    double seconds = chrono::duration<double>(end_time - start_time).count();
    size_t mem_kb = getMemoryUsageKB();

    async_log("Salvataggio risultati in " + json_filename + "...\n");
    write_sorted_json(json_results, json_filename);

    if (!verbose)
    {
        async_log("\nTotale overlap scritti: " + to_string(json_results.size()) + "\n");
    }
    async_log("\nTempo di esecuzione: " + to_string(seconds) + " s\n");
    async_log("Memoria usata (attuale): " + to_string(mem_kb) + " KB\n");

    // Segnala al thread di logging di terminare e aspetta il join
    loggingDone.store(true);
    logQueueCV.notify_one();
    loggerThread.join();

    // Libera la memoria delle read processate
    all_reads.clear();
    all_reads.shrink_to_fit();

    return 0;
}