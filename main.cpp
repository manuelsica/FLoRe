#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>        // Per il parsing delle opzioni da linea di comando
#include "util.hpp"
#include "overlap.hpp"
#include "read.hpp"        // Se necessario, per definizioni aggiuntive
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
#include <optional>        // Per la lazy evaluation

using namespace std;
using namespace std::chrono;

// -----------------------------------------------------------------------------
// Helper per il lazy evaluation: si calcola il reverse complement e il ProcessedRead solo quando necessario.
// -----------------------------------------------------------------------------

// Dichiarazione anticipata per i tipi utilizzati in seguito (la struttura ReadData è definita più avanti).
struct ReadData;

// Restituisce (con caching) il reverse complement della read originale.
const string& get_reverse(const ReadData &rd);

// Restituisce (con caching) il ProcessedRead relativo al reverse complement della read originale.
const ProcessedRead& get_pr_rev(const ReadData &rd, int k);

// -----------------------------------------------------------------------------
// Struttura ReadData: memorizza solo la sequenza originale (forward)
// e calcola on-demand i dati relativi al reverse complement.
// -----------------------------------------------------------------------------
struct ReadData {
    string forward;
    ProcessedRead pr_fwd;
    mutable std::optional<ProcessedRead> pr_rev; // Reverse ProcessedRead calcolato on-demand
    mutable std::optional<string> rev_cache;     // Reverse complement della read (lazy)
    unordered_set<long long> set_fwd;              // Elementi fingerprint unici (forward)
    unordered_set<long long> set_rev;              // Elementi fingerprint unici (reverse)
    vector<long long> sorted_fwd;                  // Vettore ordinato per filtering (forward)
    vector<long long> sorted_rev;                  // Vettore ordinato per filtering (reverse)
};

// -----------------------------------------------------------------------------
// Implementazione delle funzioni helper per il lazy evaluation
// -----------------------------------------------------------------------------
// Calcola e memorizza il reverse complement se non già presente.
const string& get_reverse(const ReadData &rd) {
    if (!rd.rev_cache.has_value()) {
         ReadData &nonConst = const_cast<ReadData&>(rd);
         nonConst.rev_cache = reverse_complement(rd.forward);
    }
    return rd.rev_cache.value();
}

// Calcola e memorizza il ProcessedRead del reverse se non già presente.
const ProcessedRead& get_pr_rev(const ReadData &rd, int k) {
    if (!rd.pr_rev.has_value()) {
         ReadData &nonConst = const_cast<ReadData&>(rd);
         nonConst.pr_rev = process_read(get_reverse(nonConst), k);
    }
    return rd.pr_rev.value();
}

// -----------------------------------------------------------------------------
// Funzioni per il processing e il confronto
// -----------------------------------------------------------------------------
// Ritorna la fingerprint overlap come stringa.
string getFingerprintRegion(const ProcessedRead &pr, int comp_start, int length) {
    ostringstream oss;
    for (int i = comp_start; i < comp_start + length; i++) {
         oss << pr.comp.comp_fp[i];
         if (i < comp_start + length - 1)
              oss << "-";
    }
    return oss.str();
}

// Ritorna un'annotazione in base alla regione e al match fingerprint.
string get_overlap_annotation(const string &region, int fingerprint_match, int min_overlap, int max_repeat_threshold) {
    if (fingerprint_match < min_overlap) {
         return " (SCARTATA)";
    }
    char current = '\0';
    int currentCount = 0, maxCount = 0;
    char maxChar = '\0';
    for (char c : region) {
         if (c == current) {
              currentCount++;
         } else {
              current = c;
              currentCount = 1;
         }
         if (currentCount > maxCount) {
              maxCount = currentCount;
              maxChar = current;
         }
    }
    if (maxCount > max_repeat_threshold) {
         ostringstream oss;
         oss << " (troppi valori consecutivi di \"" << maxChar << "\")";
         return oss.str();
    }
    return "";
}

// Funzione cross-platform per ottenere l'utilizzo massimo di memoria in KB.
size_t getMemoryUsageKB() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if ( GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc)) )
        return pmc.PeakWorkingSetSize / 1024;
    return 0;
#elif defined(__APPLE__)
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss / 1024; // macOS: ru_maxrss in bytes
#else
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // Linux: ru_maxrss in KB
#endif
}

// Calcola l'intersezione fra due vettori ordinati di long long.
int sorted_intersection_size(const vector<long long>& v1, const vector<long long>& v2) {
    int i = 0, j = 0, cnt = 0;
    while (i < v1.size() && j < v2.size()) {
        if (v1[i] < v2[j])
            ++i;
        else if (v2[j] < v1[i])
            ++j;
        else {
            ++cnt;
            ++i;
            ++j;
        }
    }
    return cnt;
}

// Struttura per i risultati, inclusa la fingerprint overlap.
struct OverlapResult {
    int overlap_len;
    int start1;
    int end1;
    int start2;
    int end2;
    string orientation1;
    string orientation2;
    string combination;
    string r1;
    string r2;
    string fingerprint_r1;
    string fingerprint_r2;
};

// Struttura per rappresentare una coppia di read (i < j).
struct Pair {
    int i, j;
    bool operator==(const Pair &other) const { return i == other.i && j == other.j; }
};

// Funzione hash per la struttura Pair.
struct PairHash {
    size_t operator()(const Pair &p) const {
        return ((size_t)p.i * 100000) ^ (size_t)p.j;
    }
};

// -----------------------------------------------------------------------------
// Costruisce gli indici invertiti (forward e reverse).
// -----------------------------------------------------------------------------
void build_inverted_index(vector<ReadData>& all_reads,
                          unordered_map<long long, vector<int>> &index_fwd,
                          unordered_map<long long, vector<int>> &index_rev,
                          int k) {
    for (size_t i = 0; i < all_reads.size(); i++) {
        unordered_set<long long> uniq;
        // Indicizzazione forward
        for (auto val : all_reads[i].pr_fwd.comp.comp_fp)
            uniq.insert(val);
        for (auto val : uniq)
            index_fwd[val].push_back(i);
        uniq.clear();
        // Indicizzazione reverse (lazy)
        const ProcessedRead &pr_rev = get_pr_rev(all_reads[i], k);
        for (auto val : pr_rev.comp.comp_fp)
            uniq.insert(val);
        for (auto val : uniq)
            index_rev[val].push_back(i);
    }
}

// -----------------------------------------------------------------------------
// Genera le coppie candidate utilizzando gli indici invertiti.
// -----------------------------------------------------------------------------
vector<Pair> generate_candidate_pairs_vector(const unordered_map<long long, vector<int>> &index_fwd,
                                               const unordered_map<long long, vector<int>> &index_rev) {
    vector<Pair> candidate_vector;
    // Coppie da index_fwd.
    for (const auto &entry : index_fwd) {
        const vector<int>& vec = entry.second;
        for (size_t i = 0; i < vec.size(); i++) {
            for (size_t j = i + 1; j < vec.size(); j++) {
                candidate_vector.push_back({vec[i], vec[j]});
            }
        }
    }
    // Coppie da index_rev.
    for (const auto &entry : index_rev) {
        const vector<int>& vec = entry.second;
        for (size_t i = 0; i < vec.size(); i++) {
            for (size_t j = i + 1; j < vec.size(); j++) {
                candidate_vector.push_back({vec[i], vec[j]});
            }
        }
    }
    // Incrocio tra index_fwd e index_rev.
    for (const auto &entry : index_fwd) {
        long long code = entry.first;
        if (index_rev.find(code) != index_rev.end()) {
            const vector<int>& vec1 = entry.second;
            const vector<int>& vec2 = index_rev.at(code);
            for (int i : vec1) {
                for (int j : vec2) {
                    if(i == j) continue;
                    candidate_vector.push_back({min(i,j), max(i,j)});
                }
            }
        }
    }
    for (const auto &entry : index_rev) {
        long long code = entry.first;
        if (index_fwd.find(code) != index_fwd.end()) {
            const vector<int>& vec1 = entry.second;
            const vector<int>& vec2 = index_fwd.at(code);
            for (int i : vec1) {
                for (int j : vec2) {
                    if(i == j) continue;
                    candidate_vector.push_back({min(i,j), max(i,j)});
                }
            }
        }
    }
    // Ordina ed elimina duplicati.
    sort(candidate_vector.begin(), candidate_vector.end(), [](const Pair &a, const Pair &b) {
        return (a.i != b.i) ? (a.i < b.i) : (a.j < b.j);
    });
    candidate_vector.erase(unique(candidate_vector.begin(), candidate_vector.end(), [](const Pair &a, const Pair &b) {
        return a.i == b.i && a.j == b.j;
    }), candidate_vector.end());
    return candidate_vector;
}

// -----------------------------------------------------------------------------
// Confronta una coppia di candidate per calcolare il longest common substring
// per le 4 combinazioni (ff, fr, rf, rr).
// -----------------------------------------------------------------------------
OverlapResult compare_candidate_pair(const ReadData &r1, const ReadData &r2, int k) {
    OverlapResult bestResult;
    bestResult.overlap_len = 0;
    int len, comp_idx1, comp_idx2;
    
    // Combinazione ff: r1 forward, r2 forward.
    tie(len, comp_idx1, comp_idx2) = longest_common_substring_suffix_automaton(r1.pr_fwd.comp.comp_fp, r2.pr_fwd.comp.comp_fp);
    if (len > bestResult.overlap_len) {
        bestResult.overlap_len = len;
        bestResult.combination = "ff";
        bestResult.orientation1 = "forward";
        bestResult.orientation2 = "forward";
        bestResult.start1 = r1.pr_fwd.comp.comp_indices[comp_idx1];
        bestResult.end1 = r1.pr_fwd.comp.comp_indices[comp_idx1 + len - 1] + k;
        bestResult.start2 = r2.pr_fwd.comp.comp_indices[comp_idx2];
        bestResult.end2 = r2.pr_fwd.comp.comp_indices[comp_idx2 + len - 1] + k;
        bestResult.r1 = r1.forward;
        bestResult.r2 = r2.forward;
        bestResult.fingerprint_r1 = getFingerprintRegion(r1.pr_fwd, comp_idx1, len);
        bestResult.fingerprint_r2 = getFingerprintRegion(r2.pr_fwd, comp_idx2, len);
    }
    // Combinazione fr: r1 forward, r2 reverse.
    const ProcessedRead &pr_rev_r2 = get_pr_rev(r2, k);
    tie(len, comp_idx1, comp_idx2) = longest_common_substring_suffix_automaton(r1.pr_fwd.comp.comp_fp, pr_rev_r2.comp.comp_fp);
    if (len > bestResult.overlap_len) {
        bestResult.overlap_len = len;
        bestResult.combination = "fr";
        bestResult.orientation1 = "forward";
        bestResult.orientation2 = "reverse";
        bestResult.start1 = r1.pr_fwd.comp.comp_indices[comp_idx1];
        bestResult.end1 = r1.pr_fwd.comp.comp_indices[comp_idx1 + len - 1] + k;
        bestResult.start2 = pr_rev_r2.comp.comp_indices[comp_idx2];
        bestResult.end2 = pr_rev_r2.comp.comp_indices[comp_idx2 + len - 1] + k;
        bestResult.r1 = r1.forward;
        bestResult.r2 = get_reverse(r2);
        bestResult.fingerprint_r1 = getFingerprintRegion(r1.pr_fwd, comp_idx1, len);
        bestResult.fingerprint_r2 = getFingerprintRegion(pr_rev_r2, comp_idx2, len);
    }
    // Combinazione rf: r1 reverse, r2 forward.
    const ProcessedRead &pr_rev_r1 = get_pr_rev(r1, k);
    tie(len, comp_idx1, comp_idx2) = longest_common_substring_suffix_automaton(pr_rev_r1.comp.comp_fp, r2.pr_fwd.comp.comp_fp);
    if (len > bestResult.overlap_len) {
        bestResult.overlap_len = len;
        bestResult.combination = "rf";
        bestResult.orientation1 = "reverse";
        bestResult.orientation2 = "forward";
        bestResult.start1 = pr_rev_r1.comp.comp_indices[comp_idx1];
        bestResult.end1 = pr_rev_r1.comp.comp_indices[comp_idx1 + len - 1] + k;
        bestResult.start2 = r2.pr_fwd.comp.comp_indices[comp_idx2];
        bestResult.end2 = r2.pr_fwd.comp.comp_indices[comp_idx2 + len - 1] + k;
        bestResult.r1 = get_reverse(r1);
        bestResult.r2 = r2.forward;
        bestResult.fingerprint_r1 = getFingerprintRegion(pr_rev_r1, comp_idx1, len);
        bestResult.fingerprint_r2 = getFingerprintRegion(r2.pr_fwd, comp_idx2, len);
    }
    // Combinazione rr: r1 reverse, r2 reverse.
    const ProcessedRead &pr_rev_r1_rr = get_pr_rev(r1, k);
    const ProcessedRead &pr_rev_r2_rr = get_pr_rev(r2, k);
    tie(len, comp_idx1, comp_idx2) = longest_common_substring_suffix_automaton(pr_rev_r1_rr.comp.comp_fp, pr_rev_r2_rr.comp.comp_fp);
    if (len > bestResult.overlap_len) {
        bestResult.overlap_len = len;
        bestResult.combination = "rr";
        bestResult.orientation1 = "reverse";
        bestResult.orientation2 = "reverse";
        bestResult.start1 = pr_rev_r1_rr.comp.comp_indices[comp_idx1];
        bestResult.end1 = pr_rev_r1_rr.comp.comp_indices[comp_idx1 + len - 1] + k;
        bestResult.start2 = pr_rev_r2_rr.comp.comp_indices[comp_idx2];
        bestResult.end2 = pr_rev_r2_rr.comp.comp_indices[comp_idx2 + len - 1] + k;
        bestResult.r1 = get_reverse(r1);
        bestResult.r2 = get_reverse(r2);
        bestResult.fingerprint_r1 = getFingerprintRegion(pr_rev_r1_rr, comp_idx1, len);
        bestResult.fingerprint_r2 = getFingerprintRegion(pr_rev_r2_rr, comp_idx2, len);
    }
    return bestResult;
}

// -----------------------------------------------------------------------------
// Stampa l'uso del programma.
// -----------------------------------------------------------------------------
void print_usage(const char* prog_name) {
    cerr << "Usage: " << prog_name << " -f <file_fasta> [options]" << endl;
    cerr << "Options:" << endl;
    cerr << "  -f, --fasta <file>               File FASTA contenente le reads (obbligatorio)" << endl;
    cerr << "  -m, --min_overlap <int>          Lunghezza minima dell'overlap [default: 13]" << endl;
    cerr << "  -r, --max_repeat_threshold <int> Soglia massima per ripetizioni consecutive [default: 10]" << endl;
    cerr << "  -j, --json <file>                File di output JSON [default: results.json]" << endl;
    cerr << "  -t, --threads <int>              Numero di thread da utilizzare [default: hardware_concurrency]" << endl;
    cerr << "  -v, --verbose                  Stampa i dettagli degli overlap (come in precedenza)" << endl;
    cerr << "  -h, --help                     Mostra questo messaggio di aiuto" << endl;
}

// -----------------------------------------------------------------------------
// Funzione main
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Variabili per le opzioni da linea di comando.
    string fasta_file;
    int min_overlap = 13;
    int max_repeat_threshold = 10;
    string json_filename = "results.json";
    unsigned int num_threads = 0;  // Se 0, verrà impostato in base all'hardware
    bool verbose = false;          // Di default non vengono stampati i dettagli

    // Definizione delle opzioni lunghe.
    struct option long_options[] = {
        {"fasta", required_argument, 0, 'f'},
        {"min_overlap", required_argument, 0, 'm'},
        {"max_repeat_threshold", required_argument, 0, 'r'},
        {"json", required_argument, 0, 'j'},
        {"threads", required_argument, 0, 't'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "f:m:r:j:t:vh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'f':
                fasta_file = optarg;
                break;
            case 'm':
                min_overlap = atoi(optarg);
                break;
            case 'r':
                max_repeat_threshold = atoi(optarg);
                break;
            case 'j':
                json_filename = optarg;
                break;
            case 't':
                num_threads = static_cast<unsigned int>(atoi(optarg));
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (fasta_file.empty()) {
        cerr << "Errore: il file FASTA deve essere specificato." << endl;
        print_usage(argv[0]);
        return 1;
    }

    if (num_threads == 0) {
        num_threads = thread::hardware_concurrency();
        if (num_threads == 0)
            num_threads = 2;
    }

    auto start_time = steady_clock::now();

    cout << "Lettura del file FASTA: " << fasta_file << endl;
    vector<string> reads = read_fasta(fasta_file);
    if (reads.empty()) {
        cerr << "Nessuna read trovata nel file." << endl;
        return 1;
    }
    cout << "Numero di read trovate: " << reads.size() << endl;

    int k = 15;
    size_t N = reads.size();

    cout << "Pre-processing delle read in corso..." << endl;
    vector<ReadData> all_reads(N);
    for (size_t i = 0; i < N; i++) {
        all_reads[i].forward = reads[i];
        all_reads[i].pr_fwd = process_read(reads[i], k);
    }

    cout << "Costruzione degli unordered_set per ciascuna read..." << endl;
    for (size_t i = 0; i < N; i++) {
        for (auto &val : all_reads[i].pr_fwd.comp.comp_fp)
            all_reads[i].set_fwd.insert(val);
        const ProcessedRead &pr_rev = get_pr_rev(all_reads[i], k);
        for (auto &val : pr_rev.comp.comp_fp)
            all_reads[i].set_rev.insert(val);
    }

    cout << "Ordinamento delle fingerprint per migliorare il filtering..." << endl;
    for (size_t i = 0; i < N; i++) {
        all_reads[i].sorted_fwd.assign(all_reads[i].set_fwd.begin(), all_reads[i].set_fwd.end());
        sort(all_reads[i].sorted_fwd.begin(), all_reads[i].sorted_fwd.end());
        all_reads[i].sorted_rev.assign(all_reads[i].set_rev.begin(), all_reads[i].set_rev.end());
        sort(all_reads[i].sorted_rev.begin(), all_reads[i].sorted_rev.end());
    }

    cout << "Costruzione degli indici invertiti..." << endl;
    unordered_map<long long, vector<int>> index_fwd, index_rev;
    build_inverted_index(all_reads, index_fwd, index_rev, k);

    cout << "Generazione delle coppie candidate..." << endl;
    vector<Pair> candidate_vector = generate_candidate_pairs_vector(index_fwd, index_rev);

    cout << "Pre-filtraggio delle coppie candidate..." << endl;
    vector<Pair> filtered_candidates;
    for (const auto &cand : candidate_vector) {
        int common_ff = sorted_intersection_size(all_reads[cand.i].sorted_fwd, all_reads[cand.j].sorted_fwd);
        int common_fr = sorted_intersection_size(all_reads[cand.i].sorted_fwd, all_reads[cand.j].sorted_rev);
        int common_rf = sorted_intersection_size(all_reads[cand.i].sorted_rev, all_reads[cand.j].sorted_fwd);
        int common_rr = sorted_intersection_size(all_reads[cand.i].sorted_rev, all_reads[cand.j].sorted_rev);
        int max_common = max({common_ff, common_fr, common_rf, common_rr});
        if (max_common >= min_overlap)
            filtered_candidates.push_back(cand);
    }
    cout << "Numero di coppie candidate filtrate: " << filtered_candidates.size() << endl;

    cout << "Elaborazione delle candidate con multi-threading (" << num_threads << " thread)..." << endl;
    mutex output_mutex;
    mutex json_mutex;
    vector<string> json_results;

    size_t total = filtered_candidates.size();
    size_t chunk = total / num_threads;
    if (chunk == 0) chunk = 1;

    auto worker = [&](size_t start_idx, size_t end_idx) {
        for (size_t idx = start_idx; idx < end_idx && idx < total; idx++) {
            Pair p = filtered_candidates[idx];
            OverlapResult best = compare_candidate_pair(all_reads[p.i], all_reads[p.j], k);
            if (best.overlap_len > 0) {
                string region_r1 = best.r1.substr(best.start1, best.end1 - best.start1);
                string region_r2 = best.r2.substr(best.start2, best.end2 - best.start2);
                string annotation = get_overlap_annotation(region_r1, best.overlap_len, min_overlap, max_repeat_threshold);
                if (verbose) {
                    lock_guard<mutex> lock(output_mutex);
                    cout << "---------------------------------------" << endl;
                    cout << "Coppia di read " << p.i+1 << " e " << p.j+1 << ":" << endl;
                    cout << "La piu' lunga sottosequenza contigua comune (graph-based su fingerprint) ha lunghezza (in match fingerprint): " 
                         << best.overlap_len << annotation << endl;
                    cout << "Regione di overlap in Read " << p.i+1 << " (" << best.orientation1 << "): " 
                         << region_r1 << endl;
                    cout << "Fingerprint overlap in Read " << p.i+1 << ": " 
                         << best.fingerprint_r1 << endl;
                    cout << "Regione di overlap in Read " << p.j+1 << " (" << best.orientation2 << "): " 
                         << region_r2 << endl;
                    cout << "Fingerprint overlap in Read " << p.j+1 << ": " 
                         << best.fingerprint_r2 << endl;
                    cout << "Tabella riassuntiva:" << endl;
                    cout << "Read1\tRead2\tOrientation\tStart1\tEnd1\tLen(Read1)\tStart2\tEnd2\tLen(Read2)" << endl;
                    int orient_flag = (best.combination == "ff" || best.combination == "rr") ? 0 : 1;
                    cout << p.i+1 << "\t" << p.j+1 << "\t" << orient_flag << "\t"
                         << best.start1 << "\t" << best.end1 << "\t" << best.r1.size() << "\t"
                         << best.start2 << "\t" << best.end2 << "\t" << best.r2.size() << endl;
                    cout << "---------------------------------------" << endl;
                }
                // Salva il risultato in formato JSON
                if (annotation.empty()) {
                    ostringstream oss;
                    oss << "{";
                    oss << "\"read1\": " << p.i+1 << ", ";
                    oss << "\"read2\": " << p.j+1 << ", ";
                    oss << "\"orientation1\": \"" << best.orientation1 << "\", ";
                    oss << "\"orientation2\": \"" << best.orientation2 << "\", ";
                    oss << "\"start1\": " << best.start1 << ", ";
                    oss << "\"end1\": " << best.end1 << ", ";
                    oss << "\"len_read1\": " << best.r1.size() << ", ";
                    oss << "\"start2\": " << best.start2 << ", ";
                    oss << "\"end2\": " << best.end2 << ", ";
                    oss << "\"len_read2\": " << best.r2.size() << ", ";
                    oss << "\"overlap_length\": " << best.overlap_len << ", ";
                    oss << "\"overlap_region_read1\": \"" << region_r1 << "\", ";
                    oss << "\"fingerprint_read1\": \"" << best.fingerprint_r1 << "\", ";
                    oss << "\"overlap_region_read2\": \"" << region_r2 << "\", ";
                    oss << "\"fingerprint_read2\": \"" << best.fingerprint_r2 << "\"";
                    oss << "}";
                    lock_guard<mutex> lock(json_mutex);
                    json_results.push_back(oss.str());
                }
            }
        }
    };

    size_t start_idx = 0;
    vector<thread> pool;
    while (start_idx < total) {
        size_t end_idx = min(total, start_idx + chunk);
        pool.emplace_back(worker, start_idx, end_idx);
        start_idx = end_idx;
    }
    for (auto &t : pool)
        t.join();

    cout << "Salvataggio dei risultati nel file JSON: " << json_filename << endl;
    ofstream json_file(json_filename);
    if (json_file.is_open()) {
        json_file << "[\n";
        for (size_t i = 0; i < json_results.size(); i++) {
            json_file << json_results[i];
            if (i < json_results.size() - 1)
                json_file << ",\n";
        }
        json_file << "\n]\n";
        json_file.close();
    } else {
        cerr << "Errore nell'apertura di " << json_filename << " per la scrittura." << endl;
    }

    if (!verbose)
        cout << "\nTotale overlap scritti in " << json_filename << ": " << json_results.size() << endl;

    auto end_time = steady_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time).count();
    size_t mem_kb = getMemoryUsageKB();

    cout << "\nTempo di esecuzione: " << duration << " ms" << endl;
    cout << "Memoria massima utilizzata: " << mem_kb << " KB" << endl;

    return 0;
}
