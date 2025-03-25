#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
using namespace std;

/*
 * compute_val:
 * Calcola il valore numerico del k‑mer usando il mapping:
 *   A → 0, C → 1, G → 2, T → 3,
 * con una rappresentazione in base (default = 5)
 */
long long compute_val(const string &factor, int base = 5);

/*
 * encode_factor:
 * Codifica il k‑mer come intero:
 *   (compute_val(factor) << shift) | len(factor)
 */
long long encode_factor(const string &factor, int shift = 8, int base = 5);

/*
 * sliding_window_read:
 * Segmenta una read in k‑mer utilizzando una finestra scorrevole (step = 1)
 */
vector<string> sliding_window_read(const string &read, int k, int step = 1);

/*
 * reverse_complement:
 * Restituisce il reverse complement di una sequenza nucleotidica.
 */
string reverse_complement(const string &seq);

/*
 * Struttura per memorizzare la fingerprint compressa:
 * - comp_fp: vettore dei codici compressi.
 * - comp_kmers: i k‑mer corrispondenti (non ridondanti).
 * - comp_indices: gli indici originali dei k‑mer nella read completa.
 */
struct CompressedFingerprint {
    vector<long long> comp_fp;
    vector<string> comp_kmers;
    vector<int> comp_indices;
};

/*
 * compress_fingerprint:
 * Comprimi la fingerprint eliminando duplicati consecutivi e registra gli indici.
 */
CompressedFingerprint compress_fingerprint(const vector<long long> &fingerprint, const vector<string> &kmers);

/*
 * Struttura per il risultato dell’elaborazione di una read:
 * - comp: fingerprint compressa (codici, kmers, comp_indices).
 * - kmers: lista completa dei k‑mer.
 * - fingerprint: vettore completo di codici (non compresso).
 * - comp_prefix_mod1, comp_prefix_mod2: prefix hash precomputati sulla fingerprint compressa.
 */
struct ProcessedRead {
    CompressedFingerprint comp;
    vector<string> kmers;
    vector<long long> fingerprint;
    vector<long long> comp_prefix_mod1;
    vector<long long> comp_prefix_mod2;
};

/*
 * process_read:
 * Data una read e il valore di k, segmenta in k‑mer, codifica ogni k‑mer,
 * comprime la fingerprint e precompila i prefix hash sulla fingerprint compressa.
 */
ProcessedRead process_read(const string &read, int k);

#endif
