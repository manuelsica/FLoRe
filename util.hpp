#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <string_view>
#include <vector>
using namespace std;

/*
 * encode_kmer_bit:
 * Codifica un k‑mer usando 2 bit per base (A->0, C->1, G->2, T->3).
 */
unsigned int encode_kmer_bit(const string &kmer);

/*
 * sliding_window_read:
 * Segmenta una read in k‑mer utilizzando una finestra scorrevole (step = 1).
 */
vector<string> sliding_window_read(const string &read, int k, int step = 1);

/*
 * reverse_complement:
 * Restituisce il reverse complement di una sequenza nucleotidica.
 */
string reverse_complement(const string &seq);

/*
 * CompressedFingerprint:
 *   - comp_fp: vettore dei codici compressi (utilizzando unsigned int)
 *   - comp_kmers: i k‑mer corrispondenti (non ridondanti) memorizzati come string_view.
 *   - comp_indices: gli indici originali dei k‑mer nella read.
 */
struct CompressedFingerprint {
    vector<unsigned int> comp_fp;
    vector<string_view> comp_kmers;
    vector<int> comp_indices;
};

/*
 * compress_fingerprint:
 * Comprimi la fingerprint eliminando duplicati consecutivi e registrando gli indici.
 */
CompressedFingerprint compress_fingerprint(const vector<unsigned int> &fingerprint, const vector<string_view> &kmers);

/*
 * ProcessedRead:
 *   - comp: fingerprint compressa.
 *   - kmers: lista completa dei k‑mer (string_view).
 *   - fingerprint: vettore completo dei codici (non compresso) usando 32 bit.
 */
struct ProcessedRead {
    CompressedFingerprint comp;
    vector<string_view> kmers;
    vector<unsigned int> fingerprint;
};

/*
 * process_read:
 * Data una read e il valore di k, segmenta in k‑mer usando codifica bit‑level.
 */
ProcessedRead process_read(const string &read, int k);

#endif
