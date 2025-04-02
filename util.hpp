#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
#include <string_view>
#include <unordered_set>
#include <unordered_map>

/*
 * encode_kmer_bit:
 *   Converte un k-mer (stringa) in un valore unsigned int (2 bit per base).
 */
unsigned int encode_kmer_bit(const std::string &kmer);

/*
 * sliding_window_read:
 *   Genera k-mer avanzando di step. (Non usato direttamente, ma disponibile)
 */
std::vector<std::string> sliding_window_read(const std::string &read, int k, int step=1);

/*
 * reverse_complement:
 *   Calcola il reverse complement di una stringa di DNA.
 */
std::string reverse_complement(const std::string &seq);

/*
 * Strutture dati:
 *   CompressedFingerprint: fingerprint compresso senza ripetizioni consecutive
 *   ProcessedRead: contiene fingerprint, kmers, e la versione compressa
 */
struct CompressedFingerprint {
    std::vector<unsigned int> comp_fp;
    std::vector<std::string_view> comp_kmers;
    std::vector<int> comp_indices;
};

struct ProcessedRead {
    std::vector<unsigned int> fingerprint;
    std::vector<std::string_view> kmers;
    CompressedFingerprint comp;
};

/*
 * compress_fingerprint:
 *   Rimuove ripetizioni consecutive all’interno del fingerprint.
 */
CompressedFingerprint compress_fingerprint(const std::vector<unsigned int> &fingerprint,
                                           const std::vector<std::string_view> &kmers);

/*
 * process_read:
 *   Versione standard (in 2 passate).
 */
ProcessedRead process_read(const std::string &read, int k);

/*
 * processReadFingerprint:
 *   Ex “process_read_optimized”, esegue fingerprint + compressione in un’unica passata.
 */
ProcessedRead processReadFingerprint(const std::string &read, int k);

/*
 * buildGlobalFingerprintFrequency (ex build_global_kmer_frequency):
 *   Conta la frequenza di tutti i fingerprint (di lunghezza k) nel dataset.
 */
std::unordered_map<unsigned int,int> buildGlobalFingerprintFrequency(const std::vector<std::string> &reads, int k);

/*
 * processReadSolidFingerprint:
 *   Ex “process_read_solid_optimized”; scarta i fingerprint non presenti nel set “solido”.
 */
ProcessedRead processReadSolidFingerprint(const std::string &read,
                                          int k,
                                          const std::unordered_set<unsigned int> &solid_fingerprint_set);

#endif