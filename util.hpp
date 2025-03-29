#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
#include <string_view>

/*
 * encode_kmer_bit:
 *   Converte A,C,G,T in 2 bit ciascuno (A->00, C->01, G->10, T->11).
 */
unsigned int encode_kmer_bit(const std::string &kmer);

/*
 * sliding_window_read:
 *   Genera k-mer avanzando di step. Non usata direttamente nel main.
 */
std::vector<std::string> sliding_window_read(const std::string &read, int k, int step = 1);

/*
 * reverse_complement:
 *   Ritorna il reverse complement di una stringa di DNA.
 */
std::string reverse_complement(const std::string &seq);

/*
 * CompressedFingerprint:
 *   comp_fp       : fingerprint compressi (senza ripetizioni consecutive)
 *   comp_kmers    : i k-mer corrispondenti
 *   comp_indices  : indici di partenza di ogni “blocco” compresso
 */
struct CompressedFingerprint {
    std::vector<unsigned int> comp_fp;
    std::vector<std::string_view> comp_kmers;
    std::vector<int> comp_indices;
};

/*
 * ProcessedRead:
 *   fingerprint, kmers, e comp (compressione)
 */
struct ProcessedRead {
    std::vector<unsigned int> fingerprint;
    std::vector<std::string_view> kmers;
    CompressedFingerprint comp;
};

/*
 * compress_fingerprint:
 *   Elimina le ripetizioni consecutive in fingerprint.
 */
CompressedFingerprint compress_fingerprint(const std::vector<unsigned int> &fingerprint,
                                           const std::vector<std::string_view> &kmers);

/*
 * process_read:
 *   Versione standard (in 2 passate).
 */
ProcessedRead process_read(const std::string &read, int k);

/*
 * process_read_optimized (OPKP):
 *   Calcola i fingerprint e la compressione in un’unica passata.
 */
ProcessedRead process_read_optimized(const std::string &read, int k);

#endif
