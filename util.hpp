// util.hpp
// Header per le funzioni utilitarie: encoding dei k-mer, reverse complement,
// compressione dei fingerprint, safe_substr, generazione dell'annotazione dell'overlap,
// e costruzione del set di fingerprint solidi.
#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

/*
 * encode_kmer_bit: Converte un k-mer (stringa) in un valore unsigned int usando 2 bit per base.
 */
unsigned int encode_kmer_bit(const std::string &kmer);

/*
 * reverse_complement: Calcola il reverse complement di una sequenza di DNA.
 */
std::string reverse_complement(const std::string &seq);

/*
 * Struttura CompressedFingerprint: Memorizza il fingerprint compresso, i k-mer e gli indici.
 */
struct CompressedFingerprint {
    std::vector<unsigned int>       comp_fp;
    std::vector<std::string_view>   comp_kmers;
    std::vector<int>                comp_indices;
};

/*
 * Struttura ProcessedRead: Contiene il fingerprint, i k-mer e il fingerprint compresso.
 */
struct ProcessedRead {
    std::vector<unsigned int>       fingerprint;
    std::vector<std::string_view>   kmers;
    CompressedFingerprint           comp;
};

/*
 * compress_fingerprint: Elimina i k-mer consecutivi duplicati per ottenere il fingerprint compresso.
 */
CompressedFingerprint compress_fingerprint(const std::vector<unsigned int> &fingerprint,
                                           const std::vector<std::string_view> &kmers);

/*
 * processReadFingerprint: Calcola il fingerprint per una read e lo comprime in una passata.
 */
ProcessedRead processReadFingerprint(const std::string &read, int k);

/*
 * buildGlobalFingerprintFrequency: Conta la frequenza di ogni fingerprint in un vettore di read.
 */
std::unordered_map<unsigned int,int> buildGlobalFingerprintFrequency(const std::vector<std::string> &reads, int k);

/*
 * processReadSolidFingerprint: Processa una read scartando i fingerprint non "solidi".
 */
ProcessedRead processReadSolidFingerprint(const std::string &read,
                                          int k,
                                          const std::unordered_set<unsigned int> &solid_fingerprint_set);

/*
 * safe_substr: Restituisce una sottostringa da una stringa in maniera sicura.
 */
std::string safe_substr(const std::string &s, size_t start, size_t length);

/*
 * get_overlap_annotation: Genera un'annotazione per l'overlap.
 * Se il numero di fingerprint match è inferiore a min_overlap restituisce "(SCARTATA)".
 * Se ci sono troppe ripetizioni consecutive, restituisce un messaggio specifico.
 */
std::string get_overlap_annotation(const std::string &region,
                                   int fingerprint_match,
                                   int min_overlap,
                                   int max_repeat_threshold);

/*
 * buildSolidFingerprintSet: Costruisce un set di fingerprint "solidi" basato sulla frequenza minima e massima.
 */
std::unordered_set<unsigned int> buildSolidFingerprintSet(const std::vector<std::string> &reads,
                                                             int k,
                                                             int min_freq,
                                                             int max_freq);

#endif