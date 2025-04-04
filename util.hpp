// util.hpp
// Header per le funzioni utilitarie: encoding dei k-mer, reverse complement,
// compressione dei fingerprint, e altre funzioni di supporto.

#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

/*
 * Funzione encode_kmer_bit:
 * Converte un k-mer (stringa) in un valore unsigned int utilizzando 2 bit per base.
 */
unsigned int encode_kmer_bit(const std::string &kmer);

/*
 * Funzione reverse_complement:
 * Calcola il reverse complement di una sequenza di DNA.
 */
std::string reverse_complement(const std::string &seq);

/*
 * Struttura CompressedFingerprint:
 * Contiene il fingerprint compresso, eliminando ripetizioni consecutive,
 * insieme alle informazioni relative ai k-mer originali e agli indici.
 */
struct CompressedFingerprint {
    std::vector<unsigned int>       comp_fp;       // Vettore dei fingerprint compressi
    std::vector<std::string_view>   comp_kmers;    // Vettore dei k-mer corrispondenti
    std::vector<int>                comp_indices;  // Indici in cui appaiono i k-mer compressi
};

/*
 * Struttura ProcessedRead:
 * Contiene il fingerprint (lista di valori unsigned int), i k-mer (come string_view),
 * e il fingerprint compresso.
 */
struct ProcessedRead {
    std::vector<unsigned int>       fingerprint;
    std::vector<std::string_view>   kmers;
    CompressedFingerprint           comp;
};

/*
 * Funzione compress_fingerprint:
 * Rimuove i k-mer ripetuti consecutivamente dal fingerprint.
 */
CompressedFingerprint compress_fingerprint(const std::vector<unsigned int> &fingerprint,
                                           const std::vector<std::string_view> &kmers);

/*
 * Funzione processReadFingerprint:
 * Calcola il fingerprint e lo comprime in una singola passata.
 */
ProcessedRead processReadFingerprint(const std::string &read, int k);

/*
 * Funzione buildGlobalFingerprintFrequency:
 * Conta la frequenza di ogni k-mer (fingerprint) in tutte le sequenze.
 */
std::unordered_map<unsigned int,int> buildGlobalFingerprintFrequency(const std::vector<std::string> &reads, int k);

/*
 * Funzione processReadSolidFingerprint:
 * Calcola il fingerprint scartando quelli che non sono "solidi"
 * (presenti in un insieme definito di fingerprint solidi).
 */
ProcessedRead processReadSolidFingerprint(const std::string &read,
                                          int k,
                                          const std::unordered_set<unsigned int> &solid_fingerprint_set);

/*
 * Funzione safe_substr:
 * Restituisce una sottostringa in maniera sicura (evita errori out_of_range).
 */
std::string safe_substr(const std::string &s, size_t start, size_t length);

/*
 * Funzione get_overlap_annotation:
 * Genera una stringa di annotazione per un overlap:
 * - Se la lunghezza dell'overlap è inferiore a min_overlap restituisce "(SCARTATA)".
 * - Se ci sono troppe ripetizioni consecutive restituisce un messaggio di scarto.
 * - Altrimenti restituisce una stringa vuota.
 */
std::string get_overlap_annotation(const std::string &region,
                                   int fingerprint_match,
                                   int min_overlap,
                                   int max_repeat_threshold);

#endif