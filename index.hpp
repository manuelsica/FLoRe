// index.hpp
// Header per le funzioni di costruzione degli indici invertiti e del pre-filtraggio delle coppie.

#ifndef INDEX_HPP
#define INDEX_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "util.hpp"
// Includiamo overlap.hpp per avere la definizione completa di SuffixAutomaton
#include "overlap.hpp"

/*
 * Struttura Pair:
 * Rappresenta una coppia di indici (i, j) delle read candidate.
 */
struct Pair
{
    int i, j;
    bool operator==(const Pair &o) const { return (i == o.i && j == o.j); }
};

/*
 * Funzione hash per Pair, necessaria per utilizzare unordered_set.
 */
struct PairHash
{
    size_t operator()(const Pair &p) const
    {
        // Combina gli hash dei due interi in modo semplice
        return ((size_t)p.i * 1315423917ULL) ^ (size_t)p.j;
    }
};

/*
 * Struttura ReadData:
 * Contiene le informazioni processate per ogni read.
 * Include le sequenze (forward e reverse), i ProcessedRead e altre strutture per il pre-filtraggio.
 */
struct ReadData
{
    std::string forward;       // Sequenza originale
    std::string reverse_seq;   // Reverse complement della sequenza

    ProcessedRead pr_fwd;      // ProcessedRead per l'orientamento forward
    ProcessedRead pr_rev;      // ProcessedRead per l'orientamento reverse

    // Strutture per il pre-filtraggio HSIF
    std::unordered_set<long long> set_fwd;  // Set dei fingerprint (forward)
    std::unordered_set<long long> set_rev;  // Set dei fingerprint (reverse)
    std::vector<long long> sorted_fwd;      // Vettore ordinato dei fingerprint (forward)
    std::vector<long long> sorted_rev;      // Vettore ordinato dei fingerprint (reverse)

    // Suffix Automaton, costruito "on demand" se necessario
    SuffixAutomaton suffix_automaton_fwd;
    SuffixAutomaton suffix_automaton_rev;
    bool sa_built_fwd = false;
    bool sa_built_rev = false;
};

/*
 * Funzione buildInvertedIndex:
 * Per ogni read, inserisce i fingerprint (forward e reverse) in due mappe.
 * Le mappe associano ogni fingerprint ad un vettore di indici delle read che lo contengono.
 */
void buildInvertedIndex(const std::vector<ReadData> &all_reads,
                        std::unordered_map<long long, std::vector<int>> &index_fwd,
                        std::unordered_map<long long, std::vector<int>> &index_rev);

/*
 * Funzione generateCandidatePairs:
 * Genera tutte le coppie candidate (senza duplicati) che condividono almeno un fingerprint.
 */
std::vector<Pair> generateCandidatePairs(const std::unordered_map<long long, std::vector<int>> &index_fwd,
                                         const std::unordered_map<long long, std::vector<int>> &index_rev);

/*
 * Funzione hybrid_sorted_intersection_size:
 * Calcola la dimensione dell'intersezione tra due vettori ordinati (utilizzato per il pre-filtraggio HSIF).
 */
int hybrid_sorted_intersection_size(const std::vector<long long> &v1,
                                    const std::vector<long long> &v2);

/*
 * Funzione buildAllReadsData:
 * Processa tutte le read (calcola il fingerprint e il reverse complement) e le salva in ReadData.
 */
void buildAllReadsData(std::vector<ReadData> &all_reads,
                       const std::vector<std::string> &reads,
                       int k,
                       bool use_solid_fingerprint,
                       const std::unordered_set<unsigned int> &solid_fingerprint_set);

/*
 * Funzione fillFingerprintSets:
 * Riempie i set e i vettori ordinati dei fingerprint per ogni ReadData.
 */
void fillFingerprintSets(std::vector<ReadData> &all_reads);

#endif