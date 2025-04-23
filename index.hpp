// index.hpp
#ifndef INDEX_HPP
#define INDEX_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "util.hpp"
#include "overlap.hpp" // Per SuffixAutomaton

// Struttura Pair: Rappresenta una coppia di indici delle read candidate.
struct Pair {
    int i, j;
    bool operator==(const Pair &o) const { return i == o.i && j == o.j; }
};

// Funzione hash per Pair (se mai servisse).
struct PairHash {
    size_t operator()(const Pair &p) const {
        // shift-xor è un semplice mix; non usato nel nuovo codice
        return (static_cast<size_t>(p.i) << 32) ^ static_cast<size_t>(p.j);
    }
};

// Struttura ReadData: contiene sequenze, fingerprint compressi, indici per pre-filtro.
struct ReadData {
    std::string forward;
    std::string reverse_seq;
    ProcessedRead pr_fwd;
    ProcessedRead pr_rev;
    std::unordered_set<long long> set_fwd;
    std::unordered_set<long long> set_rev;
    std::vector<long long> sorted_fwd;
    std::vector<long long> sorted_rev;
    SuffixAutomaton suffix_automaton_fwd;
    SuffixAutomaton suffix_automaton_rev;
    bool sa_built_fwd = false;
    bool sa_built_rev = false;
};

// Costruisce pr_fwd/pr_rev per ogni read (unchanged)
void buildAllReadsData(std::vector<ReadData> &all_reads,
                       const std::vector<std::string> &reads,
                       int k,
                       bool use_solid_fingerprint,
                       const std::unordered_set<unsigned int> &solid_fingerprint_set);

// Riempie i set e sorted_* (unchanged)
void fillFingerprintSets(std::vector<ReadData> &all_reads);

/*
 * buildInvertedIndex:
 *  - Usa sorted_fwd e sorted_rev (già privi di duplicati)
 *  - Riserva spazio basato sul totale delle chiavi
 */
void buildInvertedIndex(const std::vector<ReadData> &all_reads,
                        std::unordered_map<long long, std::vector<int>> &index_fwd,
                        std::unordered_map<long long, std::vector<int>> &index_rev);

/*
 * generateCandidatePairs:
 *  - Accumula in un vector tutte le coppie fwd-fwd, rev-rev e fwd-rev
 *  - Poi ordina + unique per deduplicare
 */
std::vector<Pair> generateCandidatePairs(const std::unordered_map<long long, std::vector<int>> &index_fwd,
                                         const std::unordered_map<long long, std::vector<int>> &index_rev);

// Intersezione ibrida tra due vettori ordinati (unchanged)
int hybrid_sorted_intersection_size(const std::vector<long long> &v1,
                                    const std::vector<long long> &v2);

#endif // INDEX_HPP
