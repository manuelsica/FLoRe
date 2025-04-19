// index.hpp
// Header per le funzioni di costruzione degli indici invertiti e per il pre-filtraggio delle coppie candidate.
#ifndef INDEX_HPP
#define INDEX_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "util.hpp"
#include "overlap.hpp" // Per la definizione completa di SuffixAutomaton

// Struttura Pair: Rappresenta una coppia di indici delle read candidate.
struct Pair {
    int i, j;
    bool operator==(const Pair &o) const { return (i == o.i && j == o.j); }
};

// Funzione hash per la struttura Pair.
struct PairHash {
    size_t operator()(const Pair &p) const {
        return ((size_t)p.i * 1315423917ULL) ^ (size_t)p.j;
    }
};

// Struttura ReadData: Contiene la sequenza originale, il reverse complement,
// i fingerprint e la versione compressa, e strutture per il pre-filtraggio.
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

// Funzione per costruire l'indice invertito.
void buildInvertedIndex(const std::vector<ReadData> &all_reads,
                        std::unordered_map<long long, std::vector<int>> &index_fwd,
                        std::unordered_map<long long, std::vector<int>> &index_rev);

// Funzione per generare le coppie candidate.
std::vector<Pair> generateCandidatePairs(const std::unordered_map<long long, std::vector<int>> &index_fwd,
                                         const std::unordered_map<long long, std::vector<int>> &index_rev);

// Funzione per calcolare l'intersezione tra due vettori ordinati.
int hybrid_sorted_intersection_size(const std::vector<long long> &v1,
                                    const std::vector<long long> &v2);

// Funzione per costruire la struttura ReadData per ogni read.
void buildAllReadsData(std::vector<ReadData> &all_reads,
                       const std::vector<std::string> &reads,
                       int k,
                       bool use_solid_fingerprint,
                       const std::unordered_set<unsigned int> &solid_fingerprint_set);

// Funzione per riempire i set e i vettori ordinati per ogni ReadData.
void fillFingerprintSets(std::vector<ReadData> &all_reads);

#endif