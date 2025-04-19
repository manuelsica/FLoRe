// overlap.hpp
// Header per il modulo di calcolo degli overlap.
// Include i metodi FGOE, adaptive_overlap_extension, KHS e core (Combined Overlap Refinement Engine)
// e il fallback tramite Suffix Automaton.
#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <vector>
#include <unordered_map>
#include <tuple>
#include <string>
#include "util.hpp"

// Struttura State per il Suffix Automaton
struct State {
    int len;
    int link;
    int first_pos;
    std::unordered_map<unsigned int, int> next;
};

// Struttura SuffixAutomaton: Rappresenta l'automa per una sequenza.
struct SuffixAutomaton {
    std::vector<State> st;
    int last;
};

// Funzioni per il Suffix Automaton
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A);
std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B);
std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int> &A,
    const std::vector<unsigned int> &B);

// Struttura OverlapResult: Contiene i dettagli dell'overlap.
struct OverlapResult {
    int overlap_len = 0;
    int start1 = 0;
    int end1 = 0;
    int start2 = 0;
    int end2 = 0;
    std::string orientation1;
    std::string orientation2;
    std::string combination;
    std::string r1;
    std::string r2;
    std::string fingerprint_r1;
    std::string fingerprint_r2;
    std::string used_algorithm;
};

/*
 * compare_candidate_pair: Dati due ReadData, esegue i metodi FGOE, adaptive_overlap_extension,
 * KHS e core (Combined Overlap Refinement Engine) per determinare il miglior overlap.
 */
OverlapResult compare_candidate_pair(struct ReadData &r1,
                                     struct ReadData &r2,
                                     int k,
                                     int min_overlap,
                                     bool verbose,
                                     int max_repeat_threshold);

#endif