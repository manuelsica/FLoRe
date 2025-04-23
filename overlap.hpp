// overlap.hpp
#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <vector>
#include <unordered_map>
#include <tuple>
#include <string>
#include "util.hpp"

struct ReadData;

// Automata di suffissi per fingerprint
struct State {
    int len;
    int link;
    int first_pos;
    std::unordered_map<unsigned int,int> next;
};

struct SuffixAutomaton {
    std::vector<State> st;
    int last;
};

// Risultato dell’overlap
struct OverlapResult {
    int overlap_len = 0;
    int start1 = 0, end1 = 0;
    int start2 = 0, end2 = 0;
    std::string orientation1, orientation2;
    std::string combination;
    std::string r1, r2;
    std::string fingerprint_r1, fingerprint_r2;
    std::string used_algorithm;
};

/*
 * compare_candidate_pair:
 *   Dati due ReadData, calcola il miglior overlap usando
 *   FGOE, AOE, KHS, core e fallback tramite suffix automaton.
 */
OverlapResult compare_candidate_pair(ReadData &r1,
                                     ReadData &r2,
                                     int k,
                                     int min_overlap,
                                     bool verbose,
                                     int max_repeat_threshold);

/*
 * build_suffix_automaton:
 *   Costruisce un suffix automaton su un vettore di fingerprint.
 */
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A);

/*
 * match_suffix_automaton:
 *   Trova longest common substring tra automa e vettore B.
 */
std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B);

/*
 * longest_common_substring_suffix_automaton:
 *   Wrapper che costruisce automa e chiama match_suffix_automaton.
 */
std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int> &A,
    const std::vector<unsigned int> &B);

#endif // OVERLAP_HPP
