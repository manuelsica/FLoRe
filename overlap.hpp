#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <vector>
#include <unordered_map>
#include <tuple>

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

SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A);
std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B);

std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int> &A,
    const std::vector<unsigned int> &B
);

#endif