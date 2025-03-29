#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <vector>
#include <unordered_map>
#include <tuple>

/*
 * State:
 *   Rappresenta uno stato del suffix automaton.
 *   - len: lunghezza massima della sottostringa rappresentata
 *   - link: suffix link
 *   - first_pos: posizione di inizio di quella sottostringa (usato per LCS)
 *   - next: mappa fingerprint -> prossimo stato
 */
struct State {
    int len;
    int link;
    int first_pos;
    std::unordered_map<unsigned int, int> next;
};

/*
 * SuffixAutomaton:
 *   - st: vettore di stati
 *   - last: indice dello stato che rappresenta l'intera stringa
 */
struct SuffixAutomaton {
    std::vector<State> st;
    int last;
};

/*
 * build_suffix_automaton:
 *   Costruisce il suffix automaton per un vettore di fingerprint.
 */
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int>& A);

/*
 * match_suffix_automaton:
 *   Calcola la Longest Common Substring fra B e l'automa di A.
 *   Ritorna (lunghezza, startInA, startInB).
 */
std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B);

/*
 * longest_common_substring_suffix_automaton:
 *   Wrapper (non sempre utilizzato).
 */
std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int>& A,
    const std::vector<unsigned int>& B
);

#endif
