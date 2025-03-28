#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <string>
#include <tuple>
#include <vector>
using namespace std;

/*
 * longest_common_substring_suffix_automaton:
 * Calcola il longest common substring tra due vettori di unsigned int (le fingerprint compresse)
 * utilizzando un suffix automaton in tempo lineare.
 * Restituisce una tuple: (lunghezza_massima, indice_in_A, indice_in_B)
 */
tuple<int, int, int> longest_common_substring_suffix_automaton(const vector<unsigned int>& A, const vector<unsigned int>& B);

#endif
