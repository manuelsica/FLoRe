// overlap.hpp
// Header per le funzioni di calcolo degli overlap tra le read, inclusi i metodi FGOE, AOE, KHS, COIN
// e l'utilizzo del Suffix Automaton per il recupero della Longest Common Substring (LCS).

#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <vector>
#include <unordered_map>
#include <tuple>
#include <string>
#include "util.hpp"  // Per ProcessedRead e CompressedFingerprint

/*
 * Struttura State:
 * Rappresenta uno stato del Suffix Automaton, con lunghezza, link, posizione e mappa delle transizioni.
 */
struct State {
    int len;                                  // Lunghezza massima della sottostringa rappresentata
    int link;                                 // Link suffix (stato di fallback)
    int first_pos;                            // Prima posizione in cui compare la sottostringa
    std::unordered_map<unsigned int, int> next; // Mappa delle transizioni: carattere -> indice dello stato
};

/*
 * Struttura SuffixAutomaton:
 * Contiene tutti gli stati e l'indice dell'ultimo stato creato.
 */
struct SuffixAutomaton {
    std::vector<State> st;  // Vettore degli stati
    int last;               // Indice dell'ultimo stato (stato corrente)
};

/*
 * Funzione build_suffix_automaton:
 * Costruisce il Suffix Automaton a partire da un vettore di fingerprint (unsigned int).
 */
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A);

/*
 * Funzione match_suffix_automaton:
 * Dato un Suffix Automaton e un vettore B, trova il longest common substring (LCS)
 * e restituisce una tupla contenente la lunghezza dell'LCS e gli indici di partenza in A e B.
 */
std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B);

/*
 * Funzione longest_common_substring_suffix_automaton:
 * Costruisce il Suffix Automaton per A e poi chiama match_suffix_automaton per B.
 */
std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int> &A,
    const std::vector<unsigned int> &B);

/*
 * Struttura OverlapResult:
 * Contiene tutti i dettagli dell'overlap trovato tra due read.
 */
struct OverlapResult
{
    int overlap_len = 0;              // Lunghezza dell'overlap
    int start1 = 0;                   // Indice di inizio dell'overlap nella prima read
    int end1 = 0;                     // Indice di fine dell'overlap nella prima read
    int start2 = 0;                   // Indice di inizio dell'overlap nella seconda read
    int end2 = 0;                     // Indice di fine dell'overlap nella seconda read
    std::string orientation1;         // Orientamento della prima read ("forward" o "reverse")
    std::string orientation2;         // Orientamento della seconda read
    std::string combination;          // Codice per la combinazione degli algoritmi usati
    std::string r1;                   // Sequenza della prima read
    std::string r2;                   // Sequenza della seconda read
    std::string fingerprint_r1;       // Rappresentazione testuale del fingerprint della prima read nell'overlap
    std::string fingerprint_r2;       // Rappresentazione testuale del fingerprint della seconda read
    std::string used_algorithm;       // Nome dell'algoritmo che ha prodotto l'overlap
};

/*
 * Funzione compare_candidate_pair:
 * Dati due ReadData, esegue diversi metodi (FGOE, AOE, KHS, COIN)
 * per trovare il miglior overlap, eventualmente integrando con il Suffix Automaton.
 */
OverlapResult compare_candidate_pair(struct ReadData &r1,
                                     struct ReadData &r2,
                                     int k,
                                     int min_overlap,
                                     bool verbose,
                                     int max_repeat_threshold);

#endif