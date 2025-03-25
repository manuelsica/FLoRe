#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <string>
#include <tuple>
#include <vector>
using namespace std;

// ---------------------------------------------------------------------------
// Funzione: graph_overlap_fp_precomputed
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola l'overlap massimo (numero di k‑mer consecutivi corrispondenti) tra due
//   fingerprint compresse (comp_fp1 e comp_fp2) utilizzando double hashing e
//   ricerca binaria. Per ogni lunghezza L, usa la funzione check_common_double per verificare
//   se esiste una sottosequenza comune. Se sì, aumenta L; altrimenti riduce L.
//   Gli hash dei prefix sono già precomputati (prefix1_mod1, prefix1_mod2, ecc.).
//
// Input:
//   - comp_fp1: Vettore dei codici compressi della prima read.
//   - prefix1_mod1, prefix1_mod2: Vettori dei prefix hash per la prima read.
//   - idx1: Vettore degli indici originali per la prima read (per convertire gli indici compressi).
//   - comp_fp2, prefix2_mod1, prefix2_mod2, idx2: Parametri analoghi per la seconda read.
//   - k: Lunghezza del k‑mer (necessaria per calcolare l'indice finale nella read originale).
//
// Output:
//   Restituisce una tuple:
//     (overlap_len, orig_start1, orig_end1, orig_start2, orig_end2, comp_idx1, comp_idx2)
//   Dove:
//     - overlap_len: Numero di k‑mer consecutivi corrispondenti (match fingerprint).
//     - orig_start1: Indice iniziale dell'overlap nella prima read (ottenuto da idx1).
//     - orig_end1: Indice finale (non incluso) nella prima read, calcolato come idx1[comp_idx1 + overlap_len - 1] + k.
//     - Gli altri parametri sono analoghi per la seconda read.
// ---------------------------------------------------------------------------
tuple<int, int, int, int, int, int, int> graph_overlap_fp_precomputed(
    const vector<long long>& comp_fp1,
    const vector<long long>& prefix1_mod1,
    const vector<long long>& prefix1_mod2,
    const vector<int>& idx1,
    const vector<long long>& comp_fp2,
    const vector<long long>& prefix2_mod1,
    const vector<long long>& prefix2_mod2,
    const vector<int>& idx2,
    int k
);

#endif
