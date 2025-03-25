#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <string>
#include <tuple>
#include <vector>
using namespace std;

/*
 * graph_overlap_fp_precomputed:
 * Calcola la più lunga sottosequenza contigua comune tra due fingerprint compresse,
 * utilizzando rolling hash (double hashing) e ricerca binaria, con prefix hash precomputati.
 *
 * Input:
 *   - comp_fp1: fingerprint compressa della prima read.
 *   - prefix1_mod1, prefix1_mod2: prefix hash precomputati per la prima read.
 *   - idx1: vettore degli indici originali (comp_indices) della prima read.
 *   - comp_fp2, prefix2_mod1, prefix2_mod2, idx2: analoghi per la seconda read.
 *   - k: lunghezza del k‑mer (necessaria per ricostruire la regione nucleotidica).
 *
 * Restituisce una tuple:
 *   (overlap_len, orig_start1, orig_end1, orig_start2, orig_end2, comp_idx1, comp_idx2)
 * Se nessun match viene trovato, restituisce 0 e -1 per gli indici.
 */
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
