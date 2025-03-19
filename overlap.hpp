#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <string>
#include <tuple>
#include <vector>
using namespace std;

/*
 * graph_overlap_fp:
 * Calcola la più lunga sottosequenza contigua comune tra due fingerprint compresse,
 * utilizzando rolling hash a doppia verifica e ricerca binaria.
 *
 * Input:
 *   - fp1, fp2: vettori dei codici della fingerprint compressa (long long)
 *   - idx1, idx2: vettori degli indici originali corrispondenti
 *   - k: lunghezza del k‑mer (necessaria per ricostruire la regione nucleotidica)
 *
 * Restituisce una tuple:
 *   (overlap_len, start1, end1, start2, end2)
 * Se nessun match viene trovato, restituisce 0 e -1 per gli indici.
 */
tuple<int, int, int, int, int> graph_overlap_fp(const vector<long long> &fp1,
                                                  const vector<long long> &fp2,
                                                  const vector<int> &idx1,
                                                  const vector<int> &idx2,
                                                  int k);

#endif
