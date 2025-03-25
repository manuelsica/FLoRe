// Includiamo il nostro header overlap.hpp e altre librerie necessarie.
#include "overlap.hpp"
#include <vector>
#include <unordered_map>
#include <tuple>
#include <algorithm>
using namespace std;

// ---------------------------------------------------------------------------
// Funzione: compute_power
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola un vettore 'power' tale che power[i] = base^i mod mod, per i = 0,...,n.
//   Questi valori sono utilizzati per il calcolo dei prefix hash.
//   Esempio: se n = 3, base = 131, mod = 1000000007, allora power[0] = 1, power[1] = 131, power[2] = 17161, ...
// ---------------------------------------------------------------------------
static vector<long long> compute_power(int n, long long base, long long mod) {
    vector<long long> power(n+1, 1);
    power.reserve(n+1);
    for (int i = 1; i <= n; i++) {
        power[i] = (power[i-1] * base) % mod;
    }
    return power;
}

// ---------------------------------------------------------------------------
// Funzione: get_hash
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola l'hash della sottosequenza del vettore 'prefix' compresa tra gli indici l (incluso) e r (escluso).
//   Usa la formula: (prefix[r] - prefix[l] * power[r-l]) mod mod.
// ---------------------------------------------------------------------------
static inline long long get_hash(const vector<long long>& prefix, const vector<long long>& power, int l, int r, long long mod) {
    long long h = (prefix[r] - (prefix[l] * power[r - l]) % mod + mod) % mod;
    return h;
}

// ---------------------------------------------------------------------------
// Funzione: check_common_double
// ---------------------------------------------------------------------------
// Descrizione:
//   Dato un valore L, verifica se esiste una sottosequenza comune di lunghezza L
//   tra due fingerprint (fp1 e fp2) usando il double hashing.
//   Costruisce una mappa hash per tutte le sottosequenze di fp1 di lunghezza L, e poi
//   per ogni sottosequenza in fp2 verifica se esiste un hash corrispondente.
//   Se trova un match, esegue un controllo elemento per elemento per confermare.
//   Restituisce una tuple: (true, i, j) se viene trovato un match, dove
//   i è l'indice in fp1 e j in fp2; altrimenti (false, -1, -1).
// ---------------------------------------------------------------------------
static tuple<bool, int, int> check_common_double(
    const vector<long long>& fp1, const vector<long long>& fp2,
    const vector<long long>& prefix1_mod1, const vector<long long>& prefix1_mod2,
    const vector<long long>& prefix2_mod1, const vector<long long>& prefix2_mod2,
    const vector<long long>& power_mod1, const vector<long long>& power_mod2,
    int L, long long mod1, long long mod2)
{
    int n1 = (int)fp1.size(), n2 = (int)fp2.size();
    if (L == 0 || L > n1 || L > n2)
        return make_tuple(false, -1, -1);

    // Costruiamo la mappa hash dalla fingerprint più piccola, così da evitare
    // overhead eccessivo in caso di dataset squilibrati.
    bool firstIsSmaller = (n1 <= n2);

    const auto &fp_small      = firstIsSmaller ? fp1           : fp2;
    const auto &fp_big        = firstIsSmaller ? fp2           : fp1;
    const auto &prefix_small1 = firstIsSmaller ? prefix1_mod1  : prefix2_mod1;
    const auto &prefix_small2 = firstIsSmaller ? prefix1_mod2  : prefix2_mod2;
    const auto &prefix_big1   = firstIsSmaller ? prefix2_mod1  : prefix1_mod1;
    const auto &prefix_big2   = firstIsSmaller ? prefix2_mod2  : prefix1_mod2;

    int nSmall = (int)fp_small.size();
    int nBig   = (int)fp_big.size();

    // Prepariamo l'unordered_map con la giusta capacità
    int numSubseqSmall = nSmall - L + 1;
    unordered_map<unsigned long long, vector<int>> hash_map;
    hash_map.reserve(numSubseqSmall);

    // Costruiamo l'hash delle sottosequenze di lunghezza L nella fingerprint più piccola
    for (int i = 0; i <= nSmall - L; i++) {
        long long h1 = get_hash(prefix_small1, power_mod1, i, i + L, mod1);
        long long h2 = get_hash(prefix_small2, power_mod2, i, i + L, mod2);
        // Combina i due hash in un singolo valore (bitwise).
        unsigned long long combined = ((unsigned long long)h1 << 32) ^ ((unsigned long long)h2);
        hash_map[combined].push_back(i);
    }

    // Ricerchiamo match nelle sottosequenze della fingerprint più grande
    for (int j = 0; j <= nBig - L; j++) {
        long long h1 = get_hash(prefix_big1, power_mod1, j, j + L, mod1);
        long long h2 = get_hash(prefix_big2, power_mod2, j, j + L, mod2);
        unsigned long long combined = ((unsigned long long)h1 << 32) ^ ((unsigned long long)h2);

        auto it = hash_map.find(combined);
        if (it != hash_map.end()) {
            // Se esiste almeno una sottosequenza con lo stesso hash, verifichiamo
            // elemento per elemento per evitare collisioni.
            for (int iSmall : it->second) {
                bool match = true;
                for (int k = 0; k < L; k++) {
                    if (fp_small[iSmall + k] != fp_big[j + k]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    // Restituiamo gli indici nella convenzione originale: i in fp1, j in fp2
                    if (firstIsSmaller) {
                        return make_tuple(true, iSmall, j); // fp1 è la 'small'
                    } else {
                        return make_tuple(true, j, iSmall); // fp2 è la 'small'
                    }
                }
            }
        }
    }
    return make_tuple(false, -1, -1); // Nessun match trovato
}

// ---------------------------------------------------------------------------
// Funzione: graph_overlap_fp_precomputed
// ---------------------------------------------------------------------------
// Descrizione:
//   Utilizza una ricerca binaria per trovare la massima lunghezza L per cui esiste
//   una sottosequenza comune tra due fingerprint compressi (comp_fp1 e comp_fp2).
//   Utilizza la funzione check_common_double per verificare l'esistenza di una sottosequenza
//   comune di lunghezza L. Restituisce la lunghezza dell'overlap e, usando il vettore degli
//   indici (idx), converte gli indici della fingerprint compressa negli indici originali.
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
    int k)
{
    int n1 = (int)comp_fp1.size();
    int n2 = (int)comp_fp2.size();
    if (n1 == 0 || n2 == 0)
        return make_tuple(0, -1, -1, -1, -1, -1, -1);

    long long mod1 = 1000000007LL;
    long long mod2 = 1000000009LL;
    long long base = 131LL;

    int maxLen = (n1 < n2) ? n1 : n2;
    vector<long long> power_mod1 = compute_power(maxLen, base, mod1);
    vector<long long> power_mod2 = compute_power(maxLen, base, mod2);

    int low = 0, high = maxLen + 1;
    int bestL = 0;
    int best_i = -1, best_j = -1;

    // Ricerca binaria per la massima lunghezza L che produce un match.
    while (low < high - 1) {
        int mid = (low + high) / 2;
        auto [found, iCandidate, jCandidate] = check_common_double(
            comp_fp1, comp_fp2,
            prefix1_mod1, prefix1_mod2,
            prefix2_mod1, prefix2_mod2,
            power_mod1, power_mod2,
            mid, mod1, mod2
        );
        if (found) {
            low = mid;
            bestL = mid;
            best_i = iCandidate;
            best_j = jCandidate;
        } else {
            high = mid;
        }
    }

    if (bestL == 0)
        return make_tuple(0, -1, -1, -1, -1, -1, -1);

    // Converte gli indici compressi in indici originali usando il vettore idx.
    int orig_start1 = idx1[best_i];
    int orig_end1   = idx1[best_i + bestL - 1] + k;
    int orig_start2 = idx2[best_j];
    int orig_end2   = idx2[best_j + bestL - 1] + k;

    return make_tuple(bestL, orig_start1, orig_end1, orig_start2, orig_end2, best_i, best_j);
}
