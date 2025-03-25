#ifndef UTIL_HPP
#define UTIL_HPP

// Includiamo le librerie standard per gestire stringhe e vettori.
#include <string>
#include <vector>
using namespace std;

// ---------------------------------------------------------------------------
// Funzione: compute_val
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola un valore numerico per il k‑mer rappresentato dalla stringa 'factor'.
//   Usa il mapping: A → 0, C → 1, G → 2, T → 3, e interpreta il k‑mer come numero
//   in una base (default = 5). Ad esempio, "AC" diventa 0*5 + 1 = 1.
// ---------------------------------------------------------------------------
long long compute_val(const string &factor, int base = 5);

// ---------------------------------------------------------------------------
// Funzione: encode_factor
// ---------------------------------------------------------------------------
// Descrizione:
//   Codifica il k‑mer combinando il valore numerico calcolato con compute_val e
//   la lunghezza del k‑mer tramite un'operazione di shift.
//   Restituisce: (compute_val(factor) << shift) | len(factor)
//   Esempio: se compute_val("AC") = 1 e len("AC") = 2, con shift = 8, si ottiene: (1 << 8) | 2 = 258.
// ---------------------------------------------------------------------------
long long encode_factor(const string &factor, int shift = 8, int base = 5);

// ---------------------------------------------------------------------------
// Funzione: sliding_window_read
// ---------------------------------------------------------------------------
// Descrizione:
//   Segmenta la read in k‑mer usando una finestra scorrevole (step = 1).
//   Se read = "ATCGT" e k = 3, ritorna {"ATC", "TCG", "CGT"}.
// ---------------------------------------------------------------------------
vector<string> sliding_window_read(const string &read, int k, int step = 1);

// ---------------------------------------------------------------------------
// Funzione: reverse_complement
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola il reverse complement di una sequenza nucleotidica.
//   Ad esempio, per "ATCG" restituisce "CGAT".
// ---------------------------------------------------------------------------
string reverse_complement(const string &seq);

// ---------------------------------------------------------------------------
// Struttura: CompressedFingerprint
// ---------------------------------------------------------------------------
// Descrizione:
//   Memorizza la fingerprint compressa di una read. Include:
//     - comp_fp: vettore dei codici compressi (senza duplicati consecutivi)
//     - comp_kmers: i k‑mer corrispondenti (la prima occorrenza di ogni codice)
//     - comp_indices: gli indici originali nella read (riferiti ai k‑mer)
// ---------------------------------------------------------------------------
struct CompressedFingerprint {
    vector<long long> comp_fp;
    vector<string> comp_kmers;
    vector<int> comp_indices;
};

// ---------------------------------------------------------------------------
// Funzione: compress_fingerprint
// ---------------------------------------------------------------------------
// Descrizione:
//   Comprimi la fingerprint eliminando duplicati consecutivi. Ad esempio, se la
//   fingerprint è {10, 10, 12, 12, 12, 14}, restituisce {10, 12, 14} con gli indici
//   corrispondenti.
// ---------------------------------------------------------------------------
CompressedFingerprint compress_fingerprint(const vector<long long> &fingerprint, const vector<string> &kmers);

// ---------------------------------------------------------------------------
// Funzione Template: duval_factorization
// ---------------------------------------------------------------------------
// Descrizione:
//   Applica l'algoritmo di Duval per la Lyndon Factorization su un vettore 's'.
//   Restituisce un vettore di vettori, dove ogni vettore rappresenta un fattore Lyndon.
//   Esempio:
//     Se s = {3, 2, 2, 3, 1} (dove gli elementi sono confrontabili),
//     la funzione dividerà s in fattori che sono le più piccole stringhe in ordine lessicografico.
// ---------------------------------------------------------------------------
template<typename T>
vector<vector<T>> duval_factorization(const vector<T>& s) {
    vector<vector<T>> factors;    // Vettore dei fattori Lyndon
    int n = (int)s.size();        // Lunghezza del vettore in input
    int i = 0;
    // Mentre non abbiamo elaborato tutta la sequenza...
    while (i < n) {
        int j = i + 1, k = i;
        // Cerca il punto in cui la sequenza smette di essere in ordine non decrescente.
        while (j < n && s[k] <= s[j]) {
            if (s[k] < s[j])
                k = i;        // Se troviamo un elemento maggiore, resettiamo k all'inizio
            else
                k++;          // Altrimenti incrementiamo k
            j++;
        }
        // Aggiunge ripetutamente il fattore Lyndon ottenuto
        while (i <= k) {
            factors.push_back(vector<T>(s.begin() + i, s.begin() + i + (j - k)));
            i += (j - k);
        }
    }
    return factors;
}

// ---------------------------------------------------------------------------
// Struttura: ProcessedRead
// ---------------------------------------------------------------------------
// Descrizione:
//   Memorizza il risultato del preprocessing di una read, che include:
//     - comp: la fingerprint compressa (con i suoi vettori associati)
//     - kmers: la lista completa dei k‑mer
//     - fingerprint: il vettore completo di codici (non compresso)
//     - comp_prefix_mod1 e comp_prefix_mod2: vettori dei prefix hash per il double hashing
//     - lyndon_factors: la Lyndon Factorization (Duval) della fingerprint compressa,
//       utile per ottimizzare ulteriormente i confronti.
// ---------------------------------------------------------------------------
struct ProcessedRead {
    CompressedFingerprint comp;
    vector<string> kmers;
    vector<long long> fingerprint;
    vector<long long> comp_prefix_mod1;
    vector<long long> comp_prefix_mod2;
    vector<vector<long long>> lyndon_factors; // Fattori Lyndon calcolati con Duval
};

// ---------------------------------------------------------------------------
// Funzione: process_read
// ---------------------------------------------------------------------------
// Descrizione:
//   Data una read e un valore di k, questa funzione esegue:
//     1. Segmentazione della read in k‑mer (finestra scorrevole).
//     2. Codifica di ogni k‑mer in un intero.
//     3. Creazione della fingerprint (vettore dei codici) e la sua compressione.
//     4. Calcolo dei prefix hash per il double hashing.
//     5. Calcolo della Duval (Lyndon) Factorization della fingerprint compressa.
// ---------------------------------------------------------------------------
ProcessedRead process_read(const string &read, int k);

#endif
