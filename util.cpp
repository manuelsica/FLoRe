// Includiamo il nostro header util.hpp e altre librerie necessarie.
#include "util.hpp"
#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// Funzione: compute_prefix_hash
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola il vettore dei prefix hash per un vettore di valori 'arr'.
//   Ogni elemento prefix[i+1] rappresenta l'hash della sottosequenza arr[0..i],
//   calcolato con la formula: prefix[i+1] = (prefix[i] * base + arr[i]) mod mod.
//   Esempio: se arr = {a, b, c}, base = 131 e mod = 1000000007, allora:
//     prefix[0] = 0
//     prefix[1] = (0*131 + a) mod mod = a
//     prefix[2] = (a*131 + b) mod mod
//     prefix[3] = (prefix[2]*131 + c) mod mod
// ---------------------------------------------------------------------------
static vector<long long> compute_prefix_hash(const vector<long long>& arr, long long base, long long mod) {
    int n = (int)arr.size();
    vector<long long> prefix(n+1, 0); // Inizializza un vettore di lunghezza n+1 con 0
    prefix.reserve(n+1);
    for (int i = 0; i < n; i++) {
        prefix[i+1] = (prefix[i] * base + (arr[i] % mod)) % mod;
    }
    return prefix;
}

// ---------------------------------------------------------------------------
// Funzione: compute_val
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola un valore numerico per la stringa 'factor' basandosi sul mapping:
//   'A'->0, 'C'->1, 'G'->2, 'T'->3.
//   Esempio: per "AC" e base=5, compute_val("AC") = 0*5 + 1 = 1.
// ---------------------------------------------------------------------------
long long compute_val(const string &factor, int base) {
    int L = (int)factor.size();  // Lunghezza della stringa
    long long val = 0;           // Variabile in cui accumulare il valore

    // Calcolo iterativo: val = (val * base) + m
    // dove m è il valore numerico del carattere.
    for (int i = 0; i < L; i++) {
        int m;
        switch(factor[i]) {
            case 'A': m = 0; break; // Esempio: 'A' → 0
            case 'C': m = 1; break; // 'C' → 1
            case 'G': m = 2; break; // 'G' → 2
            case 'T': m = 3; break; // 'T' → 3
            default: m = 0; break;  // Qualsiasi altro carattere viene trattato come 0
        }
        val = val * base + m;
    }
    return val;
}

// ---------------------------------------------------------------------------
// Funzione: encode_factor
// ---------------------------------------------------------------------------
// Descrizione:
//   Combina il valore numerico del k‑mer (ottenuto da compute_val) con la sua lunghezza,
//   spostando a sinistra il valore e usando un OR bit a bit per includere la lunghezza.
//   Esempio: per "AC", se compute_val("AC") = 1 e la lunghezza è 2, con shift=8 otteniamo:
//   (1 << 8) | 2 = 256 | 2 = 258.
// ---------------------------------------------------------------------------
long long encode_factor(const string &factor, int shift, int base) {
    int L = (int)factor.size();
    long long val = compute_val(factor, base);
    return (val << shift) | L;
}

// ---------------------------------------------------------------------------
// Funzione: sliding_window_read
// ---------------------------------------------------------------------------
// Descrizione:
//   Segmenta la read in k‑mer utilizzando una finestra scorrevole.
//   Se read = "ATCGT" e k = 3, allora:
//     i = 0: "ATC"
//     i = 1: "TCG"
//     i = 2: "CGT"
//   Restituisce un vettore con questi k‑mer.
// ---------------------------------------------------------------------------
vector<string> sliding_window_read(const string &read, int k, int step) {
    vector<string> segments;
    if (k <= 0 || (int)read.size() < k) {
        return segments;
    }
    int n = (int)read.size() - k + 1; // Numero di k‑mer possibili
    if (n <= 0) return segments;
    
    // Prealloca memoria per efficienza
    int reserveSize = (step > 0) ? (1 + (n - 1) / step) : n;
    segments.reserve(reserveSize);
    
    // Ciclo per estrarre i k‑mer
    for (int i = 0; i + k <= (int)read.size(); i += step) {
        segments.push_back(read.substr(i, k));
    }
    return segments;
}

// ---------------------------------------------------------------------------
// Funzione: reverse_complement
// ---------------------------------------------------------------------------
// Descrizione:
//   Calcola il reverse complement di una sequenza nucleotidica.
//   Ad esempio: "ATCG" -> Reverse: "GCTA", Complement: "CGAT".
// ---------------------------------------------------------------------------
string reverse_complement(const string &seq) {
    string result;
    result.reserve(seq.size()); // Prealloca spazio per la sequenza complementare
    // Itera dalla fine all'inizio della sequenza
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        // Per ogni carattere, aggiunge il complementare
        switch(*it) {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default: result.push_back('N'); break; // Se il carattere non è A, C, G o T, usa 'N'
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// Funzione: compress_fingerprint
// ---------------------------------------------------------------------------
// Descrizione:
//   Data una fingerprint (vettore di codici) e il vettore dei corrispondenti k‑mer,
//   elimina i duplicati consecutivi e salva i k‑mer e gli indici corrispondenti.
//   Esempio: se fingerprint = {100, 100, 101, 101, 102} allora verrà restituito:
//     comp_fp = {100, 101, 102}
//     comp_kmers = { primo k-mer per 100, primo per 101, primo per 102 }
//     comp_indices = { indice del primo 100, indice del primo 101, indice del primo 102 }
// ---------------------------------------------------------------------------
CompressedFingerprint compress_fingerprint(const vector<long long> &fingerprint, const vector<string> &kmers) {
    CompressedFingerprint comp;
    if (fingerprint.empty())
        return comp; // Se il vettore è vuoto, restituisce una struttura vuota

    comp.comp_fp.reserve(fingerprint.size());
    comp.comp_kmers.reserve(kmers.size());
    comp.comp_indices.reserve(kmers.size());

    // Inserisce il primo elemento di fingerprint e il relativo k-mer e indice 0.
    comp.comp_fp.push_back(fingerprint[0]);
    comp.comp_kmers.push_back(kmers[0]);
    comp.comp_indices.push_back(0);

    // Scorre la fingerprint a partire dal secondo elemento.
    for (size_t i = 1; i < fingerprint.size(); i++) {
        if (fingerprint[i] != fingerprint[i - 1]) { // Se l'elemento corrente è diverso dal precedente...
            comp.comp_fp.push_back(fingerprint[i]);    // lo aggiunge al vettore compresso
            comp.comp_kmers.push_back(kmers[i]);       // salva il corrispondente k‑mer
            comp.comp_indices.push_back((int)i);       // salva l'indice originale
        }
    }
    return comp;
}

// ---------------------------------------------------------------------------
// Funzione: process_read
// ---------------------------------------------------------------------------
// Descrizione:
//   Preprocessa una read eseguendo:
//     1. Segmentazione in k‑mer (finestra scorrevole).
//     2. Codifica di ogni k‑mer in un intero (usando encode_factor).
//     3. Compressione della fingerprint per eliminare duplicati consecutivi.
//     4. Calcolo dei prefix hash sulla fingerprint compressa (usando double hashing).
//     5. Calcolo della Duval Factorization della fingerprint compressa, che permette
//        di ottenere una scomposizione in fattori Lyndon.
//   Restituisce una struttura ProcessedRead contenente tutti questi dati.
// ---------------------------------------------------------------------------
ProcessedRead process_read(const string &read, int k) {
    ProcessedRead pr;
    if (k <= 0 || (int)read.size() < k) {
        return pr;
    }

    int n = (int)read.size() - k + 1; // Calcola quanti k‑mer possono essere estratti
    pr.kmers.reserve(n);             // Prealloca memoria per il vettore dei k‑mer
    pr.fingerprint.reserve(n);       // Prealloca memoria per il vettore della fingerprint

    pr.kmers = sliding_window_read(read, k, 1); // Divide la read in k‑mer con step 1
    pr.fingerprint.resize(pr.kmers.size());

    for (size_t i = 0; i < pr.kmers.size(); i++) {
        // Codifica ogni k‑mer in un intero
        pr.fingerprint[i] = encode_factor(pr.kmers[i], 8, 5);
    }
    // Comprimi la fingerprint per eliminare duplicati consecutivi.
    pr.comp = compress_fingerprint(pr.fingerprint, pr.kmers);

    // Parametri per il double hashing: moduli e base.
    long long mod1 = 1000000007LL, mod2 = 1000000009LL, base = 131LL;

    // Calcola i prefix hash della fingerprint compressa.
    pr.comp_prefix_mod1 = compute_prefix_hash(pr.comp.comp_fp, base, mod1);
    pr.comp_prefix_mod2 = compute_prefix_hash(pr.comp.comp_fp, base, mod2);

    // Calcola la Duval Factorization della fingerprint compressa
    pr.lyndon_factors = duval_factorization(pr.comp.comp_fp);

    return pr;
}
