// overlap.cpp
// Implementa le funzioni per il calcolo degli overlap tra le read.
// Include quattro metodi principali (FGOE, AOE, KHS, COIN) e il fallback al Suffix Automaton.
// Il file include commenti dettagliati per ogni blocco di codice.

#include "overlap.hpp"
#include "index.hpp"
#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <cmath>
#include <sstream>  // Necessario per std::ostringstream

/*
 * Metodo 1: Fingerprint-Guided Overlap Extension (FGOE)
 * Confronta il fingerprint compresso della prima read partendo dalla fine e quello della seconda partendo dall'inizio.
 * Restituisce una tupla (lunghezza_match, startA, startB).
 */
static std::tuple<int,int,int> fingerprint_guided_overlap_extension(const ProcessedRead &pr1,
                                                                     const ProcessedRead &pr2,
                                                                     int k,
                                                                     int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;  // Fingerprint compressi della prima read
    const auto &fp2 = pr2.comp.comp_fp;  // Fingerprint compressi della seconda read
    int size1 = (int)fp1.size();
    int size2 = (int)fp2.size();
    if (size1 == 0 || size2 == 0)
        return {0,0,0};
    int i = size1 - 1;
    int j = 0;
    int length_match = 0;
    // Confronta dal fondo di fp1 e dall'inizio di fp2
    while (i >= 0 && j < size2 && fp1[i] == fp2[j])
    {
        length_match++;
        i--;
        j++;
    }
    // Se la lunghezza match è sufficiente, restituisce gli indici di partenza
    if (length_match >= min_overlap)
    {
        int startA = pr1.comp.comp_indices[i+1];
        int startB = pr2.comp.comp_indices[0];
        return {length_match, startA, startB};
    }
    return {0,0,0};
}

/*
 * Metodo 2: Adaptive Overlap Extension (AOE)
 * Scorre tutte le possibili posizioni nei fingerprint per trovare il miglior match.
 */
static std::tuple<int,int,int> adaptive_overlap_extension(const ProcessedRead &pr1,
                                                          const ProcessedRead &pr2,
                                                          int k,
                                                          int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int n1 = (int)fp1.size();
    int n2 = (int)fp2.size();
    if (n1 == 0 || n2 == 0)
        return {0,0,0};
    int best_len = 0, best_x = 0, best_y = 0;
    // Scorre tutte le combinazioni di posizioni
    for (int x = 0; x < n1; x++)
    {
        for (int y = 0; y < n2; y++)
        {
            if (fp1[x] == fp2[y])
            {
                int tmp_len = 1;
                int xx = x + 1;
                int yy = y + 1;
                // Estende il match finché possibile
                while (xx < n1 && yy < n2 && fp1[xx] == fp2[yy])
                {
                    tmp_len++;
                    xx++;
                    yy++;
                }
                if (tmp_len > best_len)
                {
                    best_len = tmp_len;
                    best_x = x;
                    best_y = y;
                    if (best_len >= min_overlap)
                        break;
                }
            }
        }
        if (best_len >= min_overlap)
            break;
    }
    if (best_len >= min_overlap)
    {
        int startA = pr1.comp.comp_indices[best_x];
        int startB = pr2.comp.comp_indices[best_y];
        return {best_len, startA, startB};
    }
    return {0,0,0};
}

/*
 * Metodo 3: Kmer Hopping Search (KHS)
 * Costruisce una mappa delle posizioni per i k-mer della seconda read e cerca i match in fp1.
 */
static std::tuple<int,int,int> kmer_hopping_search(const ProcessedRead &pr1,
                                                   const ProcessedRead &pr2,
                                                   int k,
                                                   int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int n1 = (int)fp1.size();
    int n2 = (int)fp2.size();
    if (n1 == 0 || n2 == 0)
        return {0,0,0};
    // Costruisce una mappa: fingerprint -> vettore di posizioni in fp2
    std::unordered_map<unsigned int, std::vector<int>> positionsInB;
    positionsInB.reserve(n2);
    for (int i = 0; i < n2; i++)
    {
        positionsInB[fp2[i]].push_back(i);
    }
    int best_len = 0, best_x = 0, best_y = 0;
    for (int i = 0; i < n1; i++)
    {
        auto it = positionsInB.find(fp1[i]);
        if (it != positionsInB.end())
        {
            // Per ogni occorrenza in fp2, estende il match
            for (int posB : it->second)
            {
                int length_match = 1;
                int xx = i + 1;
                int yy = posB + 1;
                while (xx < n1 && yy < n2 && fp1[xx] == fp2[yy])
                {
                    length_match++;
                    xx++;
                    yy++;
                }
                if (length_match > best_len)
                {
                    best_len = length_match;
                    best_x = i;
                    best_y = posB;
                    if (best_len >= min_overlap)
                        break;
                }
            }
        }
        if (best_len >= min_overlap)
            break;
    }
    if (best_len >= min_overlap)
    {
        int startA = pr1.comp.comp_indices[best_x];
        int startB = pr2.comp.comp_indices[best_y];
        return {best_len, startA, startB};
    }
    return {0,0,0};
}

/*
 * Metodo 4: Combined Overlap IntelliSearch (COIN)
 * Combina i risultati parziali ottenuti dai metodi FGOE, AOE e KHS.
 * Se i match parziali sono vicini (entro una soglia k), li unisce.
 */
static std::tuple<int,int,int> combined_overlap_intellisearch(const ProcessedRead &pr1,
                                                              const ProcessedRead &pr2,
                                                              int k,
                                                              int min_overlap)
{
    // Calcola i match parziali con soglia ridotta (1)
    auto [fgoe_len, fgoe_startA, fgoe_startB] = fingerprint_guided_overlap_extension(pr1, pr2, k, 1);
    auto [aoe_len, aoe_startA, aoe_startB] = adaptive_overlap_extension(pr1, pr2, k, 1);
    auto [khs_len, khs_startA, khs_startB] = kmer_hopping_search(pr1, pr2, k, 1);

    if (fgoe_len == 0 && aoe_len == 0 && khs_len == 0)
        return {0,0,0};

    // Raggruppa i risultati in un vettore di tuple (lunghezza, startA, startB)
    std::vector<std::tuple<int,int,int>> partials = {
        {fgoe_len, fgoe_startA, fgoe_startB},
        {aoe_len, aoe_startA, aoe_startB},
        {khs_len, khs_startA, khs_startB}
    };
    // Ordina le tuple in base allo startA (per trovare match vicini)
    std::sort(partials.begin(), partials.end(),
              [](auto &a, auto &b) {
                  return std::get<1>(a) < std::get<1>(b);
              });

    int best_len = 0, best_startA = 0, best_startB = 0;
    // Combina ogni coppia di match parziali che siano vicini
    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 3; j++)
        {
            auto [len1, sA1, sB1] = partials[i];
            auto [len2, sA2, sB2] = partials[j];
            if (len1 > 0 && len2 > 0)
            {
                int endA1 = sA1 + len1;
                int endB1 = sB1 + len1;
                // Se la distanza tra il termine di un match e l'inizio dell'altro è entro k
                if (std::abs(sA2 - endA1) <= k && std::abs(sB2 - endB1) <= k)
                {
                    int combined_len = len1 + len2;
                    if (combined_len > best_len)
                    {
                        best_len = combined_len;
                        best_startA = std::min(sA1, sA2);
                        best_startB = std::min(sB1, sB2);
                    }
                }
            }
        }
    }
    // Fallback: prendi il match più lungo tra i parziali se non è possibile combinarli
    int naive_best = 0, naive_A = 0, naive_B = 0;
    for (auto &p : partials)
    {
        if (std::get<0>(p) > naive_best)
        {
            naive_best = std::get<0>(p);
            naive_A = std::get<1>(p);
            naive_B = std::get<2>(p);
        }
    }
    if (best_len >= min_overlap && best_len > naive_best)
        return {best_len, best_startA, best_startB};
    else if (naive_best >= min_overlap)
        return {naive_best, naive_A, naive_B};
    return {0,0,0};
}

/*
 * Funzione build_suffix_automaton:
 * Costruisce il Suffix Automaton a partire dal vettore di fingerprint A.
 * Implementa l'algoritmo classico per la costruzione del suffix automaton in O(n).
 */
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A)
{
    SuffixAutomaton sa;
    sa.st.reserve(2 * A.size()); // Prealloca spazio per almeno 2*n stati
    sa.last = 0;

    // Stato iniziale
    State initial;
    initial.len = 0;
    initial.link = -1;
    initial.first_pos = 0;
    sa.st.push_back(initial);

    // Costruzione dello automa per ogni carattere (k-mer codificato)
    for (int i = 0; i < (int)A.size(); i++)
    {
        unsigned int c = A[i];
        int cur = (int)sa.st.size();  // Nuovo stato da creare
        State stCur;
        stCur.len = sa.st[sa.last].len + 1;
        stCur.link = 0;  // Link temporaneo
        stCur.first_pos = i;
        sa.st.push_back(stCur);

        int p = sa.last;
        // Aggiorna le transizioni per tutti gli stati che non hanno la transizione 'c'
        while (p != -1 && sa.st[p].next.find(c) == sa.st[p].next.end())
        {
            sa.st[p].next[c] = cur;
            p = sa.st[p].link;
        }
        if (p == -1)
        {
            // Se non esiste alcun p valido, il link del nuovo stato è lo stato iniziale
            sa.st[cur].link = 0;
        }
        else
        {
            int q = sa.st[p].next[c];
            if (sa.st[p].len + 1 == sa.st[q].len)
            {
                sa.st[cur].link = q;
            }
            else
            {
                // Crea un clone dello stato q
                int clone = (int)sa.st.size();
                State stClone;
                stClone.len = sa.st[p].len + 1;
                stClone.next = sa.st[q].next;  // Copia le transizioni
                stClone.link = sa.st[q].link;
                stClone.first_pos = sa.st[q].first_pos;
                sa.st.push_back(stClone);
                // Aggiorna i link per gli stati che puntano a q
                while (p != -1 && sa.st[p].next[c] == q)
                {
                    sa.st[p].next[c] = clone;
                    p = sa.st[p].link;
                }
                sa.st[q].link = clone;
                sa.st[cur].link = clone;
            }
        }
        sa.last = cur;
    }
    return sa;
}

/*
 * Funzione match_suffix_automaton:
 * Dato il Suffix Automaton e il vettore B, trova il longest common substring (LCS).
 * Restituisce una tupla (lunghezza LCS, indice di partenza in A, indice di partenza in B).
 */
std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B)
{
    int v = 0;          // Stato corrente dell'automa
    int l = 0;          // Lunghezza corrente del match
    int best = 0;       // Migliore lunghezza trovata
    int bestposB = 0;   // Posizione in B in cui il match migliore termina
    int best_state = 0; // Stato corrispondente al match migliore

    for (int i = 0; i < (int)B.size(); i++)
    {
        unsigned int c = B[i];
        // Se esiste una transizione con il carattere c, prosegue il match
        if (sa.st[v].next.find(c) != sa.st[v].next.end())
        {
            v = sa.st[v].next.at(c);
            l++;
        }
        else
        {
            // Altrimenti, retrocede seguendo i link suffix fino a trovare una transizione
            while (v != -1 && sa.st[v].next.find(c) == sa.st[v].next.end())
            {
                v = sa.st[v].link;
            }
            if (v == -1)
            {
                v = 0;
                l = 0;
            }
            else
            {
                l = sa.st[v].len + 1;
                v = sa.st[v].next.at(c);
            }
        }
        if (l > best)
        {
            best = l;
            bestposB = i;
            best_state = v;
        }
    }
    int startB = bestposB - best + 1;
    int startA = sa.st[best_state].first_pos - best + 1;
    return {best, startA, startB};
}

/*
 * Funzione longest_common_substring_suffix_automaton:
 * Costruisce il suffix automaton per A e usa match_suffix_automaton per trovare l'LCS con B.
 */
std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int> &A,
    const std::vector<unsigned int> &B)
{
    SuffixAutomaton sa = build_suffix_automaton(A);
    return match_suffix_automaton(sa, B);
}

/*
 * Funzione compare_candidate_pair:
 * Per due ReadData, esegue i 4 metodi (FGOE, AOE, KHS, COIN) per determinare l'overlap migliore.
 * Se nessuno dei metodi raggiunge la lunghezza minima, prova a usare il Suffix Automaton.
 */
OverlapResult compare_candidate_pair(ReadData &r1,
                                     ReadData &r2,
                                     int k,
                                     int min_overlap,
                                     bool /*verbose*/,
                                     int /*max_repeat_threshold*/)
{
    OverlapResult best;  // Risultato inizializzato a zero

    // Funzione lambda per aggiornare "best" se si trova un overlap migliore.
    auto try_update = [&](int len, int comp_idx1, int comp_idx2,
                          const ProcessedRead &prA, const ProcessedRead &prB,
                          const std::string &readA, const std::string &readB,
                          const std::string &orientA, const std::string &orientB,
                          const std::string &comb, const std::string &algorithmName)
    {
        if (len > best.overlap_len)
        {
            best.overlap_len = len;
            best.combination = comb;
            best.orientation1 = orientA;
            best.orientation2 = orientB;
            best.start1 = prA.comp.comp_indices[comp_idx1];
            best.end1 = prA.comp.comp_indices[comp_idx1 + len - 1] + k;
            best.start2 = prB.comp.comp_indices[comp_idx2];
            best.end2 = prB.comp.comp_indices[comp_idx2 + len - 1] + k;
            best.r1 = readA;
            best.r2 = readB;
            // Utilizza una lambda per costruire la rappresentazione testuale del fingerprint
            auto getFingerprintRegion = [&](const ProcessedRead &p, int cstart, int length){
                std::ostringstream oss;
                int stop = std::min<int>(cstart + length, (int)p.comp.comp_fp.size());
                for (int i = cstart; i < stop; i++)
                {
                    oss << p.comp.comp_fp[i];
                    if (i < stop - 1)
                        oss << "-";
                }
                return oss.str();
            };
            best.fingerprint_r1 = getFingerprintRegion(prA, comp_idx1, len);
            best.fingerprint_r2 = getFingerprintRegion(prB, comp_idx2, len);
            best.used_algorithm = algorithmName;
        }
    };

    // Lambda per eseguire i 4 metodi su un particolare orientamento
    auto compute_best_overlap_4stages = [&](ProcessedRead &A, ProcessedRead &B,
                                             std::string &rA, std::string &rB,
                                             const std::string &oA, const std::string &oB,
                                             const std::string &comb)
    {
        // Metodo 1: FGOE
        {
            auto [matchLen, startA, startB] = fingerprint_guided_overlap_extension(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                // Trova l'indice del fingerprint corrispondente al valore startA in A
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++)
                {
                    if (A.comp.comp_indices[i] == startA)
                    {
                        idxA = i;
                        break;
                    }
                }
                // Trova l'indice del fingerprint corrispondente al valore startB in B
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++)
                {
                    if (B.comp.comp_indices[j] == startB)
                    {
                        idxB = j;
                        break;
                    }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "FGOE");
                return;
            }
        }
        // Metodo 2: AOE
        {
            auto [matchLen, startA, startB] = adaptive_overlap_extension(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++)
                {
                    if (A.comp.comp_indices[i] == startA)
                    {
                        idxA = i;
                        break;
                    }
                }
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++)
                {
                    if (B.comp.comp_indices[j] == startB)
                    {
                        idxB = j;
                        break;
                    }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "AOE");
                return;
            }
        }
        // Metodo 3: KHS
        {
            auto [matchLen, startA, startB] = kmer_hopping_search(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++)
                {
                    if (A.comp.comp_indices[i] == startA)
                    {
                        idxA = i;
                        break;
                    }
                }
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++)
                {
                    if (B.comp.comp_indices[j] == startB)
                    {
                        idxB = j;
                        break;
                    }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "KHS");
                return;
            }
        }
        // Metodo 4: COIN
        {
            auto [matchLen, startA, startB] = combined_overlap_intellisearch(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++)
                {
                    if (A.comp.comp_indices[i] == startA)
                    {
                        idxA = i;
                        break;
                    }
                }
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++)
                {
                    if (B.comp.comp_indices[j] == startB)
                    {
                        idxB = j;
                        break;
                    }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "COIN");
                return;
            }
        }
    };

    // Lambda per provare il Suffix Automaton come fallback se nessuno dei 4 metodi raggiunge min_overlap.
    auto compute_best_overlap_suffix_automaton = [&](ProcessedRead &A, ProcessedRead &B,
                                                      std::string &rA, std::string &rB,
                                                      const std::string &oA, const std::string &oB,
                                                      const std::string &comb,
                                                      bool &saBuilt,
                                                      SuffixAutomaton &saCache)
    {
        if (!saBuilt)
        {
            saCache = build_suffix_automaton(A.comp.comp_fp);
            saBuilt = true;
        }
        auto [length, c1, c2] = match_suffix_automaton(saCache, B.comp.comp_fp);
        if (length > best.overlap_len)
        {
            best.overlap_len = length;
            best.combination = comb;
            best.orientation1 = oA;
            best.orientation2 = oB;
            best.start1 = A.comp.comp_indices[c1];
            best.end1 = A.comp.comp_indices[c1 + length - 1] + k;
            best.start2 = B.comp.comp_indices[c2];
            best.end2 = B.comp.comp_indices[c2 + length - 1] + k;
            best.r1 = rA;
            best.r2 = rB;
            auto getFingerprintRegion = [&](const ProcessedRead &p, int cstart, int ln){
                std::ostringstream oss;
                int stop = std::min<int>(cstart + ln, (int)p.comp.comp_fp.size());
                for (int i = cstart; i < stop; i++)
                {
                    oss << p.comp.comp_fp[i];
                    if (i < stop - 1)
                        oss << "-";
                }
                return oss.str();
            };
            best.fingerprint_r1 = getFingerprintRegion(A, c1, length);
            best.fingerprint_r2 = getFingerprintRegion(B, c2, length);
            best.used_algorithm = "SuffixAutomaton";
        }
    };

    // Esegue i 4 metodi su tutte le possibili combinazioni di orientamento:
    // ff: forward-forward, fr: forward-reverse, rf: reverse-forward, rr: reverse-reverse.
    compute_best_overlap_4stages(r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward, "forward", "forward", "ff");
    compute_best_overlap_4stages(r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq, "forward", "reverse", "fr");
    compute_best_overlap_4stages(r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward, "reverse", "forward", "rf");
    compute_best_overlap_4stages(r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq, "reverse", "reverse", "rr");

    // Se non si raggiunge min_overlap, prova il fallback con il Suffix Automaton
    if (best.overlap_len < min_overlap)
    {
        compute_best_overlap_suffix_automaton(r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward, "forward", "forward", "ff", r1.sa_built_fwd, r1.suffix_automaton_fwd);
        compute_best_overlap_suffix_automaton(r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq, "forward", "reverse", "fr", r1.sa_built_fwd, r1.suffix_automaton_fwd);
        compute_best_overlap_suffix_automaton(r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward, "reverse", "forward", "rf", r1.sa_built_rev, r1.suffix_automaton_rev);
        compute_best_overlap_suffix_automaton(r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq, "reverse", "reverse", "rr", r1.sa_built_rev, r1.suffix_automaton_rev);
    }

    return best;
}