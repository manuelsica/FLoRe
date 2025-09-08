// overlap.cpp
#include "overlap.hpp"
#include "index.hpp"
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <string>
#include <sstream>
#include <cassert>

// Cerca “x” esattamente in idx; ritorna l’indice o -1 se non trovato
static int find_comp_index(const std::vector<int> &idx, int x) {
    auto it = std::lower_bound(idx.begin(), idx.end(), x);
    if (it != idx.end() && *it == x)
        return int(std::distance(idx.begin(), it));
    return -1;
}

/* =========================  METODO 1: FGOE  =========================
 * Offset-voting + guided extension su fingerprint compressi.
 * 1) Indicizza posizioni di B per valore di fingerprint.
 * 2) Per ogni match (a,b) con fa[a]==fb[b], vota la diagonale d=a-b.
 * 3) Seleziona d* con il massimo supporto; misura la run contigua più lunga
 *    su d* e, se >= min_overlap, restituisce (len, startA, startB).
 *
 * Ritorna: (len, startA, startB) con start* in coordinate di comp_indices.
 */
static std::tuple<int,int,int>
fingerprint_guided_overlap_extension(const ProcessedRead &A,
                                     const ProcessedRead &B,
                                     int /*k*/, int min_overlap)
{
    const auto &fa = A.comp.comp_fp;
    const auto &fb = B.comp.comp_fp;
    const auto &ia = A.comp.comp_indices;
    const auto &ib = B.comp.comp_indices;

    const int na = (int)fa.size();
    const int nb = (int)fb.size();
    if (na == 0 || nb == 0) return {0,0,0};

    // 1) posizioni in B per fingerprint
    std::unordered_map<unsigned int, std::vector<int>> posB;
    posB.reserve((size_t)nb);
    for (int j = 0; j < nb; ++j) posB[fb[j]].push_back(j);

    // 2) votazione diagonali d = a - b
    std::unordered_map<int,int> votes;
    votes.reserve((size_t)std::min(na, nb));
    int d_star = 0, best_votes = 0;

    for (int a = 0; a < na; ++a) {
        auto it = posB.find(fa[a]);
        if (it == posB.end()) continue;
        const auto &vec = it->second; // j in ordine crescente
        for (int b : vec) {
            const int d = a - b;
            int v = ++votes[d];
            if (v > best_votes) {
                best_votes = v;
                d_star = d;
            }
        }
    }

    // Se la diagonale migliore non ha abbastanza supporto, termina
    if (best_votes < min_overlap) return {0,0,0};

    // 3) Estrazione delle coppie sulla diagonale d* e misura della run contigua
    std::vector<std::pair<int,int>> pairs; // (a,b) con a-b=d_star
    pairs.reserve(best_votes);

    for (int a = 0; a < na; ++a) {
        auto it = posB.find(fa[a]);
        if (it == posB.end()) continue;
        for (int b : it->second) {
            if (a - b == d_star) pairs.emplace_back(a, b);
        }
    }
    if (pairs.empty()) return {0,0,0};

    // Le coppie sono già in ordine non decrescente per 'a' (a cresce, b nelle liste è crescente),
    // ma per sicurezza ordiniamo per 'a' e poi per 'b'.
    std::sort(pairs.begin(), pairs.end());

    // Scorriamo per trovare il blocco contiguo massimo (a e b crescono entrambi di 1)
    int best_len = 1, best_a0 = pairs[0].first, best_b0 = pairs[0].second;
    int cur_len = 1, cur_a0 = pairs[0].first, cur_b0 = pairs[0].second;

    for (size_t t = 1; t < pairs.size(); ++t) {
        const auto [a, b] = pairs[t];
        const auto [pa, pb] = pairs[t-1];

        if (a == pa + 1 && b == pb + 1) {
            // continuiamo la run contigua
            ++cur_len;
        } else {
            // chiudiamo la run
            if (cur_len > best_len) {
                best_len = cur_len;
                best_a0 = cur_a0;
                best_b0 = cur_b0;
            }
            cur_len = 1;
            cur_a0 = a;
            cur_b0 = b;
        }
    }
    // ultima run
    if (cur_len > best_len) {
        best_len = cur_len;
        best_a0 = cur_a0;
        best_b0 = cur_b0;
    }

    if (best_len >= min_overlap) {
        // converti in coordinate di comp_indices
        int startA = ia[best_a0];
        int startB = ib[best_b0];
        return {best_len, startA, startB};
    }
    return {0,0,0};
}
/* =========================  METODO 2: AOE  =========================
 * Cerca un'estensione adattiva: trova un seed fa[i] in fb e misura il
 * tratto contiguo uguale da lì in avanti (conta la lunghezza len).
 * Ritorna la miglior tripletta (len, startA, startB).
 */
static std::tuple<int,int,int>
adaptive_overlap_extension(const ProcessedRead &A,
                           const ProcessedRead &B,
                           int /*k*/, int min_overlap)
{
    const auto &fa = A.comp.comp_fp;
    const auto &fb = B.comp.comp_fp;
    const auto &ia = A.comp.comp_indices;
    const auto &ib = B.comp.comp_indices;

    int best_len = 0, bestA = 0, bestB = 0;

    for (int i = 0; i < (int)fa.size(); ++i) {
        auto it = std::find(fb.begin(), fb.end(), fa[i]);
        if (it == fb.end()) continue;
        int j = int(it - fb.begin());
        // misura contiguità
        int len = 0;
        while (i + len < (int)fa.size() && j + len < (int)fb.size()
               && fa[i + len] == fb[j + len]) {
            ++len;
        }
        if (len > best_len) {
            best_len = len;
            bestA = ia[i];
            bestB = ib[j];
        }
    }

    if (best_len >= min_overlap) return {best_len, bestA, bestB};
    return {0,0,0};
}

/* =========================  METODO 3: PSH  =========================
 * Progressive Seed Hopping: indicizza le posizioni di fb per valore e
 * estende contiguamente da ogni seed. Restituisce il best (len, startA, startB).
 */
static std::tuple<int,int,int>
progressive_seed_hopping(const ProcessedRead &A,
                         const ProcessedRead &B,
                         int /*k*/, int min_overlap)
{
    const auto &fa = A.comp.comp_fp;
    const auto &fb = B.comp.comp_fp;
    const auto &ia = A.comp.comp_indices;
    const auto &ib = B.comp.comp_indices;

    int n1 = (int)fa.size(), n2 = (int)fb.size();
    if (!n1 || !n2) return {0,0,0};

    static thread_local std::unordered_map<unsigned int, std::vector<int>> posB;
    posB.clear();
    posB.reserve(n2);
    for (int j = 0; j < n2; ++j) posB[fb[j]].push_back(j);

    int best_len = 0, best_x = 0, best_y = 0;

    for (int i = 0; i < n1; ++i) {
        auto it = posB.find(fa[i]);
        if (it == posB.end()) continue;
        for (int j0 : it->second) {
            int len = 1, x = i + 1, y = j0 + 1;
            while (x < n1 && y < n2 && fa[x] == fb[y]) { ++len; ++x; ++y; }
            if (len > best_len) {
                best_len = len;
                best_x = i;
                best_y = j0;
                if (best_len >= min_overlap) break;
            }
        }
        if (best_len >= min_overlap) break;
    }

    if (best_len >= min_overlap)
        return {best_len, ia[best_x], ib[best_y]};
    return {0,0,0};
}

/* =========================  SUFFIX AUTOMATON  =========================
    * Costruisce il suffix automaton di A e lo usa per trovare la LCS con B.
 */
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A)
{
    SuffixAutomaton sa;
    sa.st.reserve(2 * A.size());
    sa.st.push_back({0, -1, 0, {}});
    sa.last = 0;
    for (int i = 0; i < (int)A.size(); ++i) {
        unsigned int c = A[i];
        int cur = int(sa.st.size());
        sa.st.push_back({sa.st[sa.last].len + 1, 0, i, {}});
        int p = sa.last;
        while (p != -1 && !sa.st[p].next.count(c)) {
            sa.st[p].next[c] = cur;
            p = sa.st[p].link;
        }
        if (p == -1) {
            sa.st[cur].link = 0;
        } else {
            int q = sa.st[p].next[c];
            if (sa.st[p].len + 1 == sa.st[q].len) {
                sa.st[cur].link = q;
            } else {
                int clone = int(sa.st.size());
                sa.st.push_back(sa.st[q]);
                sa.st[clone].len = sa.st[p].len + 1;
                while (p != -1 && sa.st[p].next[c] == q) {
                    sa.st[p].next[c] = clone;
                    p = sa.st[p].link;
                }
                sa.st[q].link = sa.st[cur].link = clone;
            }
        }
        sa.last = cur;
    }
    return sa;
}

std::tuple<int, int, int>
match_suffix_automaton(const SuffixAutomaton &sa,
                       const std::vector<unsigned int> &B)
{
    int v = 0, l = 0, best = 0, bestpos = 0, bestst = 0;
    for (int i = 0; i < (int)B.size(); ++i) {
        unsigned int c = B[i];
        if (sa.st[v].next.count(c)) {
            v = sa.st[v].next.at(c);
            ++l;
        } else {
            while (v != -1 && !sa.st[v].next.count(c))
                v = sa.st[v].link;
            if (v == -1) {
                v = 0; l = 0;
            } else {
                l = sa.st[v].len + 1;
                v = sa.st[v].next.at(c);
            }
        }
        if (l > best) {
            best = l; bestpos = i; bestst = v;
        }
    }
    int startB = bestpos - best + 1;
    int startA = sa.st[bestst].first_pos - best + 1;
    return {best, startA, startB};
}

/* =========================  METODO 4: CORE  =========================
 * Chaining semplice: ricava tre seed (da FGOE/AOE/PSH con soglia 1),
 * fa un DP su 3 punti con banda k e restituisce il best.
 * Manteniamo la stessa logica del tuo codice per compatibilità.
 */
static std::tuple<int,int,int>
core(const ProcessedRead &pr1, const ProcessedRead &pr2, int k, int min_overlap)
{
    auto [l1, s1a, s1b] = fingerprint_guided_overlap_extension(pr1, pr2, k, 1);
    auto [l2, s2a, s2b] = adaptive_overlap_extension(pr1, pr2, k, 1);
    auto [l3, s3a, s3b] = progressive_seed_hopping(pr1, pr2, k, 1);

    thread_local std::vector<std::tuple<int, int, int>> seeds;
    seeds.clear();
    if (l1 > 0) seeds.emplace_back(l1, s1a, s1b);
    if (l2 > 0) seeds.emplace_back(l2, s2a, s2b);
    if (l3 > 0) seeds.emplace_back(l3, s3a, s3b);
    if (seeds.empty()) return {0,0,0};

    std::sort(seeds.begin(), seeds.end(),
              [](auto &a, auto &b){ return std::get<1>(a) < std::get<1>(b); });

    int n = (int)seeds.size();
    static thread_local std::vector<int> dp, pred;
    dp.assign(n, 0);
    pred.assign(n, -1);

    int best = 0, bestIdx = 0;
    for (int i = 0; i < n; ++i) {
        dp[i] = std::get<0>(seeds[i]);
        for (int j = 0; j < i; ++j) {
            int la = std::get<0>(seeds[j]);
            int sa = std::get<1>(seeds[j]) + la;
            int sb = std::get<2>(seeds[j]) + la;
            int ea = std::get<1>(seeds[i]), eb = std::get<2>(seeds[i]);
            if (ea >= sa && eb >= sb && ea - sa <= k && eb - sb <= k) {
                int cand = dp[j] + std::get<0>(seeds[i]);
                if (cand > dp[i]) { dp[i] = cand; pred[i] = j; }
            }
        }
        if (dp[i] > best) { best = dp[i]; bestIdx = i; }
    }

    if (best >= min_overlap) {
        int cur = bestIdx;
        while (pred[cur] != -1) cur = pred[cur];
        return {best, std::get<1>(seeds[cur]), std::get<2>(seeds[cur])};
    }

    // fallback sul seed più lungo
    int naive_best = 0, na = 0, nb = 0;
    for (auto &s : seeds) {
        int l = std::get<0>(s);
        if (l > naive_best) { naive_best = l; na = std::get<1>(s); nb = std::get<2>(s); }
    }
    if (naive_best >= min_overlap) return {naive_best, na, nb};
    return {0,0,0};
}

#include <array>   // assicurati che sia incluso

OverlapResult compare_candidate_pair(ReadData &r1,
                                     ReadData &r2,
                                     int k,
                                     int min_overlap,
                                     bool /*verbose*/,          // non più usato
                                     int /*max_repeat_threshold*/)
{
    OverlapResult best; // overlap_len=0 iniziale

    auto build_fingerprint_string = [](const auto &p, int c, int ln) -> std::string {
        std::ostringstream os;
        for (int i = c; i < c + ln; ++i) {
            os << p.comp.comp_fp[i];
            if (i + 1 < c + ln) os << "-";
        }
        return os.str();
    };

    auto try_update = [&](int len, int cA, int cB,
                          const ProcessedRead &pA, const ProcessedRead &pB,
                          const std::string &sA, const std::string &sB,
                          const char *oA, const char *oB,
                          const char *comb, const char *algo)
    {
        if (len <= best.overlap_len || cA < 0 || cB < 0) return;

        const int maxLA = int(pA.comp.comp_indices.size()) - cA;
        const int maxLB = int(pB.comp.comp_indices.size()) - cB;
        if (maxLA <= 0 || maxLB <= 0) return;

        if (len > maxLA || len > maxLB)
            len = std::min({len, maxLA, maxLB});

        if (len < min_overlap) return;

        best.overlap_len  = len;
        best.combination  = comb;
        best.orientation1 = oA;
        best.orientation2 = oB;
        best.start1 = pA.comp.comp_indices[cA];
        best.end1   = pA.comp.comp_indices[cA + len - 1] + k;
        best.start2 = pB.comp.comp_indices[cB];
        best.end2   = pB.comp.comp_indices[cB + len - 1] + k;
        best.r1 = sA;
        best.r2 = sB;
        best.used_algorithm = algo;

        // Costruzione incondizionata per il JSON
        best.fingerprint_r1 = build_fingerprint_string(pA, cA, len);
        best.fingerprint_r2 = build_fingerprint_string(pB, cB, len);
    };

    struct OrInfo {
        const ProcessedRead* pA;
        const ProcessedRead* pB;
        const std::string*   sA;
        const std::string*   sB;
        const char* oA;
        const char* oB;
        const char* comb;
    };

    const std::array<OrInfo,4> ORS = {{
        { &r1.pr_fwd, &r2.pr_fwd, &r1.forward,     &r2.forward,     "forward","forward","ff" },
        { &r1.pr_fwd, &r2.pr_rev, &r1.forward,     &r2.reverse_seq,  "forward","reverse","fr" },
        { &r1.pr_rev, &r2.pr_fwd, &r1.reverse_seq, &r2.forward,      "reverse","forward","rf" },
        { &r1.pr_rev, &r2.pr_rev, &r1.reverse_seq, &r2.reverse_seq,  "reverse","reverse","rr" }
    }};

    // Cascata per orientamento: FGOE -> AOE -> PSH -> CORE
    for (const auto &ornt : ORS) {
        if (ornt.pA->comp.comp_fp.empty() || ornt.pB->comp.comp_fp.empty())
            continue;

        // FGOE
        {
            auto [len, a, b] = fingerprint_guided_overlap_extension(*ornt.pA, *ornt.pB, k, min_overlap);
            if (len >= min_overlap) {
                const int cA = find_comp_index(ornt.pA->comp.comp_indices, a);
                const int cB = find_comp_index(ornt.pB->comp.comp_indices, b);
                if (cA >= 0 && cB >= 0) {
                    try_update(len, cA, cB, *ornt.pA, *ornt.pB, *ornt.sA, *ornt.sB, ornt.oA, ornt.oB, ornt.comb, "FGOE");
                    continue;
                }
            }
        }
        // AOE
        {
            auto [len, a, b] = adaptive_overlap_extension(*ornt.pA, *ornt.pB, k, min_overlap);
            if (len >= min_overlap) {
                const int cA = find_comp_index(ornt.pA->comp.comp_indices, a);
                const int cB = find_comp_index(ornt.pB->comp.comp_indices, b);
                if (cA >= 0 && cB >= 0) {
                    try_update(len, cA, cB, *ornt.pA, *ornt.pB, *ornt.sA, *ornt.sB, ornt.oA, ornt.oB, ornt.comb, "AOE");
                    continue;
                }
            }
        }
        // PSH
        {
            auto [len, a, b] = progressive_seed_hopping(*ornt.pA, *ornt.pB, k, min_overlap);
            if (len >= min_overlap) {
                const int cA = find_comp_index(ornt.pA->comp.comp_indices, a);
                const int cB = find_comp_index(ornt.pB->comp.comp_indices, b);
                if (cA >= 0 && cB >= 0) {
                    try_update(len, cA, cB, *ornt.pA, *ornt.pB, *ornt.sA, *ornt.sB, ornt.oA, ornt.oB, ornt.comb, "PSH");
                    continue;
                }
            }
        }
        // CORE
        {
            auto [len, a, b] = core(*ornt.pA, *ornt.pB, k, min_overlap);
            if (len >= min_overlap) {
                const int cA = find_comp_index(ornt.pA->comp.comp_indices, a);
                const int cB = find_comp_index(ornt.pB->comp.comp_indices, b);
                if (cA >= 0 && cB >= 0) {
                    try_update(len, cA, cB, *ornt.pA, *ornt.pB, *ornt.sA, *ornt.sB, ornt.oA, ornt.oB, ornt.comb, "CORE");
                }
            }
        }
    }

    // Fallback (Suffix Automaton) se necessario
    if (best.overlap_len < min_overlap) {
        bool sa_fwd_built = false, sa_rev_built = false;
        SuffixAutomaton sa_fwd, sa_rev;

        auto get_sa = [&](const ProcessedRead* pA) -> const SuffixAutomaton& {
            if (pA == &r1.pr_fwd) {
                if (!sa_fwd_built) { sa_fwd = build_suffix_automaton(pA->comp.comp_fp); sa_fwd_built = true; }
                return sa_fwd;
            } else {
                if (!sa_rev_built) { sa_rev = build_suffix_automaton(pA->comp.comp_fp); sa_rev_built = true; }
                return sa_rev;
            }
        };

        const int fallback_min = std::max(min_overlap, best.overlap_len + 1);

        for (const auto &ornt : ORS) {
            if (ornt.pA->comp.comp_fp.empty() || ornt.pB->comp.comp_fp.empty()) continue;

            const auto &sa = get_sa(ornt.pA);
            auto [len, a, b] = match_suffix_automaton(sa, ornt.pB->comp.comp_fp);

            if (len >= fallback_min) {
                const int cA = find_comp_index(ornt.pA->comp.comp_indices, a);
                const int cB = find_comp_index(ornt.pB->comp.comp_indices, b);
                if (cA >= 0 && cB >= 0) {
                    try_update(len, cA, cB, *ornt.pA, *ornt.pB, *ornt.sA, *ornt.sB, ornt.oA, ornt.oB, ornt.comb, "SuffixAutomaton");
                }
            }
        }
    }

    return best;
}