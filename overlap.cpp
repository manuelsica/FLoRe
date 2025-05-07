// overlap.cpp
#include "overlap.hpp"
#include "index.hpp"
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <cassert>

// Cerca “x” esattamente in idx; ritorna l’indice o -1 se non trovato
static int find_comp_index(const std::vector<int> &idx, int x)
{
    auto it = std::lower_bound(idx.begin(), idx.end(), x);
    if (it != idx.end() && *it == x)
        return int(std::distance(idx.begin(), it));
    return -1;
}

//----------------------------------------------------------------------------
// Metodo 1: Fingerprint-Guided Overlap Extension (FGOE)
static std::tuple<int, int, int>
fingerprint_guided_overlap_extension(const ProcessedRead &pr1,
                                     const ProcessedRead &pr2,
                                     int /*k*/, int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;

    // 1) Se la soglia richiesta è > 3, non possiamo soddisfarla qui
    if (min_overlap > 3)
        return {0, 0, 0};

    // 2) Serve almeno 3 valori in ciascun vettore
    if (fp1.size() < 3 || fp2.size() < 3)
        return {0, 0, 0};

    // 3) Confronto degli ultimi 3 di fp1 con i primi 3 di fp2
    int n1 = (int)fp1.size();
    for (int t = 0; t < 3; ++t)
    {
        if (fp1[n1 - 3 + t] != fp2[t])
            return {0, 0, 0};
    }

    // 4) Se siamo arrivati fin qui, len = 3 >= min_overlap
    int start1 = pr1.comp.comp_indices[n1 - 3];
    int start2 = pr2.comp.comp_indices[0];
    return {3, start1, start2};
}
//----------------------------------------------------------------------------
// Metodo 2: Adaptive Overlap Extension (AOE)
static std::tuple<int, int, int>
adaptive_overlap_extension(const ProcessedRead &pr1,
                           const ProcessedRead &pr2,
                           int /*k*/, int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int best_len = 0, bestA = 0, bestB = 0;
    for (int i = 0; i < (int)fp1.size(); ++i)
    {
        auto it = std::find(fp2.begin(), fp2.end(), fp1[i]);
        if (it == fp2.end())
            continue;
        int j = int(it - fp2.begin());
        auto mis = std::mismatch(fp1.begin() + i, fp1.end(), fp2.begin() + j);
        int len = int(mis.first - (fp1.begin() + i));
        if (len > best_len)
        {
            best_len = len;
            bestA = pr1.comp.comp_indices[i];
            bestB = pr2.comp.comp_indices[j];
        }
    }
    if (best_len >= min_overlap)
        return {best_len, bestA, bestB};
    return {0, 0, 0};
}

//----------------------------------------------------------------------------
// Metodo 3: Progressive Seed Hopping (PSH)
static std::tuple<int, int, int>
progressive_seed_hopping(const ProcessedRead &pr1,
                         const ProcessedRead &pr2,
                         int /*k*/, int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int n1 = int(fp1.size()), n2 = int(fp2.size());
    if (!n1 || !n2)
        return {0, 0, 0};

    static thread_local std::unordered_map<unsigned int, std::vector<int>> posB;
    posB.clear();
    posB.reserve(n2);
    for (int j = 0; j < n2; ++j)
        posB[fp2[j]].push_back(j);

    int best_len = 0, best_x = 0, best_y = 0;
    for (int i = 0; i < n1; ++i)
    {
        auto it = posB.find(fp1[i]);
        if (it == posB.end())
            continue;
        for (int j0 : it->second)
        {
            int len = 1, x = i + 1, y = j0 + 1;
            while (x < n1 && y < n2 && fp1[x] == fp2[y])
            {
                ++len;
                ++x;
                ++y;
            }
            if (len > best_len)
            {
                best_len = len;
                best_x = i;
                best_y = j0;
                if (best_len >= min_overlap)
                    break;
            }
        }
        if (best_len >= min_overlap)
            break;
    }
    if (best_len >= min_overlap)
    {
        return {best_len,
                pr1.comp.comp_indices[best_x],
                pr2.comp.comp_indices[best_y]};
    }
    return {0, 0, 0};
}

//----------------------------------------------------------------------------
// Metodo 4: Combined Overlap Refinement Engine (CORE)
static std::tuple<int, int, int>
core(const ProcessedRead &pr1, const ProcessedRead &pr2, int k, int min_overlap)
{
    auto [l1, s1a, s1b] = fingerprint_guided_overlap_extension(pr1, pr2, k, 1);
    auto [l2, s2a, s2b] = adaptive_overlap_extension(pr1, pr2, k, 1);
    auto [l3, s3a, s3b] = progressive_seed_hopping(pr1, pr2, k, 1);

    thread_local std::vector<std::tuple<int, int, int>> seeds;
    seeds.clear();
    if (l1 > 0)
        seeds.emplace_back(l1, s1a, s1b);
    if (l2 > 0)
        seeds.emplace_back(l2, s2a, s2b);
    if (l3 > 0)
        seeds.emplace_back(l3, s3a, s3b);
    if (seeds.empty())
        return {0, 0, 0};

    std::sort(seeds.begin(), seeds.end(),
              [](auto &a, auto &b)
              {
                  return std::get<1>(a) < std::get<1>(b);
              });

    int n = int(seeds.size());
    static thread_local std::vector<int> dp, pred;
    dp.assign(n, 0);
    pred.assign(n, -1);

    int best = 0, bestIdx = 0;
    for (int i = 0; i < n; ++i)
    {
        dp[i] = std::get<0>(seeds[i]);
        for (int j = 0; j < i; ++j)
        {
            int la = std::get<0>(seeds[j]);
            int sa = std::get<1>(seeds[j]) + la;
            int sb = std::get<2>(seeds[j]) + la;
            int ea = std::get<1>(seeds[i]), eb = std::get<2>(seeds[i]);
            if (ea >= sa && eb >= sb && ea - sa <= k && eb - sb <= k)
            {
                int cand = dp[j] + std::get<0>(seeds[i]);
                if (cand > dp[i])
                {
                    dp[i] = cand;
                    pred[i] = j;
                }
            }
        }
        if (dp[i] > best)
        {
            best = dp[i];
            bestIdx = i;
        }
    }

    if (best >= min_overlap)
    {
        int cur = bestIdx;
        while (pred[cur] != -1)
            cur = pred[cur];
        return {best,
                std::get<1>(seeds[cur]),
                std::get<2>(seeds[cur])};
    }

    // fallback sul seed più lungo
    int naive_best = 0, na = 0, nb = 0;
    for (auto &s : seeds)
    {
        int l = std::get<0>(s);
        if (l > naive_best)
        {
            naive_best = l;
            na = std::get<1>(s);
            nb = std::get<2>(s);
        }
    }
    if (naive_best >= min_overlap)
        return {naive_best, na, nb};

    return {0, 0, 0};
}

//----------------------------------------------------------------------------
// Costruzione del suffix automaton
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A)
{
    SuffixAutomaton sa;
    sa.st.reserve(2 * A.size());
    sa.st.push_back({0, -1, 0, {}});
    sa.last = 0;
    for (int i = 0; i < (int)A.size(); ++i)
    {
        unsigned int c = A[i];
        int cur = int(sa.st.size());
        sa.st.push_back({sa.st[sa.last].len + 1, 0, i, {}});
        int p = sa.last;
        while (p != -1 && !sa.st[p].next.count(c))
        {
            sa.st[p].next[c] = cur;
            p = sa.st[p].link;
        }
        if (p == -1)
        {
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
                int clone = int(sa.st.size());
                sa.st.push_back(sa.st[q]);
                sa.st[clone].len = sa.st[p].len + 1;
                while (p != -1 && sa.st[p].next[c] == q)
                {
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

//----------------------------------------------------------------------------
// match_suffix_automaton: unica definizione
std::tuple<int, int, int>
match_suffix_automaton(const SuffixAutomaton &sa,
                       const std::vector<unsigned int> &B)
{
    int v = 0, l = 0, best = 0, bestpos = 0, bestst = 0;
    for (int i = 0; i < (int)B.size(); ++i)
    {
        unsigned int c = B[i];
        if (sa.st[v].next.count(c))
        {
            v = sa.st[v].next.at(c);
            ++l;
        }
        else
        {
            while (v != -1 && !sa.st[v].next.count(c))
                v = sa.st[v].link;
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
            bestpos = i;
            bestst = v;
        }
    }
    int startB = bestpos - best + 1;
    int startA = sa.st[bestst].first_pos - best + 1;
    return {best, startA, startB};
}

//----------------------------------------------------------------------------
// Wrapper LCS
std::tuple<int, int, int>
longest_common_substring_suffix_automaton(const std::vector<unsigned int> &A,
                                          const std::vector<unsigned int> &B)
{
    return match_suffix_automaton(build_suffix_automaton(A), B);
}

//----------------------------------------------------------------------------
// compare_candidate_pair: coordina tutti i metodi
OverlapResult compare_candidate_pair(ReadData &r1,
                                     ReadData &r2,
                                     int k,
                                     int min_overlap,
                                     bool verbose,
                                     int max_repeat_threshold)
{
    OverlapResult best;

    auto try_update = [&](int len, int cA, int cB,
                          const ProcessedRead &pA, const ProcessedRead &pB,
                          const std::string &sA, const std::string &sB,
                          const std::string &oA, const std::string &oB,
                          const std::string &comb, const std::string &algo)
    {
        if (len <= best.overlap_len || cA < 0 || cB < 0)
            return;

        /* ------------------------------------------------------------------
         * 1.  Verifica che cA,cB siano indici validi.
         * 2.  “Clippa” len al numero reale di fingerprint rimasti
         *     nei due vettori (così non possiamo mai uscire dal range).
         * ----------------------------------------------------------------- */
        int maxLA = int(pA.comp.comp_indices.size()) - cA;
        int maxLB = int(pB.comp.comp_indices.size()) - cB;
        if (maxLA <= 0 || maxLB <= 0)
            return; // non ci sono fingerprint utilizzabili

        if (len > maxLA || len > maxLB)
            len = std::min({len, maxLA, maxLB});

        /*  Se ora non arriva alla soglia richiesta, scartiamo.               */
        if (len < min_overlap)
            return;

        /*  Da qui in poi siamo sicuri che (cA + len − 1) e (cB + len − 1)
         *  rientrano nei rispettivi vettori.                                 */

        best.overlap_len = len;
        best.combination = comb;
        best.orientation1 = oA;
        best.orientation2 = oB;
        best.start1 = pA.comp.comp_indices[cA];
        best.end1 = pA.comp.comp_indices[cA + len - 1] + k;
        best.start2 = pB.comp.comp_indices[cB];
        best.end2 = pB.comp.comp_indices[cB + len - 1] + k;
        best.r1 = sA;
        best.r2 = sB;

        auto getFP = [&](const auto &p, int c, int ln)
        {
            std::ostringstream os;
            for (int i = c; i < c + ln; ++i)
                os << p.comp.comp_fp[i] << (i < c + ln - 1 ? "-" : "");
            return os.str();
        };
        best.fingerprint_r1 = getFP(pA, cA, len);
        best.fingerprint_r2 = getFP(pB, cB, len);
        best.used_algorithm = algo;
    };

    // 4 orientamenti: ff, fr, rf, rr
    // forward-forward
    {
        auto [len, a, b] = fingerprint_guided_overlap_extension(r1.pr_fwd, r2.pr_fwd, k, min_overlap);
        if (len >= min_overlap)
        {
            int cA = find_comp_index(r1.pr_fwd.comp.comp_indices, a);
            int cB = find_comp_index(r2.pr_fwd.comp.comp_indices, b);
            if (cA >= 0 && cB >= 0)
                try_update(len, cA, cB, r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward, "forward", "forward", "ff", "FGOE");
            else
                return best;
        }
    }
    // forward-reverse
    {
        auto [len, a, b] = adaptive_overlap_extension(r1.pr_fwd, r2.pr_rev, k, min_overlap);
        if (len >= min_overlap)
        {
            int cA = find_comp_index(r1.pr_fwd.comp.comp_indices, a);
            int cB = find_comp_index(r2.pr_rev.comp.comp_indices, b);
            if (cA >= 0 && cB >= 0)
                try_update(len, cA, cB, r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq, "forward", "reverse", "fr", "AOE");
            else
                return best;
        }
    }
    // reverse-forward
    {
        auto [len, a, b] = progressive_seed_hopping(r1.pr_rev, r2.pr_fwd, k, min_overlap);
        if (len >= min_overlap)
        {
            int cA = find_comp_index(r1.pr_rev.comp.comp_indices, a);
            int cB = find_comp_index(r2.pr_fwd.comp.comp_indices, b);
            if (cA >= 0 && cB >= 0)
                try_update(len, cA, cB, r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward, "reverse", "forward", "rf", "PSH");
            else
                return best;
        }
    }
    // reverse-reverse
    {
        auto [len, a, b] = core(r1.pr_rev, r2.pr_rev, k, min_overlap);
        if (len >= min_overlap)
        {
            int cA = find_comp_index(r1.pr_rev.comp.comp_indices, a);
            int cB = find_comp_index(r2.pr_rev.comp.comp_indices, b);
            if (cA >= 0 && cB >= 0)
                try_update(len, cA, cB, r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq, "reverse", "reverse", "rr", "CORE");
            else
                return best;
        }
    }

    // se nessuno ha superato min_overlap, fallback suffix automaton
    if (best.overlap_len < min_overlap)
    {
        // ff
        {
            thread_local bool built = false; // <-- was static
            thread_local SuffixAutomaton sa; // <-- was static
            if (!built)
            {
                sa = build_suffix_automaton(r1.pr_fwd.comp.comp_fp);
                built = true;
            }
            auto [len, a, b] = match_suffix_automaton(sa, r2.pr_fwd.comp.comp_fp);
            if (len >= min_overlap)
            {
                int cA = find_comp_index(r1.pr_fwd.comp.comp_indices, a);
                int cB = find_comp_index(r2.pr_fwd.comp.comp_indices, b);
                if (cA >= 0 && cB >= 0)
                    try_update(len, cA, cB, r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward, "forward", "forward", "ff", "SuffixAutomaton");
            }
        }
        // fr
        {
            thread_local bool built = false;
            thread_local SuffixAutomaton sa;
            if (!built)
            {
                sa = build_suffix_automaton(r1.pr_fwd.comp.comp_fp);
                built = true;
            }
            auto [len, a, b] = match_suffix_automaton(sa, r2.pr_rev.comp.comp_fp);
            if (len >= min_overlap)
            {
                int cA = find_comp_index(r1.pr_fwd.comp.comp_indices, a);
                int cB = find_comp_index(r2.pr_rev.comp.comp_indices, b);
                if (cA >= 0 && cB >= 0)
                    try_update(len, cA, cB, r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq, "forward", "reverse", "fr", "SuffixAutomaton");
            }
        }
        // rf
        {
            thread_local bool built = false;
            thread_local SuffixAutomaton sa;
            if (!built)
            {
                sa = build_suffix_automaton(r1.pr_fwd.comp.comp_fp);
                built = true;
            }
            auto [len, a, b] = match_suffix_automaton(sa, r2.pr_rev.comp.comp_fp);
            if (len >= min_overlap)
            {
                int cA = find_comp_index(r1.pr_rev.comp.comp_indices, a);
                int cB = find_comp_index(r2.pr_fwd.comp.comp_indices, b);
                if (cA >= 0 && cB >= 0)
                    try_update(len, cA, cB, r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward, "reverse", "forward", "rf", "SuffixAutomaton");
            }
        }
        // rr
        {
            thread_local bool built = false;
            thread_local SuffixAutomaton sa;
            if (!built)
            {
                sa = build_suffix_automaton(r1.pr_fwd.comp.comp_fp);
                built = true;
            }
            auto [len, a, b] = match_suffix_automaton(sa, r2.pr_rev.comp.comp_fp);
            if (len >= min_overlap)
            {
                int cA = find_comp_index(r1.pr_rev.comp.comp_indices, a);
                int cB = find_comp_index(r2.pr_rev.comp.comp_indices, b);
                if (cA >= 0 && cB >= 0)
                    try_update(len, cA, cB, r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq, "reverse", "reverse", "rr", "SuffixAutomaton");
            }
        }
    }

    return best;
}
