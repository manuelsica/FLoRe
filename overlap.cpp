// overlap.cpp
// Implementa il calcolo degli overlap tramite:
// - FGOE: Fingerprint-Guided Overlap Extension
// - adaptive_overlap_extension: Adaptive Overlap Extension migliorato con std::mismatch
// - KHS: Kmer Hopping Search
// - core: Combined Overlap Refinement Engine (precedentemente COIN) che utilizza chaining dinamico
// - fallback: Suffix Automaton per il longest common substring
#include "overlap.hpp"
#include "index.hpp"
#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <cmath>
#include <sstream>
#include <iterator>

//----------------------------------------------------------------------------
// Metodo 1: Fingerprint-Guided Overlap Extension (FGOE)
//----------------------------------------------------------------------------
static std::tuple<int,int,int> fingerprint_guided_overlap_extension(const ProcessedRead &pr1,
                                                                     const ProcessedRead &pr2,
                                                                     int k,
                                                                     int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int size1 = (int)fp1.size(), size2 = (int)fp2.size();
    if (size1 == 0 || size2 == 0)
        return {0,0,0};
    int i = size1 - 1, j = 0, length_match = 0;
    while (i >= 0 && j < size2 && fp1[i] == fp2[j])
    {
        length_match++;
        i--;
        j++;
    }
    if (length_match >= min_overlap)
    {
        int startA = pr1.comp.comp_indices[i+1];
        int startB = pr2.comp.comp_indices[0];
        return {length_match, startA, startB};
    }
    return {0,0,0};
}

//----------------------------------------------------------------------------
// Metodo 2: adaptive_overlap_extension (AOE) - Versione migliorata
//----------------------------------------------------------------------------
static std::tuple<int,int,int> adaptive_overlap_extension(const ProcessedRead &pr1,
                                                          const ProcessedRead &pr2,
                                                          int k,
                                                          int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp, &fp2 = pr2.comp.comp_fp;
    int n1 = (int)fp1.size(), n2 = (int)fp2.size();
    int best_len = 0;
    int best_startA = 0;
    int best_startB = 0;
    for (int i = 0; i < n1; i++) {
        auto it = std::find(fp2.begin(), fp2.end(), fp1[i]);
        if (it == fp2.end())
            continue;
        int j = it - fp2.begin();
        auto pair_it = std::mismatch(fp1.begin() + i, fp1.end(), fp2.begin() + j, fp2.end());
        int match_len = pair_it.first - (fp1.begin() + i);
        if (match_len > best_len)
        {
            best_len = match_len;
            best_startA = pr1.comp.comp_indices[i];
            best_startB = pr2.comp.comp_indices[j];
        }
    }
    if (best_len >= min_overlap)
        return {best_len, best_startA, best_startB};
    return {0,0,0};
}

//----------------------------------------------------------------------------
// Metodo 3: Kmer Hopping Search (KHS)
//----------------------------------------------------------------------------
static std::tuple<int,int,int> kmer_hopping_search(const ProcessedRead &pr1,
                                                   const ProcessedRead &pr2,
                                                   int k,
                                                   int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp, &fp2 = pr2.comp.comp_fp;
    int n1 = (int)fp1.size(), n2 = (int)fp2.size();
    if (n1 == 0 || n2 == 0)
        return {0,0,0};
    std::unordered_map<unsigned int, std::vector<int>> positionsInB;
    positionsInB.reserve(n2);
    for (int i = 0; i < n2; i++)
        positionsInB[fp2[i]].push_back(i);
    int best_len = 0, best_x = 0, best_y = 0;
    for (int i = 0; i < n1; i++)
    {
        auto it = positionsInB.find(fp1[i]);
        if (it != positionsInB.end())
        {
            for (int posB : it->second)
            {
                int length_match = 1, xx = i + 1, yy = posB + 1;
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

//----------------------------------------------------------------------------
// Metodo 4: core (Combined Overlap Refinement Engine)
// Utilizza chaining dinamico per combinare i seed ottenuti dai metodi FGOE, AOE e KHS.
//----------------------------------------------------------------------------
struct Seed {
    int startA;
    int startB;
    int length;
};

static std::tuple<int,int,int> core(const ProcessedRead &pr1,
                                     const ProcessedRead &pr2,
                                     int k,
                                     int min_overlap)
{
    auto [fgoe_len, fgoe_startA, fgoe_startB] = fingerprint_guided_overlap_extension(pr1, pr2, k, 1);
    auto [aoe_len, aoe_startA, aoe_startB] = adaptive_overlap_extension(pr1, pr2, k, 1);
    auto [khs_len, khs_startA, khs_startB] = kmer_hopping_search(pr1, pr2, k, 1);
    
    std::vector<Seed> seeds;
    if (fgoe_len > 0)
        seeds.push_back({fgoe_startA, fgoe_startB, fgoe_len});
    if (aoe_len > 0)
        seeds.push_back({aoe_startA, aoe_startB, aoe_len});
    if (khs_len > 0)
        seeds.push_back({khs_startA, khs_startB, khs_len});
    if (seeds.empty())
        return {0,0,0};
    
    std::sort(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) {
        return a.startA < b.startA;
    });
    
    int n = seeds.size();
    std::vector<int> dp(n, 0);
    std::vector<int> pred(n, -1);
    
    for (int i = 0; i < n; i++)
        dp[i] = seeds[i].length;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            int gapA = seeds[i].startA - (seeds[j].startA + seeds[j].length);
            int gapB = seeds[i].startB - (seeds[j].startB + seeds[j].length);
            if (gapA >= 0 && gapB >= 0 && gapA <= k && gapB <= k)
            {
                if (dp[j] + seeds[i].length > dp[i])
                {
                    dp[i] = dp[j] + seeds[i].length;
                    pred[i] = j;
                }
            }
        }
    }
    
    int best = 0, bestIdx = 0;
    for (int i = 0; i < n; i++)
    {
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
        int best_startA = seeds[cur].startA;
        int best_startB = seeds[cur].startB;
        return {best, best_startA, best_startB};
    }
    
    int naive_best = 0, naive_A = 0, naive_B = 0;
    for (auto &s : seeds)
    {
        if (s.length > naive_best)
        {
            naive_best = s.length;
            naive_A = s.startA;
            naive_B = s.startB;
        }
    }
    if (naive_best >= min_overlap)
        return {naive_best, naive_A, naive_B};
    return {0,0,0};
}

//----------------------------------------------------------------------------
// Suffix Automaton: Costruisce l'automa per un vettore di fingerprint
//----------------------------------------------------------------------------
SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int> &A)
{
    SuffixAutomaton sa;
    sa.st.reserve(2 * A.size());
    sa.last = 0;
    State initial;
    initial.len = 0;
    initial.link = -1;
    initial.first_pos = 0;
    sa.st.push_back(initial);
    for (int i = 0; i < (int)A.size(); i++)
    {
        unsigned int c = A[i];
        int cur = (int)sa.st.size();
        State stCur;
        stCur.len = sa.st[sa.last].len + 1;
        stCur.link = 0;
        stCur.first_pos = i;
        sa.st.push_back(stCur);
        int p = sa.last;
        while (p != -1 && sa.st[p].next.find(c) == sa.st[p].next.end())
        {
            sa.st[p].next[c] = cur;
            p = sa.st[p].link;
        }
        if (p == -1)
            sa.st[cur].link = 0;
        else
        {
            int q = sa.st[p].next[c];
            if (sa.st[p].len + 1 == sa.st[q].len)
                sa.st[cur].link = q;
            else
            {
                int clone = (int)sa.st.size();
                State stClone;
                stClone.len = sa.st[p].len + 1;
                stClone.next = sa.st[q].next;
                stClone.link = sa.st[q].link;
                stClone.first_pos = sa.st[q].first_pos;
                sa.st.push_back(stClone);
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

//----------------------------------------------------------------------------
// match_suffix_automaton: Trova il longest common substring tra il Suffix Automaton e il vettore B.
//----------------------------------------------------------------------------
std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B)
{
    int v = 0, l = 0, best = 0, bestposB = 0, best_state = 0;
    for (int i = 0; i < (int)B.size(); i++)
    {
        unsigned int c = B[i];
        if (sa.st[v].next.find(c) != sa.st[v].next.end())
        {
            v = sa.st[v].next.at(c);
            l++;
        }
        else
        {
            while (v != -1 && sa.st[v].next.find(c) == sa.st[v].next.end())
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
            bestposB = i;
            best_state = v;
        }
    }
    int startB = bestposB - best + 1;
    int startA = sa.st[best_state].first_pos - best + 1;
    return {best, startA, startB};
}

std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int> &A,
    const std::vector<unsigned int> &B)
{
    SuffixAutomaton sa = build_suffix_automaton(A);
    return match_suffix_automaton(sa, B);
}

//----------------------------------------------------------------------------
// compare_candidate_pair: Combina i metodi per tutti gli orientamenti e usa il fallback con Suffix Automaton se necessario.
// I metodi usati sono FGOE, adaptive_overlap_extension, KHS e core (Combined Overlap Refinement Engine).
//----------------------------------------------------------------------------
OverlapResult compare_candidate_pair(ReadData &r1,
                                     ReadData &r2,
                                     int k,
                                     int min_overlap,
                                     bool /*verbose*/,
                                     int /*max_repeat_threshold*/)
{
    OverlapResult best;
    
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
            auto getFingerprintRegion = [&](const ProcessedRead &p, int cstart, int length) {
                std::ostringstream oss;
                int stop = std::min<int>(cstart + length, (int)p.comp.comp_fp.size());
                for (int i = cstart; i < stop; i++) {
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
    
    auto compute_best_overlap_4stages = [&](ProcessedRead &A, ProcessedRead &B,
                                             std::string &rA, std::string &rB,
                                             const std::string &oA, const std::string &oB,
                                             const std::string &comb)
    {
        // Metodo FGOE
        {
            auto [matchLen, startA, startB] = fingerprint_guided_overlap_extension(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++) {
                    if (A.comp.comp_indices[i] == startA) { idxA = i; break; }
                }
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++) {
                    if (B.comp.comp_indices[j] == startB) { idxB = j; break; }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "FGOE");
                return;
            }
        }
        // Metodo adaptive_overlap_extension (AOE)
        {
            auto [matchLen, startA, startB] = adaptive_overlap_extension(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++) {
                    if (A.comp.comp_indices[i] == startA) { idxA = i; break; }
                }
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++) {
                    if (B.comp.comp_indices[j] == startB) { idxB = j; break; }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "AOE");
                return;
            }
        }
        // Metodo KHS
        {
            auto [matchLen, startA, startB] = kmer_hopping_search(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++) {
                    if (A.comp.comp_indices[i] == startA) { idxA = i; break; }
                }
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++) {
                    if (B.comp.comp_indices[j] == startB) { idxB = j; break; }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "KHS");
                return;
            }
        }
        // Metodo core (Combined Overlap Refinement Engine)
        {
            auto [matchLen, startA, startB] = core(A, B, k, min_overlap);
            if (matchLen >= min_overlap)
            {
                int idxA = 0, idxB = 0;
                for (int i = 0; i < (int)A.comp.comp_indices.size(); i++) {
                    if (A.comp.comp_indices[i] == startA) { idxA = i; break; }
                }
                for (int j = 0; j < (int)B.comp.comp_indices.size(); j++) {
                    if (B.comp.comp_indices[j] == startB) { idxB = j; break; }
                }
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "Combined Overlap Refinement Engine");
                return;
            }
        }
    };
    
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
            auto getFingerprintRegion = [&](const ProcessedRead &p, int cstart, int ln) {
                std::ostringstream oss;
                int stop = std::min<int>(cstart + ln, (int)p.comp.comp_fp.size());
                for (int i = cstart; i < stop; i++) {
                    oss << p.comp.comp_fp[i];
                    if (i < stop - 1) oss << "-";
                }
                return oss.str();
            };
            best.fingerprint_r1 = getFingerprintRegion(A, c1, length);
            best.fingerprint_r2 = getFingerprintRegion(B, c2, length);
            best.used_algorithm = "SuffixAutomaton";
        }
    };

    // Esegue i metodi per i quattro orientamenti: ff, fr, rf, rr
    compute_best_overlap_4stages(r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward, "forward", "forward", "ff");
    compute_best_overlap_4stages(r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq, "forward", "reverse", "fr");
    compute_best_overlap_4stages(r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward, "reverse", "forward", "rf");
    compute_best_overlap_4stages(r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq, "reverse", "reverse", "rr");
    
    if (best.overlap_len < min_overlap)
    {
        compute_best_overlap_suffix_automaton(r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward, "forward", "forward", "ff", r1.sa_built_fwd, r1.suffix_automaton_fwd);
        compute_best_overlap_suffix_automaton(r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq, "forward", "reverse", "fr", r1.sa_built_fwd, r1.suffix_automaton_fwd);
        compute_best_overlap_suffix_automaton(r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward, "reverse", "forward", "rf", r1.sa_built_rev, r1.suffix_automaton_rev);
        compute_best_overlap_suffix_automaton(r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq, "reverse", "reverse", "rr", r1.sa_built_rev, r1.suffix_automaton_rev);
    }
    
    return best;
}