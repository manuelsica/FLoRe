// filter.cpp
#include "filter.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <string_view>
#include <numeric>

namespace pseudo_overlap
{

    // --- già viste in precedenza ---
    double entropy(const std::string &region)
    {
        if (region.empty())
            return 0.0;
        int count[256] = {};
        for (unsigned char c : region)
            count[c]++;
        double H = 0.0;
        int N = region.size();
        for (char b : {'A', 'C', 'G', 'T', 'N'})
        {
            int f = count[(unsigned char)b];
            if (!f)
                continue;
            double p = double(f) / N;
            H -= p * std::log2(p);
        }
        return H;
    }

    bool low_complexity(const std::string &region, double min_entropy)
    {
        return entropy(region) < min_entropy;
    }

    static std::pair<int, int> banded_sw_stats(const std::string &r1,
                                               const std::string &r2,
                                               int band)
    {
        int m = r1.size(), n = r2.size();
        if (m == 0 || n == 0)
            return {0, 0};

        std::vector<std::vector<int>> H(m + 1, std::vector<int>(n + 1, 0));
        int max_i = 0, max_j = 0, max_score = 0;

        for (int i = 1; i <= m; i++)
        {
            int j0 = std::max(1, i - band), j1 = std::min(n, i + band);
            for (int j = j0; j <= j1; j++)
            {
                int match_score = (r1[i - 1] == r2[j - 1] ? 2 : -1);
                int s = std::max({0,
                                  H[i - 1][j - 1] + match_score,
                                  H[i - 1][j] - 1,
                                  H[i][j - 1] - 1});
                H[i][j] = s;
                if (s > max_score)
                {
                    max_score = s;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        // traceback
        int i = max_i, j = max_j;
        int matches = 0, length = 0;
        while (i > 0 && j > 0 && H[i][j] > 0)
        {
            int s = H[i][j];
            int diag = H[i - 1][j - 1];
            int up = H[i - 1][j];
            int left = H[i][j - 1];
            if (s == diag + (r1[i - 1] == r2[j - 1] ? 2 : -1))
            {
                if (r1[i - 1] == r2[j - 1])
                    matches++;
                length++;
                i--;
                j--;
            }
            else if (s == up - 1)
            {
                length++;
                i--;
            }
            else
            {
                length++;
                j--;
            }
        }
        return {matches, length};
    }

    // --- FCLA: fingerprint‑chaining + local extension sui gap ---
    bool fingerprint_chained_local_align(const std::string &r1,
                                         const std::string &r2,
                                         int k,
                                         double min_identity,
                                         int max_gap)
    {
        // 1) fingerprint delle region
        auto pr1 = processReadFingerprint(r1, k);
        auto pr2 = processReadFingerprint(r2, k);

        // 2) estrazione seed comuni
        std::unordered_map<unsigned int, std::vector<int>> pos2;
        pos2.reserve(pr2.comp.comp_fp.size());
        for (int j = 0; j < (int)pr2.comp.comp_fp.size(); j++)
            pos2[pr2.comp.comp_fp[j]].push_back(j);

        std::vector<std::pair<int, int>> seeds;
        seeds.reserve(pr1.comp.comp_fp.size());
        for (int i = 0; i < (int)pr1.comp.comp_fp.size(); i++)
        {
            auto it = pos2.find(pr1.comp.comp_fp[i]);
            if (it != pos2.end())
            {
                for (int j : it->second)
                    seeds.emplace_back(i, j);
            }
        }

        // se nessun seed, fallback su intera regione
        if (seeds.empty())
        {
            auto [m, l] = banded_sw_stats(r1, r2, max_gap);
            double identity = l > 0 ? double(m) / l : 0.0;
            return identity >= min_identity;
        }

        // 3) chaining dinamico dei seed
        std::sort(seeds.begin(), seeds.end(),
                  [](auto &a, auto &b)
                  { return a.first < b.first; });
        int S = seeds.size();
        std::vector<int> dp(S), pred(S, -1);
        for (int i = 0; i < S; i++)
            dp[i] = k;

        for (int i = 0; i < S; i++)
        {
            for (int j = 0; j < i; j++)
            {
                // mantieni ordine e distanze entro max_gap
                auto [i1, j1] = seeds[j];
                auto [i2, j2] = seeds[i];
                int di = i2 - i1 - k;
                int dj = j2 - j1 - k;
                if (i1 < i2 && j1 < j2 && di >= 0 && dj >= 0 && di <= max_gap && dj <= max_gap)
                {
                    if (dp[j] + k > dp[i])
                    {
                        dp[i] = dp[j] + k;
                        pred[i] = j;
                    }
                }
            }
        }

        // 4) ricostruzione del best‐chain
        int best_idx = std::max_element(dp.begin(), dp.end()) - dp.begin();
        std::vector<std::pair<int, int>> chain;
        for (int cur = best_idx; cur != -1; cur = pred[cur])
        {
            chain.push_back(seeds[cur]);
        }
        std::reverse(chain.begin(), chain.end());

        // 5) estensione locale sui gap
        int matches_total = k * chain.size();
        int length_total = matches_total;
        int ptr1 = 0, ptr2 = 0;

        for (int t = 0; t <= (int)chain.size(); t++)
        {
            int end1, end2;
            if (t < (int)chain.size())
            {
                end1 = pr1.comp.comp_indices[chain[t].first];
                end2 = pr2.comp.comp_indices[chain[t].second];
            }
            else
            {
                end1 = r1.size();
                end2 = r2.size();
            }
            if (end1 > ptr1 && end2 > ptr2)
            {
                auto sub1 = r1.substr(ptr1, end1 - ptr1);
                auto sub2 = r2.substr(ptr2, end2 - ptr2);
                auto [m, l] = banded_sw_stats(sub1, sub2, max_gap);
                matches_total += m;
                length_total += l;
            }
            if (t < (int)chain.size())
            {
                ptr1 = pr1.comp.comp_indices[chain[t].first] + k;
                ptr2 = pr2.comp.comp_indices[chain[t].second] + k;
            }
        }

        double identity = length_total > 0
                              ? double(matches_total) / length_total
                              : 0.0;
        return identity >= min_identity;
    }

    // --- NUOVO: k-mer frequency + JS divergence ---
    static std::unordered_map<std::string, int>
    compute_kmer_freq(const std::string &s, int k)
    {
        std::unordered_map<std::string, int> freq;
        int n = s.size();
        for (int i = 0; i + k <= n; i++)
        {
            freq[s.substr(i, k)]++;
        }
        return freq;
    }

    static double kldiv(const std::vector<double> &P,
                        const std::vector<double> &Q)
    {
        double d = 0.0;
        for (size_t i = 0; i < P.size(); i++)
        {
            if (P[i] > 0)
                d += P[i] * std::log2(P[i] / (Q[i] > 0 ? Q[i] : 1e-12));
        }
        return d;
    }

    bool spectrum_similarity(const std::string &r1,
                             const std::string &r2,
                             int kmer_len,
                             double max_js)
    {
        auto f1 = compute_kmer_freq(r1, kmer_len);
        auto f2 = compute_kmer_freq(r2, kmer_len);
        // unione chiavi
        std::vector<std::string> keys;
        keys.reserve(f1.size() + f2.size());
        for (auto &kv : f1)
            keys.push_back(kv.first);
        for (auto &kv : f2)
            if (!f1.count(kv.first))
                keys.push_back(kv.first);
        // compongo distribuzioni
        int N1 = r1.size() - kmer_len + 1, N2 = r2.size() - kmer_len + 1;
        std::vector<double> P, Q, M;
        P.reserve(keys.size());
        Q.reserve(keys.size());
        M.reserve(keys.size());
        for (auto &k : keys)
        {
            double p = f1.count(k) ? double(f1[k]) / N1 : 0.0;
            double q = f2.count(k) ? double(f2[k]) / N2 : 0.0;
            P.push_back(p);
            Q.push_back(q);
            M.push_back(0.5 * (p + q));
        }
        double js = 0.5 * kldiv(P, M) + 0.5 * kldiv(Q, M);
        return js <= max_js;
    }

    // --- NUOVO: entropia su blocchi ---
    bool block_entropy_consistency(const std::string &region,
                                   int block_size,
                                   int block_step,
                                   double max_var)
    {
        std::vector<double> H;
        int n = region.size();
        for (int i = 0; i + block_size <= n; i += block_step)
        {
            H.push_back(entropy(region.substr(i, block_size)));
        }
        if (H.empty())
            return true;
        double mean = std::accumulate(H.begin(), H.end(), 0.0) / H.size();
        double var = 0;
        for (double h : H)
            var += (h - mean) * (h - mean);
        var /= H.size();
        return var <= max_var;
    }

} // namespace pseudo_overlap
