// index.cpp
// Implementa le funzioni per la costruzione degli indici invertiti e il pre-filtraggio delle coppie candidate.
#include "index.hpp"
#include "util.hpp"
#include "read.hpp"
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

void buildAllReadsData(std::vector<ReadData> &all_reads,
                       const std::vector<std::string> &reads,
                       int k,
                       bool use_solid_fingerprint,
                       const std::unordered_set<unsigned int> &solid_fingerprint_set)
{
    for (size_t i = 0; i < all_reads.size(); i++)
    {
        all_reads[i].forward = reads[i];
        all_reads[i].reverse_seq = reverse_complement(reads[i]);
        if (use_solid_fingerprint)
        {
            all_reads[i].pr_fwd = processReadSolidFingerprint(all_reads[i].forward, k, solid_fingerprint_set);
            all_reads[i].pr_rev = processReadSolidFingerprint(all_reads[i].reverse_seq, k, solid_fingerprint_set);
        }
        else
        {
            all_reads[i].pr_fwd = processReadFingerprint(all_reads[i].forward, k);
            all_reads[i].pr_rev = processReadFingerprint(all_reads[i].reverse_seq, k);
        }
    }
}

void fillFingerprintSets(std::vector<ReadData> &all_reads)
{
    for (auto &rd : all_reads)
    {
        for (auto val : rd.pr_fwd.comp.comp_fp)
            rd.set_fwd.insert((long long)val);
        for (auto val : rd.pr_rev.comp.comp_fp)
            rd.set_rev.insert((long long)val);
        rd.sorted_fwd.assign(rd.set_fwd.begin(), rd.set_fwd.end());
        rd.sorted_rev.assign(rd.set_rev.begin(), rd.set_rev.end());
        std::sort(rd.sorted_fwd.begin(), rd.sorted_fwd.end());
        std::sort(rd.sorted_rev.begin(), rd.sorted_rev.end());
    }
}

void buildInvertedIndex(const std::vector<ReadData> &all_reads,
                        std::unordered_map<long long, std::vector<int>> &index_fwd,
                        std::unordered_map<long long, std::vector<int>> &index_rev)
{
    for (int i = 0; i < (int)all_reads.size(); i++)
    {
        {
            std::unordered_set<long long> uniq(all_reads[i].pr_fwd.comp.comp_fp.begin(),
                                               all_reads[i].pr_fwd.comp.comp_fp.end());
            for (auto val : uniq)
                index_fwd[val].push_back(i);
        }
        {
            std::unordered_set<long long> uniq(all_reads[i].pr_rev.comp.comp_fp.begin(),
                                               all_reads[i].pr_rev.comp.comp_fp.end());
            for (auto val : uniq)
                index_rev[val].push_back(i);
        }
    }
}

std::vector<Pair> generateCandidatePairs(const std::unordered_map<long long, std::vector<int>> &index_fwd,
                                         const std::unordered_map<long long, std::vector<int>> &index_rev)
{
    std::unordered_set<Pair, PairHash> pairs;
    for (const auto &kv : index_fwd)
    {
        const std::vector<int> &vec = kv.second;
        for (size_t i = 0; i < vec.size(); i++)
            for (size_t j = i+1; j < vec.size(); j++)
            {
                int a = vec[i], b = vec[j];
                if (a > b) std::swap(a, b);
                pairs.insert({a, b});
            }
    }
    for (const auto &kv : index_rev)
    {
        const std::vector<int> &vec = kv.second;
        for (size_t i = 0; i < vec.size(); i++)
            for (size_t j = i+1; j < vec.size(); j++)
            {
                int a = vec[i], b = vec[j];
                if (a > b) std::swap(a, b);
                pairs.insert({a, b});
            }
    }
    for (const auto &kv : index_fwd)
    {
        auto it = index_rev.find(kv.first);
        if (it != index_rev.end())
        {
            const std::vector<int> &vf = kv.second;
            const std::vector<int> &vr = it->second;
            for (int f : vf)
                for (int r : vr)
                {
                    if (f == r) continue;
                    int a = f, b = r;
                    if (a > b) std::swap(a, b);
                    pairs.insert({a, b});
                }
        }
    }
    return std::vector<Pair>(pairs.begin(), pairs.end());
}

int hybrid_sorted_intersection_size(const std::vector<long long> &v1,
                                    const std::vector<long long> &v2)
{
    if (v1.size() < 20 || v2.size() < 20)
    {
        int cnt = 0;
        if (v1.size() < v2.size())
            for (auto &x : v1)
                if (std::binary_search(v2.begin(), v2.end(), x))
                    cnt++;
        else
            for (auto &x : v2)
                if (std::binary_search(v1.begin(), v1.end(), x))
                    cnt++;
        return cnt;
    }
    int i = 0, j = 0, cnt = 0;
    while (i < (int)v1.size() && j < (int)v2.size())
    {
        if (v1[i] < v2[j])
            i++;
        else if (v2[j] < v1[i])
            j++;
        else { cnt++; i++; j++; }
    }
    return cnt;
}