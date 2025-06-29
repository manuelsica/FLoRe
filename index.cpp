// index.cpp
#include "index.hpp"
#include <algorithm>
#include "util.hpp" 
#include "cfl.hpp"
#include "icfl.hpp"
#include "cfl_icfl.hpp"
#include <cassert>

void buildAllReadsData(std::vector<ReadData>               &all_reads,
    const std::vector<std::string>      &reads,
    int                                  k,
    bool                                 use_solid_fp,
    const std::unordered_set<unsigned int> &solid_set,
    bool                                 use_cfl,
    int                                  cfl_thr,
    bool                                 use_cfl_icfl,
    int                                  cfl_icfl_thr,
    bool                                 use_icfl,
    int                                  icfl_thr)
{
assert(all_reads.size() == reads.size());

for (size_t i = 0; i < reads.size(); ++i)
{
ReadData &rd = all_reads[i];

// 1) forward / reverse
rd.forward     = reads[i];
rd.reverse_seq = reverse_complement(reads[i]);

// 2) calcolo fingerprint
if (use_cfl_icfl) {
rd.pr_fwd = processReadCFLICFL(rd.forward,     k, cfl_icfl_thr);
rd.pr_rev = processReadCFLICFL(rd.reverse_seq, k, cfl_icfl_thr);
}
else if (use_icfl) {
rd.pr_fwd = processReadICFL(rd.forward,     k, icfl_thr);
rd.pr_rev = processReadICFL(rd.reverse_seq, k, icfl_thr);
}
else if (use_cfl) {
rd.pr_fwd = processReadCFL(rd.forward,     k, cfl_thr);
rd.pr_rev = processReadCFL(rd.reverse_seq, k, cfl_thr);
}
else if (use_solid_fp) {
rd.pr_fwd = processReadSolidFingerprint(rd.forward,     k, solid_set);
rd.pr_rev = processReadSolidFingerprint(rd.reverse_seq, k, solid_set);
}
else {
rd.pr_fwd = processReadFingerprint(rd.forward,     k);
rd.pr_rev = processReadFingerprint(rd.reverse_seq, k);
}

// 3) primo fingerprint compresso (Compatibility Property filter)
rd.first_fp_fwd = rd.pr_fwd.comp.comp_fp.empty()
      ? 0
      : rd.pr_fwd.comp.comp_fp.front();
rd.first_fp_rev = rd.pr_rev.comp.comp_fp.empty()
      ? 0
      : rd.pr_rev.comp.comp_fp.front();
}
}
void fillFingerprintSets(std::vector<ReadData> &all_reads)
{
    for (auto &rd : all_reads)
    {
        for (auto val : rd.pr_fwd.comp.comp_fp)
            rd.set_fwd.insert(val);
        for (auto val : rd.pr_rev.comp.comp_fp)
            rd.set_rev.insert(val);

        rd.sorted_fwd.assign(rd.set_fwd.begin(), rd.set_fwd.end());
        rd.sorted_rev.assign(rd.set_rev.begin(), rd.set_rev.end());
        std::sort(rd.sorted_fwd.begin(), rd.sorted_fwd.end());
        std::sort(rd.sorted_rev.begin(), rd.sorted_rev.end());
    }
}

// buildInvertedIndex ottimizzato
void buildInvertedIndex(const std::vector<ReadData> &all_reads,
                        std::unordered_map<long long, std::vector<int>> &index_fwd,
                        std::unordered_map<long long, std::vector<int>> &index_rev)
{
    // Stima il numero totale di chiavi per riservare spazio
    size_t total_keys = 0;
    for (const auto &rd : all_reads)
        total_keys += rd.sorted_fwd.size() + rd.sorted_rev.size();
    index_fwd.reserve(total_keys);
    index_rev.reserve(total_keys);

    // Costruisci il posting list usando solo le chiavi uniche in sorted_*
    for (int i = 0; i < (int)all_reads.size(); ++i)
    {
        const auto &rd = all_reads[i];
        for (auto key : rd.sorted_fwd)
            index_fwd[key].push_back(i);
        for (auto key : rd.sorted_rev)
            index_rev[key].push_back(i);
    }
}

// generateCandidatePairs ottimizzato
std::vector<Pair> generateCandidatePairs(
    const std::unordered_map<long long, std::vector<int>> &index_fwd,
    const std::unordered_map<long long, std::vector<int>> &index_rev)
{
    std::vector<Pair> pairs;
    // Accumula fwd–fwd
    for (const auto &kv : index_fwd)
    {
        const auto &v = kv.second;
        for (size_t a = 0; a < v.size(); ++a)
            for (size_t b = a + 1; b < v.size(); ++b)
                pairs.push_back({std::min(v[a], v[b]), std::max(v[a], v[b])});
    }
    // Accumula rev–rev
    for (const auto &kv : index_rev)
    {
        const auto &v = kv.second;
        for (size_t a = 0; a < v.size(); ++a)
            for (size_t b = a + 1; b < v.size(); ++b)
                pairs.push_back({std::min(v[a], v[b]), std::max(v[a], v[b])});
    }
    // Accumula cross fwd–rev
    for (const auto &kv : index_fwd)
    {
        auto it = index_rev.find(kv.first);
        if (it == index_rev.end()) continue;
        const auto &vf = kv.second;
        const auto &vr = it->second;
        for (int x : vf)
            for (int y : vr)
                if (x != y)
                    pairs.push_back({std::min(x, y), std::max(x, y)});
    }

    // Ordina e rimuovi duplicati
    std::sort(pairs.begin(), pairs.end(), [](const Pair &a, const Pair &b) {
        return (a.i < b.i) || (a.i == b.i && a.j < b.j);
    });
    pairs.erase(std::unique(pairs.begin(), pairs.end(),
                             [](const Pair &a, const Pair &b) {
                                 return a.i == b.i && a.j == b.j;
                             }),
                pairs.end());
    return pairs;
}

// identical to original
int hybrid_sorted_intersection_size(const std::vector<long long> &v1,
                                    const std::vector<long long> &v2)
{
    if (v1.size() < 20 || v2.size() < 20) {
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
    while (i < (int)v1.size() && j < (int)v2.size()) {
        if (v1[i] < v2[j]) ++i;
        else if (v2[j] < v1[i]) ++j;
        else { ++cnt; ++i; ++j; }
    }
    return cnt;
}