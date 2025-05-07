#ifndef INDEX_HPP
#define INDEX_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "util.hpp"
#include "overlap.hpp"

struct Pair {
    int i, j;
    bool operator==(const Pair &o) const noexcept { return i == o.i && j == o.j; }
};

struct PairHash {
    size_t operator()(const Pair &p) const noexcept {
        constexpr uint64_t phi = 0x9E3779B97F4A7C15ULL;                 // mix ϕ
        uint64_t v = (static_cast<uint64_t>(p.i) << 32) | uint32_t(p.j);
        v ^= (v >> 30); v *= phi; v ^= (v >> 27); v *= phi; v ^= (v >> 31);
        return static_cast<size_t>(v);
    }
};

struct ReadData {
    unsigned int first_fp_fwd = 0;
    unsigned int first_fp_rev = 0;
    std::string forward, reverse_seq;
    ProcessedRead pr_fwd, pr_rev;
    std::unordered_set<long long> set_fwd, set_rev;
    std::vector<long long> sorted_fwd, sorted_rev;
    SuffixAutomaton suffix_automaton_fwd, suffix_automaton_rev;
    bool sa_built_fwd{false}, sa_built_rev{false};
};

// ↓ funzioni – firme invariate
void buildAllReadsData(std::vector<ReadData> &all_reads,
    const std::vector<std::string> &reads,
    int k,
    bool use_solid_fp,
    const std::unordered_set<unsigned int> &solid_set,
    bool use_cfl,
    int  cfl_thr,
    bool use_icfl,
    int  icfl_thr);


void fillFingerprintSets(std::vector<ReadData>&);

void buildInvertedIndex(const std::vector<ReadData>&,
                        std::unordered_map<long long, std::vector<int>>&,
                        std::unordered_map<long long, std::vector<int>>&);

std::vector<Pair> generateCandidatePairs(
        const std::unordered_map<long long, std::vector<int>>&,
        const std::unordered_map<long long, std::vector<int>>&);

int hybrid_sorted_intersection_size(const std::vector<long long>&,
                                    const std::vector<long long>&);

#endif
