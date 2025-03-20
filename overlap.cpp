#include "overlap.hpp"
#include <vector>
#include <unordered_map>
#include <tuple>
#include <algorithm>
using namespace std;

static vector<long long> compute_power(int n, long long base, long long mod) {
    vector<long long> power(n+1, 1);
    for (int i = 1; i <= n; i++) {
        power[i] = (power[i-1] * base) % mod;
    }
    return power;
}

static inline long long get_hash(const vector<long long>& prefix, const vector<long long>& power, int l, int r, long long mod) {
    long long h = (prefix[r] - (prefix[l] * power[r - l]) % mod + mod) % mod;
    return h;
}

static tuple<bool, int, int> check_common_double(
    const vector<long long>& fp1, const vector<long long>& fp2,
    const vector<long long>& prefix1_mod1, const vector<long long>& prefix1_mod2,
    const vector<long long>& prefix2_mod1, const vector<long long>& prefix2_mod2,
    const vector<long long>& power_mod1, const vector<long long>& power_mod2,
    int L, long long mod1, long long mod2)
{
    unordered_map<unsigned long long, vector<int>> hash_map;
    int n1 = fp1.size(), n2 = fp2.size();
    for (int i = 0; i <= n1 - L; i++) {
        long long h1 = get_hash(prefix1_mod1, power_mod1, i, i+L, mod1);
        long long h2 = get_hash(prefix1_mod2, power_mod2, i, i+L, mod2);
        unsigned long long combined = ((unsigned long long)h1 << 32) ^ ((unsigned long long)h2);
        hash_map[combined].push_back(i);
    }
    for (int j = 0; j <= n2 - L; j++) {
        long long h1 = get_hash(prefix2_mod1, power_mod1, j, j+L, mod1);
        long long h2 = get_hash(prefix2_mod2, power_mod2, j, j+L, mod2);
        unsigned long long combined = ((unsigned long long)h1 << 32) ^ ((unsigned long long)h2);
        if (hash_map.find(combined) != hash_map.end()) {
            for (int i : hash_map[combined]) {
                bool match = true;
                for (int k = 0; k < L; k++) {
                    if (fp1[i+k] != fp2[j+k]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    return make_tuple(true, i, j);
                }
            }
        }
    }
    return make_tuple(false, -1, -1);
}

tuple<int, int, int, int, int, int, int> graph_overlap_fp_precomputed(
    const vector<long long>& comp_fp1,
    const vector<long long>& prefix1_mod1,
    const vector<long long>& prefix1_mod2,
    const vector<int>& idx1,
    const vector<long long>& comp_fp2,
    const vector<long long>& prefix2_mod1,
    const vector<long long>& prefix2_mod2,
    const vector<int>& idx2,
    int k)
{
    int n1 = comp_fp1.size();
    int n2 = comp_fp2.size();
    if(n1 == 0 || n2 == 0)
        return make_tuple(0, -1, -1, -1, -1, -1, -1);
    long long mod1 = 1000000007LL, mod2 = 1000000009LL, base = 131LL;
    int maxLen = min(n1, n2);
    vector<long long> power_mod1 = compute_power(maxLen, base, mod1);
    vector<long long> power_mod2 = compute_power(maxLen, base, mod2);
    
    int low = 0, high = maxLen + 1;
    int bestL = 0;
    int best_i = -1, best_j = -1;
    while(low < high - 1) {
        int mid = (low + high) / 2;
        bool found;
        int iCandidate, jCandidate;
        tie(found, iCandidate, jCandidate) = check_common_double(
            comp_fp1, comp_fp2,
            prefix1_mod1, prefix1_mod2,
            prefix2_mod1, prefix2_mod2,
            power_mod1, power_mod2,
            mid, mod1, mod2
        );
        if(found) {
            low = mid;
            bestL = mid;
            best_i = iCandidate;
            best_j = jCandidate;
        } else {
            high = mid;
        }
    }
    if(bestL == 0)
        return make_tuple(0, -1, -1, -1, -1, -1, -1);
    int orig_start1 = idx1[best_i];
    int orig_end1 = idx1[best_i + bestL - 1] + k;
    int orig_start2 = idx2[best_j];
    int orig_end2 = idx2[best_j + bestL - 1] + k;
    return make_tuple(bestL, orig_start1, orig_end1, orig_start2, orig_end2, best_i, best_j);
}
