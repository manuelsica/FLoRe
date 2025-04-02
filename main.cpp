#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <tuple>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <thread>
#include <mutex>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <iterator>
#include <cassert>
#include <sys/resource.h>
#include <queue>
#include <condition_variable>
#include <atomic>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#include "util.hpp"
#include "overlap.hpp"
#include "read.hpp"

using namespace std;

/*
 * =====================================================================================
 * Logging asincrono: definizione della coda, mutex, condizione e flag di terminazione
 * =====================================================================================
 */
mutex logQueueMutex;
condition_variable logQueueCV;
queue<string> logQueue;
atomic<bool> loggingDone(false);

void async_log(const string &msg)
{
    {
        lock_guard<mutex> lock(logQueueMutex);
        logQueue.push(msg);
    }
    logQueueCV.notify_one();
}

void loggingThreadFunction()
{
    while (true)
    {
        unique_lock<mutex> lock(logQueueMutex);
        logQueueCV.wait(lock, [] { return !logQueue.empty() || loggingDone.load(); });
        while (!logQueue.empty())
        {
            string msg = logQueue.front();
            logQueue.pop();
            lock.unlock();
            cout << msg << flush;
            lock.lock();
        }
        if (loggingDone.load() && logQueue.empty())
            break;
    }
}

/*
 * =====================================================================================
 * Strutture di supporto e dichiarazioni (come nell'implementazione originale)
 * =====================================================================================
 */
struct ReadData
{
    string forward;
    string reverse_seq;

    ProcessedRead pr_fwd;
    ProcessedRead pr_rev;

    unordered_set<long long> set_fwd;
    unordered_set<long long> set_rev;
    vector<long long> sorted_fwd;
    vector<long long> sorted_rev;

    bool sa_built_fwd = false;
    bool sa_built_rev = false;
    SuffixAutomaton suffix_automaton_fwd;
    SuffixAutomaton suffix_automaton_rev;
};

struct OverlapResult
{
    int overlap_len = 0;
    int start1 = 0;
    int end1 = 0;
    int start2 = 0;
    int end2 = 0;
    string orientation1;
    string orientation2;
    string combination;
    string r1;
    string r2;
    string fingerprint_r1;
    string fingerprint_r2;
    string used_algorithm;
};

struct Pair
{
    int i, j;
    bool operator==(const Pair &o) const { return (i == o.i && j == o.j); }
};
struct PairHash
{
    size_t operator()(const Pair &p) const
    {
        return ((size_t)p.i * 1315423917ULL) ^ (size_t)p.j;
    }
};

/*
 * Struttura per memorizzare i risultati JSON ordinati
 */
struct JsonResult
{
    int read1;
    int read2;
    string json;
};

/*
 * Forward declaration
 */
static OverlapResult compare_candidate_pair(ReadData &r1,
                                            ReadData &r2,
                                            int k,
                                            int min_overlap,
                                            bool verbose,
                                            int max_repeat_threshold);

/*
 * =====================================================================================
 * Funzioni di utilità
 * =====================================================================================
 */
// La funzione ora utilizza le API Mach su macOS per ottenere la memoria residente attuale.
static size_t getMemoryUsageKB()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc, sizeof(pmc)))
    {
        return pmc.PeakWorkingSetSize / 1024;
    }
    return 0;
#elif defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    kern_return_t kr = task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                                 (task_info_t)&info, &infoCount);
    if (kr != KERN_SUCCESS)
    {
        return 0;
    }
    return info.resident_size / 1024; // Memoria attuale in KB
#else
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // Su Linux questo restituisce il picco, eventualmente va adattato
#endif
}

static string getFingerprintRegion(const ProcessedRead &pr, int comp_start, int length)
{
    ostringstream oss;
    int stop = min<int>(comp_start + length, (int)pr.comp.comp_fp.size());
    for (int i = comp_start; i < stop; i++)
    {
        oss << pr.comp.comp_fp[i];
        if (i < stop - 1)
        {
            oss << "-";
        }
    }
    return oss.str();
}

/*
 * Se overlap < min_overlap => "SCARTATA"
 * Se troppi valori consecutivi => "SCARTATA - troppi valori..."
 * Altrimenti stringa vuota => "non scartato"
 */
static string get_overlap_annotation(const string &region,
                                     int fingerprint_match,
                                     int min_overlap,
                                     int max_repeat_threshold)
{
    if (fingerprint_match < min_overlap)
    {
        return "(SCARTATA)";
    }
    char current = '\0';
    int currentCount = 0, maxCount = 0;
    char maxChar = '\0';
    for (char c : region)
    {
        if (c == current)
            currentCount++;
        else
        {
            current = c;
            currentCount = 1;
        }
        if (currentCount > maxCount)
        {
            maxCount = currentCount;
            maxChar = c;
        }
    }
    if (maxCount > max_repeat_threshold)
    {
        ostringstream oss;
        oss << "(SCARTATA - troppi valori consecutivi di '" << maxChar << "')";
        return oss.str();
    }
    return "";
}

static string safe_substr(const string &s, size_t start, size_t length)
{
    if (start >= s.size())
        return "";
    size_t max_len = s.size() - start;
    if (length > max_len)
        length = max_len;
    return s.substr(start, length);
}

/*
 * =====================================================================================
 * Algoritmi Overlap
 * (I metodi FGOE, AOE, KHS e COIN rimangono invariati)
 * =====================================================================================
 */
// 1) FGOE
static tuple<int, int, int> fingerprint_guided_overlap_extension(const ProcessedRead &pr1,
                                                                 const ProcessedRead &pr2,
                                                                 int k,
                                                                 int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int size1 = (int)fp1.size();
    int size2 = (int)fp2.size();
    if (size1 == 0 || size2 == 0)
    {
        return {0, 0, 0};
    }
    int i = size1 - 1;
    int j = 0;
    int length_match = 0;
    while (i >= 0 && j < size2 && fp1[i] == fp2[j])
    {
        length_match++;
        i--;
        j++;
    }
    if (length_match >= min_overlap)
    {
        int startA = pr1.comp.comp_indices[i + 1];
        int startB = pr2.comp.comp_indices[0];
        return {length_match, startA, startB};
    }
    return {0, 0, 0};
}

// 2) AOE
static tuple<int, int, int> adaptive_overlap_extension(const ProcessedRead &pr1,
                                                       const ProcessedRead &pr2,
                                                       int k,
                                                       int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int n1 = (int)fp1.size();
    int n2 = (int)fp2.size();
    if (n1 == 0 || n2 == 0)
    {
        return {0, 0, 0};
    }
    int best_len = 0, best_x = 0, best_y = 0;
    for (int x = 0; x < n1; x++)
    {
        for (int y = 0; y < n2; y++)
        {
            if (fp1[x] == fp2[y])
            {
                int tmp_len = 1;
                int xx = x + 1, yy = y + 1;
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
                    {
                        break;
                    }
                }
            }
        }
        if (best_len >= min_overlap)
        {
            break;
        }
    }
    if (best_len >= min_overlap)
    {
        int startA = pr1.comp.comp_indices[best_x];
        int startB = pr2.comp.comp_indices[best_y];
        return {best_len, startA, startB};
    }
    return {0, 0, 0};
}

// 3) KHS
static tuple<int, int, int> kmer_hopping_search(const ProcessedRead &pr1,
                                                const ProcessedRead &pr2,
                                                int k,
                                                int min_overlap)
{
    const auto &fp1 = pr1.comp.comp_fp;
    const auto &fp2 = pr2.comp.comp_fp;
    int n1 = (int)fp1.size();
    int n2 = (int)fp2.size();
    if (n1 == 0 || n2 == 0)
    {
        return {0, 0, 0};
    }
    unordered_map<unsigned int, vector<int>> positionsInB;
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
            for (int posB : it->second)
            {
                int length_match = 1;
                int xx = i + 1, yy = posB + 1;
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
                    {
                        break;
                    }
                }
            }
        }
        if (best_len >= min_overlap)
        {
            break;
        }
    }
    if (best_len >= min_overlap)
    {
        int startA = pr1.comp.comp_indices[best_x];
        int startB = pr2.comp.comp_indices[best_y];
        return {best_len, startA, startB};
    }
    return {0, 0, 0};
}

// 4) COIN (Combined Overlap IntelliSearch)
static tuple<int, int, int> combined_overlap_intellisearch(const ProcessedRead &pr1,
                                                           const ProcessedRead &pr2,
                                                           int k,
                                                           int min_overlap)
{
    auto [fgoe_len, fgoe_startA, fgoe_startB] = fingerprint_guided_overlap_extension(pr1, pr2, k, 1);
    auto [aoe_len, aoe_startA, aoe_startB] = adaptive_overlap_extension(pr1, pr2, k, 1);
    auto [khs_len, khs_startA, khs_startB] = kmer_hopping_search(pr1, pr2, k, 1);

    if (fgoe_len == 0 && aoe_len == 0 && khs_len == 0)
    {
        return {0, 0, 0};
    }
    vector<tuple<int, int, int>> partials = {
        {fgoe_len, fgoe_startA, fgoe_startB},
        {aoe_len, aoe_startA, aoe_startB},
        {khs_len, khs_startA, khs_startB}};
    sort(partials.begin(), partials.end(),
         [](auto &a, auto &b)
         { return get<1>(a) < get<1>(b); });

    int best_len = 0, best_startA = 0, best_startB = 0;
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
                if (abs(sA2 - endA1) <= k && abs(sB2 - endB1) <= k)
                {
                    int combined_len = len1 + len2;
                    if (combined_len > best_len)
                    {
                        best_len = combined_len;
                        best_startA = min(sA1, sA2);
                        best_startB = min(sB1, sB2);
                    }
                }
            }
        }
    }
    int naive_best = 0, naive_A = 0, naive_B = 0;
    for (auto &p : partials)
    {
        if (get<0>(p) > naive_best)
        {
            naive_best = get<0>(p);
            naive_A = get<1>(p);
            naive_B = get<2>(p);
        }
    }
    if (best_len >= min_overlap && best_len > naive_best)
    {
        return {best_len, best_startA, best_startB};
    }
    else if (naive_best >= min_overlap)
    {
        return {naive_best, naive_A, naive_B};
    }
    return {0, 0, 0};
}

/*
 * =====================================================================================
 * OIIB, HCPG, HSIF, building
 * =====================================================================================
 */
void build_inverted_index_optimized(const vector<ReadData> &all_reads,
                                    unordered_map<long long, vector<int>> &index_fwd,
                                    unordered_map<long long, vector<int>> &index_rev)
{
    for (int i = 0; i < (int)all_reads.size(); i++)
    {
        {
            unordered_set<long long> uniq(all_reads[i].pr_fwd.comp.comp_fp.begin(),
                                          all_reads[i].pr_fwd.comp.comp_fp.end());
            for (auto val : uniq)
            {
                index_fwd[val].push_back(i);
            }
        }
        {
            unordered_set<long long> uniq(all_reads[i].pr_rev.comp.comp_fp.begin(),
                                          all_reads[i].pr_rev.comp.comp_fp.end());
            for (auto val : uniq)
            {
                index_rev[val].push_back(i);
            }
        }
    }
}

vector<Pair> generate_candidate_pairs_hash(const unordered_map<long long, vector<int>> &index_fwd,
                                           const unordered_map<long long, vector<int>> &index_rev)
{
    unordered_set<Pair, PairHash> pairs;
    for (const auto &kv : index_fwd)
    {
        const vector<int> &vec = kv.second;
        for (size_t i = 0; i < vec.size(); i++)
        {
            for (size_t j = i + 1; j < vec.size(); j++)
            {
                int a = vec[i], b = vec[j];
                if (a > b)
                    std::swap(a, b);
                pairs.insert({a, b});
            }
        }
    }
    for (const auto &kv : index_rev)
    {
        const vector<int> &vec = kv.second;
        for (size_t i = 0; i < vec.size(); i++)
        {
            for (size_t j = i + 1; j < vec.size(); j++)
            {
                int a = vec[i], b = vec[j];
                if (a > b)
                    std::swap(a, b);
                pairs.insert({a, b});
            }
        }
    }
    for (const auto &kv : index_fwd)
    {
        auto it = index_rev.find(kv.first);
        if (it != index_rev.end())
        {
            const vector<int> &vf = kv.second;
            const vector<int> &vr = it->second;
            for (int f : vf)
            {
                for (int r : vr)
                {
                    if (f == r)
                        continue;
                    int a = f, b = r;
                    if (a > b)
                        std::swap(a, b);
                    pairs.insert({a, b});
                }
            }
        }
    }
    return vector<Pair>(pairs.begin(), pairs.end());
}

int hybrid_sorted_intersection_size(const vector<long long> &v1, const vector<long long> &v2)
{
    if (v1.size() < 20 || v2.size() < 20)
    {
        int cnt = 0;
        if (v1.size() < v2.size())
        {
            for (auto &x : v1)
            {
                if (binary_search(v2.begin(), v2.end(), x))
                {
                    cnt++;
                }
            }
        }
        else
        {
            for (auto &x : v2)
            {
                if (binary_search(v1.begin(), v1.end(), x))
                {
                    cnt++;
                }
            }
        }
        return cnt;
    }
    int i = 0, j = 0, cnt = 0;
    while (i < (int)v1.size() && j < (int)v2.size())
    {
        if (v1[i] < v2[j])
        {
            i++;
        }
        else if (v2[j] < v1[i])
        {
            j++;
        }
        else
        {
            cnt++;
            i++;
            j++;
        }
    }
    return cnt;
}

static void build_read_data(vector<ReadData> &all_reads, const vector<string> &reads, int k)
{
    for (size_t i = 0; i < all_reads.size(); i++)
    {
        all_reads[i].forward = reads[i];
        all_reads[i].reverse_seq = reverse_complement(reads[i]);
        all_reads[i].pr_fwd = process_read_optimized(all_reads[i].forward, k);
        all_reads[i].pr_rev = process_read_optimized(all_reads[i].reverse_seq, k);
    }
}

static void fill_fingerprint_sets(vector<ReadData> &all_reads)
{
    for (auto &rd : all_reads)
    {
        for (auto val : rd.pr_fwd.comp.comp_fp)
        {
            rd.set_fwd.insert((long long)val);
        }
        for (auto val : rd.pr_rev.comp.comp_fp)
        {
            rd.set_rev.insert((long long)val);
        }
        rd.sorted_fwd.assign(rd.set_fwd.begin(), rd.set_fwd.end());
        rd.sorted_rev.assign(rd.set_rev.begin(), rd.set_rev.end());
        sort(rd.sorted_fwd.begin(), rd.sorted_fwd.end());
        sort(rd.sorted_rev.begin(), rd.sorted_rev.end());
    }
}

static void build_inverted_index(const vector<ReadData> &all_reads,
                                 unordered_map<long long, vector<int>> &index_fwd,
                                 unordered_map<long long, vector<int>> &index_rev)
{
    build_inverted_index_optimized(all_reads, index_fwd, index_rev);
}

static vector<Pair> generate_candidate_pairs_vector(
    const unordered_map<long long, vector<int>> &index_fwd,
    const unordered_map<long long, vector<int>> &index_rev)
{
    return generate_candidate_pairs_hash(index_fwd, index_rev);
}

/*
 * =====================================================================================
 * compare_candidate_pair
 * (Rimane invariata la logica, come nelle versioni precedenti)
 * =====================================================================================
 */
static OverlapResult compare_candidate_pair(ReadData &r1,
                                            ReadData &r2,
                                            int k,
                                            int min_overlap,
                                            bool verbose,
                                            int max_repeat_threshold)
{
    OverlapResult best;

    auto try_update = [&](int len, int comp_idx1, int comp_idx2,
                          const ProcessedRead &prA, const ProcessedRead &prB,
                          const string &readA, const string &readB,
                          const string &orientA, const string &orientB,
                          const string &comb, const string &algorithmName)
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
            best.fingerprint_r1 = getFingerprintRegion(prA, comp_idx1, len);
            best.fingerprint_r2 = getFingerprintRegion(prB, comp_idx2, len);
            best.used_algorithm = algorithmName;
        }
    };

    auto compute_best_overlap_4stages = [&](ProcessedRead &A, ProcessedRead &B,
                                            string &rA, string &rB,
                                            const string &oA, const string &oB,
                                            const string &comb)
    {
        // 1) FGOE
        {
            auto [matchLen, startA, startB] = fingerprint_guided_overlap_extension(A, B, k, min_overlap);
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
                try_update(matchLen, idxA, idxB, A, B, rA, rB, oA, oB, comb, "FGOE");
                return;
            }
        }
        // 2) AOE
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
        // 3) KHS
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
        // 4) COIN
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

    auto compute_best_overlap_suffix_automaton = [&](ProcessedRead &A, ProcessedRead &B,
                                                     string &rA, string &rB,
                                                     const string &oA, const string &oB,
                                                     const string &comb,
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
            best.fingerprint_r1 = getFingerprintRegion(A, c1, length);
            best.fingerprint_r2 = getFingerprintRegion(B, c2, length);
            best.used_algorithm = "SuffixAutomaton";
        }
    };

    compute_best_overlap_4stages(r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward, "forward", "forward", "ff");
    compute_best_overlap_4stages(r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq, "forward", "reverse", "fr");
    compute_best_overlap_4stages(r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward, "reverse", "forward", "rf");
    compute_best_overlap_4stages(r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq, "reverse", "reverse", "rr");

    if (best.overlap_len < min_overlap)
    {
        compute_best_overlap_suffix_automaton(r1.pr_fwd, r2.pr_fwd, r1.forward, r2.forward,
                                              "forward", "forward", "ff",
                                              r1.sa_built_fwd, r1.suffix_automaton_fwd);
        compute_best_overlap_suffix_automaton(r1.pr_fwd, r2.pr_rev, r1.forward, r2.reverse_seq,
                                              "forward", "reverse", "fr",
                                              r1.sa_built_fwd, r1.suffix_automaton_fwd);
        compute_best_overlap_suffix_automaton(r1.pr_rev, r2.pr_fwd, r1.reverse_seq, r2.forward,
                                              "reverse", "forward", "rf",
                                              r1.sa_built_rev, r1.suffix_automaton_rev);
        compute_best_overlap_suffix_automaton(r1.pr_rev, r2.pr_rev, r1.reverse_seq, r2.reverse_seq,
                                              "reverse", "reverse", "rr",
                                              r1.sa_built_rev, r1.suffix_automaton_rev);
    }
    return best;
}

/*
 * =====================================================================================
 * Funzione per scrivere il JSON ordinato
 * =====================================================================================
 */
void write_sorted_json(vector<JsonResult>& json_results, const string& filename) {
    // Ordina il vettore in base ai campi "read1" e "read2"
    sort(json_results.begin(), json_results.end(), [](const JsonResult &a, const JsonResult &b) {
        if (a.read1 != b.read1)
            return a.read1 < b.read1;
        return a.read2 < b.read2;
    });

    // Scrive il file JSON ordinato
    ofstream jf(filename);
    if (!jf.is_open()) {
        cerr << "Errore apertura file " << filename << "\n";
        return;
    }
    jf << "[\n";
    for (size_t i = 0; i < json_results.size(); i++) {
        jf << json_results[i].json;
        if (i < json_results.size() - 1) {
            jf << ",\n";
        }
    }
    jf << "\n]\n";
    jf.close();
}

/*
 * =====================================================================================
 * print_usage
 * =====================================================================================
 */
static void print_usage(const char *prog_name)
{
    cerr << "Usage: " << prog_name << " -f <file_fasta> [options]\n"
         << "  -f, --fasta <file>\n"
         << "  -m, --min_overlap <int>          [default: 13]\n"
         << "  -r, --max_repeat_threshold <int> [default: 10]\n"
         << "  -k, --kmer <int>                 [default: 15]\n"
         << "  -j, --json <file>                [default: results.json]\n"
         << "  -t, --threads <int>              [default: hardware_concurrency]\n"
         << "  -v, --verbose\n"
         << "  -h, --help\n";
}

/*
 * =====================================================================================
 * main
 * =====================================================================================
 */
int main(int argc, char *argv[])
{
    // Avvio del thread di logging asincrono
    thread loggerThread(loggingThreadFunction);

    string fasta_file;
    int min_overlap = 13;
    int max_repeat_threshold = 10;
    int k = 15;
    string json_filename = "results.json";
    unsigned int num_threads = 0;
    bool verbose = false;

    struct option longopts[] = {
        {"fasta", required_argument, 0, 'f'},
        {"min_overlap", required_argument, 0, 'm'},
        {"max_repeat_threshold", required_argument, 0, 'r'},
        {"kmer", required_argument, 0, 'k'},
        {"json", required_argument, 0, 'j'},
        {"threads", required_argument, 0, 't'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};
    int opt;
    while ((opt = getopt_long(argc, argv, "f:m:r:k:j:t:vh", longopts, NULL)) != -1)
    {
        switch (opt)
        {
        case 'f':
            fasta_file = optarg;
            break;
        case 'm':
            min_overlap = atoi(optarg);
            break;
        case 'r':
            max_repeat_threshold = atoi(optarg);
            break;
        case 'k':
            k = atoi(optarg);
            break;
        case 'j':
            json_filename = optarg;
            break;
        case 't':
            num_threads = (unsigned int)atoi(optarg);
            break;
        case 'v':
            verbose = true;
            break;
        case 'h':
            print_usage(argv[0]);
            return 0;
        default:
            print_usage(argv[0]);
            return 1;
        }
    }
    if (fasta_file.empty())
    {
        cerr << "Errore: specificare il file FASTA\n";
        print_usage(argv[0]);
        return 1;
    }
    if (num_threads == 0)
    {
        num_threads = thread::hardware_concurrency();
        if (num_threads == 0)
        {
            num_threads = 2;
        }
    }

    // Operazioni iniziali (non incluse nella misurazione)
    async_log("Lettura del file FASTA: " + fasta_file + "\n");
    vector<string> reads = read_fasta_buffered(fasta_file);
    if (reads.empty())
    {
        cerr << "Nessuna read trovata.\n";
        return 1;
    }
    async_log("Numero di read: " + to_string(reads.size()) + "\n");
    async_log("k-mer: " + to_string(k) + "\n");

    vector<ReadData> all_reads(reads.size());
    build_read_data(all_reads, reads, k);
    fill_fingerprint_sets(all_reads);

    // Svuota il vettore delle reads, non più necessario (risparmia memoria)
    reads.clear();
    reads.shrink_to_fit();

    async_log("Costruzione indici invertiti (OIIB)...\n");
    unordered_map<long long, vector<int>> index_fwd, index_rev;
    index_fwd.reserve(all_reads.size() * 4);
    index_rev.reserve(all_reads.size() * 4);
    build_inverted_index(all_reads, index_fwd, index_rev);

    // Dopo la generazione delle coppie, liberiamo le mappe per ridurre l'uso di memoria
    async_log("Generazione coppie candidate (HCPG)...\n");
    vector<Pair> candidate_pairs = generate_candidate_pairs_vector(index_fwd, index_rev);
    index_fwd.clear();
    index_rev.clear();
    unordered_map<long long, vector<int>>().swap(index_fwd);
    unordered_map<long long, vector<int>>().swap(index_rev);

    async_log("Pre-filtraggio delle coppie candidate (HSIF)...\n");
    vector<Pair> filtered;
    filtered.reserve(candidate_pairs.size());
    for (auto &c : candidate_pairs)
    {
        auto &rd1 = all_reads[c.i];
        auto &rd2 = all_reads[c.j];
        int ff = hybrid_sorted_intersection_size(rd1.sorted_fwd, rd2.sorted_fwd);
        int fr = hybrid_sorted_intersection_size(rd1.sorted_fwd, rd2.sorted_rev);
        int rf = hybrid_sorted_intersection_size(rd1.sorted_rev, rd2.sorted_fwd);
        int rr = hybrid_sorted_intersection_size(rd1.sorted_rev, rd2.sorted_rev);
        int best_val = max({ff, fr, rf, rr});
        if (best_val >= min_overlap)
        {
            filtered.push_back(c);
        }
    }
    // Liberiamo candidate_pairs, non più necessario
    candidate_pairs.clear();
    candidate_pairs.shrink_to_fit();
    filtered.shrink_to_fit();

    async_log("Coppie candidate dopo filtro: " + to_string(filtered.size()) + "\n");

    // Avvio della misurazione del tempo (I/O iniziali esclusi)
    using Clock = chrono::steady_clock;
    auto start_time = Clock::now();

    async_log("Elaborazione candidate in multi-thread (" + to_string(num_threads) + " thread)...\n");

    mutex json_mutex;
    vector<JsonResult> json_results;
    json_results.reserve(filtered.size());
    atomic<size_t> candidateIndex(0);

    auto worker = [&]()
    {
        while (true)
        {
            size_t idx = candidateIndex.fetch_add(1);
            if (idx >= filtered.size())
                break;
            auto &p = filtered[idx];
            OverlapResult best_ov = compare_candidate_pair(all_reads[p.i],
                                                           all_reads[p.j],
                                                           k,
                                                           min_overlap,
                                                           verbose,
                                                           max_repeat_threshold);
            if (best_ov.overlap_len > 0)
            {
                string region_r1 = safe_substr(best_ov.r1, best_ov.start1, best_ov.end1 - best_ov.start1);
                string region_r2 = safe_substr(best_ov.r2, best_ov.start2, best_ov.end2 - best_ov.start2);
                string annotation = get_overlap_annotation(region_r1, best_ov.overlap_len, min_overlap, max_repeat_threshold);
                if (verbose)
                {
                    async_log("-----------------------------------\n");
                    async_log("Coppia di read " + to_string(p.i + 1) + " - " + to_string(p.j + 1) + "\n");
                    async_log("Overlap = " + to_string(best_ov.overlap_len) + " " + annotation + "\n");
                    async_log("Algoritmo usato: " + best_ov.used_algorithm + "\n");
                    async_log("Orientation: " + best_ov.orientation1 + " - " + best_ov.orientation2 + "\n");
                    async_log("Region r1: " + region_r1 + "\n");
                    async_log("Fingerprint r1: " + best_ov.fingerprint_r1 + "\n");
                    async_log("Region r2: " + region_r2 + "\n");
                    async_log("Fingerprint r2: " + best_ov.fingerprint_r2 + "\n");
                }
                if (annotation.find("SCARTATA") == string::npos)
                {
                    ostringstream oss;
                    oss << "{";
                    oss << "\"read1\":" << (p.i + 1) << ","
                        << "\"read2\":" << (p.j + 1) << ","
                        << "\"orientation1\":\"" << best_ov.orientation1 << "\","
                        << "\"orientation2\":\"" << best_ov.orientation2 << "\","
                        << "\"start1\":" << best_ov.start1 << ",\"end1\":" << best_ov.end1 << ","
                        << "\"len_read1\":" << best_ov.r1.size() << ","
                        << "\"start2\":" << best_ov.start2 << ",\"end2\":" << best_ov.end2 << ","
                        << "\"len_read2\":" << best_ov.r2.size() << ","
                        << "\"overlap_length\":" << best_ov.overlap_len << ","
                        << "\"overlap_region_read1\":\"" << region_r1 << "\","
                        << "\"fingerprint_read1\":\"" << best_ov.fingerprint_r1 << "\","
                        << "\"overlap_region_read2\":\"" << region_r2 << "\","
                        << "\"fingerprint_read2\":\"" << best_ov.fingerprint_r2 << "\","
                        << "\"used_algorithm\":\"" << best_ov.used_algorithm << "\"";
                    oss << "}";
                    {
                        lock_guard<mutex> lk(json_mutex);
                        json_results.push_back(JsonResult{p.i + 1, p.j + 1, oss.str()});
                    }
                }
            }
        }
    };

    vector<thread> pool;
    pool.reserve(num_threads);
    for (unsigned int th = 0; th < num_threads; th++)
    {
        pool.emplace_back(worker);
    }
    for (auto &t : pool)
    {
        t.join();
    }
    // Liberiamo il vettore delle coppie filtrate
    filtered.clear();
    filtered.shrink_to_fit();

    // Fine elaborazione: calcolo tempo di esecuzione e memoria usata
    auto end_time = Clock::now();
    double seconds = chrono::duration<double>(end_time - start_time).count();
    size_t mem_kb = getMemoryUsageKB();

    async_log("Salvataggio risultati in " + json_filename + "...\n");
    write_sorted_json(json_results, json_filename);

    if (!verbose)
    {
        async_log("\nTotale overlap scritti: " + to_string(json_results.size()) + "\n");
    }
    async_log("\nTempo di esecuzione: " + to_string(seconds) + " s\n");
    async_log("Memoria usata (attuale): " + to_string(mem_kb) + " KB\n");

    // Segnala al thread di logging di terminare e attendi la sua terminazione
    loggingDone.store(true);
    logQueueCV.notify_one();
    loggerThread.join();

    // Infine, se non è più necessario, liberiamo anche all_reads
    all_reads.clear();
    all_reads.shrink_to_fit();

    return 0;
}
