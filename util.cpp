// util.cpp
#include "util.hpp"
#include <cmath>

long long compute_val(const string &factor, int base) {
    int L = factor.size();
    long long val = 0;
    for (int i = 0; i < L; i++) {
        int m;
        switch(factor[i]) {
            case 'A': m = 0; break;
            case 'C': m = 1; break;
            case 'G': m = 2; break;
            case 'T': m = 3; break;
            default: m = 0; break;
        }
        long long power = 1;
        for (int j = 0; j < L - 1 - i; j++) {
            power *= base;
        }
        val += m * power;
    }
    return val;
}

long long encode_factor(const string &factor, int shift, int base) {
    int L = factor.size();
    long long val = compute_val(factor, base);
    return (val << shift) | L;
}

vector<string> sliding_window_read(const string &read, int k, int step) {
    vector<string> segments;
    int n = read.size() - k + 1;
    segments.reserve(n);
    for (size_t i = 0; i + k <= read.size(); i += step) {
        segments.push_back(read.substr(i, k));
    }
    return segments;
}

string reverse_complement(const string &seq) {
    string result;
    result.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch(*it) {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default: result.push_back('N'); break;
        }
    }
    return result;
}

CompressedFingerprint compress_fingerprint(const vector<long long> &fingerprint, const vector<string> &kmers) {
    CompressedFingerprint comp;
    if (fingerprint.empty())
        return comp;
    comp.comp_fp.reserve(fingerprint.size());
    comp.comp_kmers.reserve(kmers.size());
    comp.comp_indices.reserve(kmers.size());
    comp.comp_fp.push_back(fingerprint[0]);
    comp.comp_kmers.push_back(kmers[0]);
    comp.comp_indices.push_back(0);
    for (size_t i = 1; i < fingerprint.size(); i++) {
        if (fingerprint[i] != fingerprint[i - 1]) {
            comp.comp_fp.push_back(fingerprint[i]);
            comp.comp_kmers.push_back(kmers[i]);
            comp.comp_indices.push_back(i);
        }
    }
    return comp;
}

ProcessedRead process_read(const string &read, int k) {
    ProcessedRead pr;
    int n = read.size() - k + 1;
    pr.kmers.reserve(n);
    pr.fingerprint.reserve(n);
    pr.kmers = sliding_window_read(read, k, 1);
    pr.fingerprint.resize(pr.kmers.size());
    for (size_t i = 0; i < pr.kmers.size(); i++) {
        pr.fingerprint[i] = encode_factor(pr.kmers[i], 8, 5);
    }
    pr.comp = compress_fingerprint(pr.fingerprint, pr.kmers);
    return pr;
}
