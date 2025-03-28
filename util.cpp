#include "util.hpp"
#include <vector>
#include <string>
#include <string_view>
#include <algorithm>

// Codifica un k‑mer usando 2 bit per base.
unsigned int encode_kmer_bit(const string &kmer) {
    unsigned int value = 0;
    for (char c : kmer) {
        value <<= 2;
        switch(c) {
            case 'A': case 'a': value |= 0; break;
            case 'C': case 'c': value |= 1; break;
            case 'G': case 'g': value |= 2; break;
            case 'T': case 't': value |= 3; break;
            default: value |= 0; break;
        }
    }
    return value;
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
            default:  result.push_back('N'); break;
        }
    }
    return result;
}

CompressedFingerprint compress_fingerprint(const vector<unsigned int> &fingerprint, const vector<string_view> &kmers) {
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
    int n = read.size();
    if (n < k)
        return pr;
    int numKmers = n - k + 1;
    pr.fingerprint.resize(numKmers);
    pr.kmers.reserve(numKmers);
    
    // Calcola il primo k‑mer usando la codifica bit‑level.
    string_view first_kmer(read.data(), k);
    unsigned int val = encode_kmer_bit(string(first_kmer));
    pr.fingerprint[0] = val;
    pr.kmers.push_back(first_kmer);
    
    int total_bits = 2 * k;
    unsigned int mask = (1u << total_bits) - 1;
    for (int i = 1; i < numKmers; i++) {
        // Aggiornamento in O(1): shift a sinistra, mascheratura e aggiunta della nuova base.
        val = (val << 2) & mask;
        char c = read[i + k - 1];
        unsigned int bits = 0;
        switch(c) {
            case 'A': case 'a': bits = 0; break;
            case 'C': case 'c': bits = 1; break;
            case 'G': case 'g': bits = 2; break;
            case 'T': case 't': bits = 3; break;
            default: bits = 0; break;
        }
        val |= bits;
        pr.fingerprint[i] = val;
        pr.kmers.push_back(string_view(read.data() + i, k));
    }
    
    pr.comp = compress_fingerprint(pr.fingerprint, pr.kmers);
    return pr;
}
