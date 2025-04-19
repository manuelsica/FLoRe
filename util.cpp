// util.cpp
// Implementa le funzioni utilitarie per il tool FLoRe.
#include "util.hpp"
#include <algorithm>
#include <cctype>
#include <unordered_map>
#include <sstream>

// encode_kmer_bit: Converte un k-mer in un intero utilizzando 2 bit per base.
unsigned int encode_kmer_bit(const std::string &kmer)
{
    unsigned int value = 0;
    for (char c : kmer)
    {
        value <<= 2;  // Sposta a sinistra di 2 bit per fare spazio alla nuova base
        switch(std::toupper(static_cast<unsigned char>(c)))
        {
            case 'A': break;              // A -> 00
            case 'C': value |= 1; break;    // C -> 01
            case 'G': value |= 2; break;    // G -> 10
            case 'T': value |= 3; break;    // T -> 11
            default: break;
        }
    }
    return value;
}

// reverse_complement: Calcola il reverse complement di una sequenza di DNA.
std::string reverse_complement(const std::string &seq)
{
    std::string result;
    result.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
    {
        switch(std::toupper(static_cast<unsigned char>(*it)))
        {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default:  result.push_back('N'); break;
        }
    }
    return result;
}

// compress_fingerprint: Elimina i k-mer consecutivi duplicati per ottenere il fingerprint compresso.
CompressedFingerprint compress_fingerprint(const std::vector<unsigned int> &fingerprint,
                                           const std::vector<std::string_view> &kmers)
{
    CompressedFingerprint comp;
    if (fingerprint.empty())
        return comp;
    comp.comp_fp.reserve(fingerprint.size());
    comp.comp_kmers.reserve(kmers.size());
    comp.comp_indices.reserve(kmers.size());
    comp.comp_fp.push_back(fingerprint[0]);
    comp.comp_kmers.push_back(kmers[0]);
    comp.comp_indices.push_back(0);
    for (size_t i = 1; i < fingerprint.size(); i++)
    {
        if (fingerprint[i] != fingerprint[i-1])
        {
            comp.comp_fp.push_back(fingerprint[i]);
            comp.comp_kmers.push_back(kmers[i]);
            comp.comp_indices.push_back((int)i);
        }
    }
    return comp;
}

// processReadFingerprint: Calcola il fingerprint per una read e lo comprime in una singola passata.
ProcessedRead processReadFingerprint(const std::string &read, int k)
{
    ProcessedRead pr;
    int n = (int)read.size();
    if (n < k)
        return pr;
    int numKmers = n - k + 1;
    pr.fingerprint.reserve(numKmers);
    pr.kmers.reserve(numKmers);
    unsigned int mask = (1u << (2 * k)) - 1;
    unsigned int val = 0;
    for (int i = 0; i < k; i++)
    {
        val <<= 2;
        char c = read[i];
        switch(std::toupper(static_cast<unsigned char>(c)))
        {
            case 'A': break;
            case 'C': val |= 1; break;
            case 'G': val |= 2; break;
            case 'T': val |= 3; break;
            default: break;
        }
    }
    pr.fingerprint.push_back(val);
    pr.kmers.push_back(std::string_view(read.data(), k));
    
    std::vector<unsigned int> comp_fp;
    std::vector<std::string_view> comp_kmers;
    std::vector<int> comp_indices;
    
    comp_fp.push_back(val);
    comp_kmers.push_back(std::string_view(read.data(), k));
    comp_indices.push_back(0);
    
    for (int i = 1; i < numKmers; i++)
    {
        val = ((val << 2) & mask);
        char c = read[i + k - 1];
        switch(std::toupper(static_cast<unsigned char>(c)))
        {
            case 'A': break;
            case 'C': val |= 1; break;
            case 'G': val |= 2; break;
            case 'T': val |= 3; break;
            default: break;
        }
        pr.fingerprint.push_back(val);
        pr.kmers.push_back(std::string_view(read.data() + i, k));
        if (val != comp_fp.back())
        {
            comp_fp.push_back(val);
            comp_kmers.push_back(std::string_view(read.data() + i, k));
            comp_indices.push_back(i);
        }
    }
    pr.comp.comp_fp = std::move(comp_fp);
    pr.comp.comp_kmers = std::move(comp_kmers);
    pr.comp.comp_indices = std::move(comp_indices);
    
    return pr;
}

// buildGlobalFingerprintFrequency: Conta la frequenza di ogni fingerprint in tutte le read.
std::unordered_map<unsigned int,int> buildGlobalFingerprintFrequency(const std::vector<std::string> &reads, int k)
{
    std::unordered_map<unsigned int,int> freq_map;
    for (auto &read : reads)
    {
        int n = (int)read.size();
        if (n < k) continue;
        unsigned int mask = (1u << (2 * k)) - 1;
        unsigned int val = 0;
        for (int i = 0; i < k; i++)
        {
            val <<= 2;
            char c = read[i];
            switch(std::toupper(static_cast<unsigned char>(c)))
            {
                case 'A': break;
                case 'C': val |= 1; break;
                case 'G': val |= 2; break;
                case 'T': val |= 3; break;
                default: break;
            }
        }
        freq_map[val]++;
        int numKmers = n - k + 1;
        for (int i = 1; i < numKmers; i++)
        {
            val = ((val << 2) & mask);
            char c = read[i + k - 1];
            switch(std::toupper(static_cast<unsigned char>(c)))
            {
                case 'A': break;
                case 'C': val |= 1; break;
                case 'G': val |= 2; break;
                case 'T': val |= 3; break;
                default: break;
            }
            freq_map[val]++;
        }
    }
    return freq_map;
}

// processReadSolidFingerprint: Processa una read scartando i fingerprint non "solidi".
ProcessedRead processReadSolidFingerprint(const std::string &read,
                                          int k,
                                          const std::unordered_set<unsigned int> &solid_fingerprint_set)
{
    ProcessedRead pr;
    int n = (int)read.size();
    if (n < k)
        return pr;
    int numKmers = n - k + 1;
    pr.fingerprint.reserve(numKmers);
    pr.kmers.reserve(numKmers);
    unsigned int mask = (1u << (2 * k)) - 1;
    unsigned int val = 0;
    for (int i = 0; i < k; i++)
    {
        val <<= 2;
        char c = read[i];
        switch(std::toupper(static_cast<unsigned char>(c)))
        {
            case 'A': break;
            case 'C': val |= 1; break;
            case 'G': val |= 2; break;
            case 'T': val |= 3; break;
            default: break;
        }
    }
    std::vector<unsigned int> comp_fp;
    std::vector<std::string_view> comp_kmers;
    std::vector<int> comp_indices;
    if (solid_fingerprint_set.find(val) != solid_fingerprint_set.end())
    {
        pr.fingerprint.push_back(val);
        pr.kmers.push_back(std::string_view(read.data(), k));
        comp_fp.push_back(val);
        comp_kmers.push_back(std::string_view(read.data(), k));
        comp_indices.push_back(0);
    }
    for (int i = 1; i < numKmers; i++)
    {
        val = ((val << 2) & mask);
        char c = read[i + k - 1];
        switch(std::toupper(static_cast<unsigned char>(c)))
        {
            case 'A': break;
            case 'C': val |= 1; break;
            case 'G': val |= 2; break;
            case 'T': val |= 3; break;
            default: break;
        }
        if (solid_fingerprint_set.find(val) != solid_fingerprint_set.end())
        {
            pr.fingerprint.push_back(val);
            pr.kmers.push_back(std::string_view(read.data() + i, k));
            if (comp_fp.empty() || (val != comp_fp.back()))
            {
                comp_fp.push_back(val);
                comp_kmers.push_back(std::string_view(read.data() + i, k));
                comp_indices.push_back(i);
            }
        }
    }
    pr.comp.comp_fp = std::move(comp_fp);
    pr.comp.comp_kmers = std::move(comp_kmers);
    pr.comp.comp_indices = std::move(comp_indices);
    return pr;
}

// safe_substr: Restituisce una sottostringa in maniera sicura.
std::string safe_substr(const std::string &s, size_t start, size_t length)
{
    if (start >= s.size())
        return "";
    size_t max_len = s.size() - start;
    if (length > max_len)
        length = max_len;
    return s.substr(start, length);
}

// get_overlap_annotation: Genera un'annotazione per l'overlap.
std::string get_overlap_annotation(const std::string &region,
                                   int fingerprint_match,
                                   int min_overlap,
                                   int max_repeat_threshold)
{
    if (fingerprint_match < min_overlap)
        return "(SCARTATA)";
    char current = '\0';
    int currentCount = 0, maxCount = 0;
    char maxChar = '\0';
    for (char c : region)
    {
        if (c == current)
            currentCount++;
        else {
            current = c;
            currentCount = 1;
        }
        if (currentCount > maxCount) {
            maxCount = currentCount;
            maxChar = c;
        }
    }
    if (maxCount > max_repeat_threshold)
    {
        std::ostringstream oss;
        oss << "(SCARTATA - troppi valori consecutivi di '" << maxChar << "')";
        return oss.str();
    }
    return "";
}

// buildSolidFingerprintSet: Costruisce un set di fingerprint "solidi" basato sulla frequenza minima e massima.
std::unordered_set<unsigned int> buildSolidFingerprintSet(const std::vector<std::string> &reads,
                                                             int k,
                                                             int min_freq,
                                                             int max_freq)
{
    std::unordered_set<unsigned int> solid_set;
    auto freq_map = buildGlobalFingerprintFrequency(reads, k);
    for (auto &kv : freq_map)
    {
        if (kv.second >= min_freq && kv.second <= max_freq)
            solid_set.insert(kv.first);
    }
    return solid_set;
}