// util.cpp
// Implementa le funzioni utilitarie per il tool FLoRe.

#include "util.hpp"
#include <algorithm>    // per std::min, std::max, std::transform
#include <cctype>       // per std::toupper
#include <unordered_map>
#include <sstream>      // per std::ostringstream

/*
 * Funzione encode_kmer_bit:
 * Per ogni carattere del k-mer, sposta i bit verso sinistra e aggiunge il valore codificato.
 * A -> 00, C -> 01, G -> 10, T -> 11.
 */
unsigned int encode_kmer_bit(const std::string &kmer)
{
    unsigned int value = 0;
    for (char c : kmer)
    {
        value <<= 2;  // Sposta di 2 bit a sinistra
        switch(std::toupper(static_cast<unsigned char>(c)))
        {
            case 'A': break;              // A viene codificato come 00
            case 'C': value |= 1; break;    // C: 01
            case 'G': value |= 2; break;    // G: 10
            case 'T': value |= 3; break;    // T: 11
            default: break;
        }
    }
    return value;
}

/*
 * Funzione reverse_complement:
 * Calcola il complementare inverso di una sequenza di DNA.
 * Legge la stringa al contrario e sostituisce ogni base con il suo complementare.
 */
std::string reverse_complement(const std::string &seq)
{
    std::string result;
    result.reserve(seq.size()); // Prealloca lo spazio necessario
    // Scorre la sequenza in ordine inverso
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
    {
        // Converte la base in maiuscolo e sostituisce con il complementare
        switch(std::toupper(static_cast<unsigned char>(*it)))
        {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default:  result.push_back('N'); break;  // N per basi sconosciute
        }
    }
    return result;
}

/*
 * Funzione compress_fingerprint:
 * Crea una versione compressa del fingerprint eliminando k-mer consecutivi uguali.
 */
CompressedFingerprint compress_fingerprint(const std::vector<unsigned int> &fingerprint,
                                           const std::vector<std::string_view> &kmers)
{
    CompressedFingerprint comp;
    if (fingerprint.empty())
        return comp;
    // Prealloca lo spazio per i vettori compressi
    comp.comp_fp.reserve(fingerprint.size());
    comp.comp_kmers.reserve(kmers.size());
    comp.comp_indices.reserve(kmers.size());

    // Aggiunge il primo elemento
    comp.comp_fp.push_back(fingerprint[0]);
    comp.comp_kmers.push_back(kmers[0]);
    comp.comp_indices.push_back(0);

    // Scorre il fingerprint e aggiunge solo se diverso dal precedente
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

/*
 * Funzione processReadFingerprint:
 * Calcola il fingerprint per una sequenza e produce il fingerprint compresso.
 * Utilizza una finestra mobile per ottenere i k-mer e il loro valore codificato.
 */
ProcessedRead processReadFingerprint(const std::string &read, int k)
{
    ProcessedRead pr;
    int n = (int)read.size();
    if (n < k)
        return pr;  // Se la sequenza è troppo corta, ritorna vuoto

    int numKmers = n - k + 1;
    pr.fingerprint.reserve(numKmers);
    pr.kmers.reserve(numKmers);

    unsigned int mask = (1u << (2 * k)) - 1;  // Maschera per mantenere solo 2*k bit
    unsigned int val = 0;

    // Calcola il valore per il primo k-mer
    for (int i = 0; i < k; i++)
    {
        val <<= 2;
        char c = read[i];
        switch (std::toupper(static_cast<unsigned char>(c)))
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

    // Vettori per la compressione
    std::vector<unsigned int> comp_fp;
    std::vector<std::string_view> comp_kmers;
    std::vector<int> comp_indices;
    comp_fp.push_back(val);
    comp_kmers.push_back(std::string_view(read.data(), k));
    comp_indices.push_back(0);

    // Processa i successivi k-mer spostandosi di 1 carattere
    for (int i = 1; i < numKmers; i++)
    {
        // Aggiorna il valore usando lo shift e applica la maschera
        val = ((val << 2) & mask);
        char c = read[i + k - 1];
        switch (std::toupper(static_cast<unsigned char>(c)))
        {
            case 'A': break;
            case 'C': val |= 1; break;
            case 'G': val |= 2; break;
            case 'T': val |= 3; break;
            default: break;
        }
        pr.fingerprint.push_back(val);
        pr.kmers.push_back(std::string_view(read.data() + i, k));
        // Se il valore è diverso dall'ultimo compresso, lo aggiunge
        if (val != comp_fp.back())
        {
            comp_fp.push_back(val);
            comp_kmers.push_back(std::string_view(read.data() + i, k));
            comp_indices.push_back(i);
        }
    }
    // Salva i vettori compressi nella struttura ProcessedRead
    pr.comp.comp_fp = std::move(comp_fp);
    pr.comp.comp_kmers = std::move(comp_kmers);
    pr.comp.comp_indices = std::move(comp_indices);

    return pr;
}

/*
 * Funzione buildGlobalFingerprintFrequency:
 * Conta la frequenza di ogni k-mer (codificato) in tutte le sequenze.
 */
std::unordered_map<unsigned int,int> buildGlobalFingerprintFrequency(const std::vector<std::string> &reads, int k)
{
    std::unordered_map<unsigned int,int> freq_map;
    for (auto &read : reads)
    {
        int n = (int)read.size();
        if (n < k)
            continue;
        unsigned int mask = (1u << (2 * k)) - 1;
        unsigned int val = 0;
        // Calcola il primo k-mer
        for (int i = 0; i < k; i++)
        {
            val <<= 2;
            char c = read[i];
            switch (std::toupper(static_cast<unsigned char>(c)))
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
        // Processa gli altri k-mer
        for (int i = 1; i < numKmers; i++)
        {
            val = ((val << 2) & mask);
            char c = read[i + k - 1];
            switch (std::toupper(static_cast<unsigned char>(c)))
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

/*
 * Funzione processReadSolidFingerprint:
 * Come processReadFingerprint, ma scarta i k-mer che non sono nel set dei fingerprint solidi.
 */
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

    // Calcola il primo k-mer
    for (int i = 0; i < k; i++)
    {
        val <<= 2;
        char c = read[i];
        switch (std::toupper(static_cast<unsigned char>(c)))
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

    // Aggiunge il primo k-mer solo se è presente nel set dei fingerprint solidi
    if (solid_fingerprint_set.find(val) != solid_fingerprint_set.end())
    {
        pr.fingerprint.push_back(val);
        pr.kmers.push_back(std::string_view(read.data(), k));
        comp_fp.push_back(val);
        comp_kmers.push_back(std::string_view(read.data(), k));
        comp_indices.push_back(0);
    }

    // Processa gli altri k-mer
    for (int i = 1; i < numKmers; i++)
    {
        val = ((val << 2) & mask);
        char c = read[i + k - 1];
        switch (std::toupper(static_cast<unsigned char>(c)))
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
            if (comp_fp.empty())
            {
                comp_fp.push_back(val);
                comp_kmers.push_back(std::string_view(read.data() + i, k));
                comp_indices.push_back(i);
            }
            else
            {
                if (val != comp_fp.back())
                {
                    comp_fp.push_back(val);
                    comp_kmers.push_back(std::string_view(read.data() + i, k));
                    comp_indices.push_back(i);
                }
            }
        }
    }
    pr.comp.comp_fp = std::move(comp_fp);
    pr.comp.comp_kmers = std::move(comp_kmers);
    pr.comp.comp_indices = std::move(comp_indices);

    return pr;
}

/*
 * Funzione safe_substr:
 * Restituisce una sottostringa da 's' partendo da 'start' per 'length' caratteri in modo sicuro.
 */
std::string safe_substr(const std::string &s, size_t start, size_t length)
{
    if (start >= s.size())
        return "";
    size_t max_len = s.size() - start;
    if (length > max_len)
        length = max_len;
    return s.substr(start, length);
}

/*
 * Funzione get_overlap_annotation:
 * Se il match è inferiore a min_overlap, restituisce "(SCARTATA)".
 * Se ci sono troppe ripetizioni consecutive, restituisce un messaggio di scarto.
 * Altrimenti restituisce una stringa vuota.
 */
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
    // Scorre la regione per contare ripetizioni consecutive
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
        std::ostringstream oss;
        oss << "(SCARTATA - troppi valori consecutivi di '" << maxChar << "')";
        return oss.str();
    }
    return "";
}