// read.cpp
// Implementa le funzioni per leggere file FASTA in modalità tradizionale e buffered,
// e per pre-elaborare le sequenze (rimozione di spazi e conversione in uppercase).

#include "read.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>

/*
 * Helper interno: streaming linea per linea.
 */
static void read_fasta_stream(std::istream &in, std::vector<std::string> &reads) {
    std::string line, current_seq;
    reads.clear();
    while (std::getline(in, line)) {
        if (line.empty()) 
            continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                reads.push_back(std::move(current_seq));
                current_seq.clear();
            }
        } else {
            current_seq += line;
        }
    }
    if (!current_seq.empty())
        reads.push_back(std::move(current_seq));
}

std::vector<std::string> read_fasta(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Errore nell'apertura del file " << filename << std::endl;
        return {};
    }
    std::vector<std::string> reads;
    read_fasta_stream(infile, reads);
    return reads;
}

std::vector<std::string> read_fasta_buffered(const std::string &filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Errore nell'apertura del file " << filename << std::endl;
        return {};
    }
    std::vector<std::string> reads;
    read_fasta_stream(infile, reads);
    return reads;
}

std::vector<std::string> preprocess_reads(const std::vector<std::string> &raw_reads) {
    std::vector<std::string> processed;
    processed.reserve(raw_reads.size());
    for (const auto &s : raw_reads) {
        std::string seq;
        seq.reserve(s.size());
        for (unsigned char uc : s) {
            if (!std::isspace(uc))
                seq.push_back(static_cast<char>(std::toupper(uc)));
        }
        processed.push_back(std::move(seq));
    }
    return processed;
}
