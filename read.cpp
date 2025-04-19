// read.cpp
// Implementa le funzioni per leggere file FASTA in modalità tradizionale e buffered,
// e una funzione per pre-elaborare le sequenze (rimozione di spazi e conversione in uppercase).
#include "read.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>

std::vector<std::string> read_fasta(const std::string &filename)
{
    std::vector<std::string> reads;
    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Errore nell'apertura del file " << filename << std::endl;
        return reads;
    }
    std::string line, current_seq;
    while (std::getline(infile, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '>')
        {
            if (!current_seq.empty())
            {
                reads.push_back(current_seq);
                current_seq.clear();
            }
        }
        else
            current_seq.append(line);
    }
    if (!current_seq.empty())
        reads.push_back(current_seq);
    return reads;
}

std::vector<std::string> read_fasta_buffered(const std::string &filename)
{
    std::vector<std::string> reads;
    std::ifstream infile(filename, std::ios::binary);
    if (!infile)
    {
        std::cerr << "Errore nell'apertura del file " << filename << std::endl;
        return reads;
    }
    std::string fileContents((std::istreambuf_iterator<char>(infile)),
                             std::istreambuf_iterator<char>());
    size_t pos = 0;
    std::string current_seq;
    while (pos < fileContents.size())
    {
        size_t line_end = fileContents.find('\n', pos);
        if (line_end == std::string::npos)
            line_end = fileContents.size();
        std::string line = fileContents.substr(pos, line_end - pos);
        pos = line_end + 1;
        if (line.empty())
            continue;
        if (line[0] == '>')
        {
            if (!current_seq.empty())
            {
                reads.push_back(current_seq);
                current_seq.clear();
            }
        }
        else
            current_seq.append(line);
    }
    if (!current_seq.empty())
        reads.push_back(current_seq);
    return reads;
}

std::vector<std::string> preprocess_reads(const std::vector<std::string> &raw_reads)
{
    std::vector<std::string> processed;
    processed.reserve(raw_reads.size());
    for (const auto &s : raw_reads)
    {
        std::string seq;
        seq.reserve(s.size());
        for (char c : s)
        {
            if (!isspace(static_cast<unsigned char>(c)))
                seq.push_back(toupper(static_cast<unsigned char>(c)));
        }
        processed.push_back(seq);
    }
    return processed;
}