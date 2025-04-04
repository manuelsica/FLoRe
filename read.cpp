// read.cpp
// Implementa le funzioni per la lettura dei file FASTA.

#include "read.hpp"
#include <fstream>      // per ifstream
#include <iostream>     // per std::cerr
#include <sstream>      // per istringstream

/*
 * Funzione read_fasta:
 * Legge un file FASTA in maniera tradizionale (riga per riga).
 */
std::vector<std::string> read_fasta(const std::string &filename)
{
    std::vector<std::string> reads;   // Vettore per contenere le sequenze
    std::ifstream infile(filename);   // Apre il file in modalità testo
    if (!infile)
    {
        std::cerr << "Errore nell'apertura del file " << filename << std::endl;
        return reads;
    }
    std::string line, current_seq;
    // Legge il file riga per riga
    while (std::getline(infile, line))
    {
        if (line.empty())
            continue;  // Salta righe vuote
        if (line[0] == '>')
        {
            // Se la riga inizia con '>', indica un header
            if (!current_seq.empty())
            {
                reads.push_back(current_seq);  // Salva la sequenza accumulata
                current_seq.clear();             // Resetta la sequenza
            }
        }
        else
        {
            // Accumula le righe che contengono la sequenza
            current_seq.append(line);
        }
    }
    // Aggiunge l'ultima sequenza se presente
    if (!current_seq.empty())
    {
        reads.push_back(current_seq);
    }
    return reads;
}

/*
 * Funzione read_fasta_buffered:
 * Legge l'intero file in un unico buffer per maggiore efficienza.
 */
std::vector<std::string> read_fasta_buffered(const std::string &filename)
{
    std::vector<std::string> reads;
    std::ifstream infile(filename, std::ios::binary);  // Apre il file in modalità binaria
    if (!infile)
    {
        std::cerr << "Errore nell'apertura del file " << filename << std::endl;
        return reads;
    }
    // Legge l'intero file in una stringa
    std::string fileContents((std::istreambuf_iterator<char>(infile)),
                             std::istreambuf_iterator<char>());
    size_t pos = 0;
    std::string current_seq;
    // Processa il buffer riga per riga
    while (pos < fileContents.size())
    {
        size_t line_end = fileContents.find('\n', pos);
        if (line_end == std::string::npos)
        {
            line_end = fileContents.size();
        }
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
        {
            current_seq.append(line);
        }
    }
    if (!current_seq.empty())
    {
        reads.push_back(current_seq);
    }
    return reads;
}