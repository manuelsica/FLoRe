// read.hpp
// Header per le funzioni di lettura dei file FASTA.

#ifndef READ_HPP
#define READ_HPP

#include <string>
#include <vector>

/*
 * Funzione read_fasta:
 * Legge un file FASTA riga per riga e restituisce un vettore di sequenze.
 */
std::vector<std::string> read_fasta(const std::string &filename);

/*
 * Funzione read_fasta_buffered:
 * Legge l'intero file FASTA in un buffer per prestazioni migliori,
 * e restituisce un vettore di sequenze.
 */
std::vector<std::string> read_fasta_buffered(const std::string &filename);

#endif