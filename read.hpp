// read.hpp
// Header per le funzioni di lettura dei file FASTA.
// Supporta esclusivamente il formato FASTA.
#ifndef READ_HPP
#define READ_HPP

#include <string>
#include <vector>

/*
 * read_fasta: Legge un file FASTA riga per riga e restituisce un vettore di sequenze.
 */
std::vector<std::string> read_fasta(const std::string &filename);

/*
 * read_fasta_buffered: Legge un file FASTA in modalità binary, ma in streaming
 * per prestazioni migliori senza duplicare il buffer completo.
 */
std::vector<std::string> read_fasta_buffered(const std::string &filename);

/*
 * preprocess_reads: Esegue la pre-elaborazione delle sequenze, rimuovendo spazi
 * e convertendo in uppercase.
 */
std::vector<std::string> preprocess_reads(const std::vector<std::string> &raw_reads);

#endif // READ_HPP
