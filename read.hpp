#ifndef READ_HPP
#define READ_HPP

#include <string>
#include <vector>

/*
 * Funzioni di lettura FASTA:
 *   read_fasta
 *   read_fasta_buffered (BFR)
 */
std::vector<std::string> read_fasta(const std::string &filename);
std::vector<std::string> read_fasta_buffered(const std::string &filename);

#endif