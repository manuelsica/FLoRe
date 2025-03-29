#ifndef READ_HPP
#define READ_HPP

#include <string>
#include <vector>

/*
 * read_fasta:
 *   Lettura standard di un file FASTA.
 * read_fasta_buffered:
 *   Versione “Buffered FASTA Reader (BFR)” per file grandi.
 */
std::vector<std::string> read_fasta(const std::string &filename);
std::vector<std::string> read_fasta_buffered(const std::string &filename);

#endif
