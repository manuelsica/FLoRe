#ifndef READ_HPP
#define READ_HPP

#include <string>
#include <vector>
using namespace std;

/*
 * read_fasta:
 * Legge un file FASTA e restituisce un vettore di stringhe, ciascuna contenente la sequenza
 * (concatenando eventuali righe multiple per ogni record).
 */
vector<string> read_fasta(const string &filename);

#endif
