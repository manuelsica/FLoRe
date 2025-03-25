#ifndef READ_HPP
#define READ_HPP

// Includiamo le librerie standard necessarie
#include <string>
#include <vector>
using namespace std;

// ---------------------------------------------------------------------------
// Funzione: read_fasta
// ---------------------------------------------------------------------------
// Descrizione:
//   Legge un file FASTA e restituisce un vettore di stringhe, in cui ogni
//   stringa corrisponde alla sequenza (concatenando tutte le righe di un record).
//
// Esempio di file FASTA:
//   >read1
//   ATCG
//   GTCA
//   >read2
//   TTAA
//
// L'output sarà: {"ATCGGTCA", "TTAA"}
// ---------------------------------------------------------------------------
vector<string> read_fasta(const string &filename);

#endif
