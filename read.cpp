// Includiamo il header del modulo e le librerie necessarie
#include "read.hpp"       // Dichiarazione di read_fasta
#include <fstream>        // Per ifstream
#include <sstream>        // Per istringstream (se necessario)
#include <iostream>       // Per cout, cerr
#include <vector>         // Per vector
#include <string>         // Per string
using namespace std;

// ---------------------------------------------------------------------------
// Funzione: read_fasta
// ---------------------------------------------------------------------------
// Descrizione:
//   Legge il file FASTA specificato dal parametro filename, concatenando
//   le righe non vuote di ciascun record (ignora le righe che iniziano con '>')
//   e restituisce un vettore di stringhe, ognuna contenente una sequenza.
// ---------------------------------------------------------------------------
vector<string> read_fasta(const string &filename) {
    vector<string> reads;           // Vettore in cui verranno salvate le sequenze lette.

    // Apro il file in modalità binaria e posiziono il cursore alla fine (ios::ate)
    ifstream infile(filename, ios::binary | ios::ate);
    if (!infile) {
        cerr << "Errore nell'apertura del file " << filename << endl;
        return reads;  // Ritorno un vettore vuoto in caso di errore.
    }

    // Determino la dimensione del file
    streamsize fileSize = infile.tellg();
    // Torno all'inizio del file per poter leggere tutto
    infile.seekg(0, ios::beg);

    // Alloco una stringa grande quanto il file
    string buffer(fileSize, '\0');
    // Leggo l'intero file in un unico passaggio
    if (!infile.read(&buffer[0], fileSize)) {
        cerr << "Errore durante la lettura del file " << filename << endl;
        return reads;
    }

    // Variabili per il parsing
    istringstream iss(buffer);
    string line;
    string current_seq;

    // Leggo riga per riga dal buffer (che contiene tutto il file)
    while (std::getline(iss, line)) {
        if (line.empty()) // Se la riga è vuota, la ignoro.
            continue;
        if (line[0] == '>') {  // Se la riga inizia con '>', è un header FASTA.
            // Se current_seq non è vuota, significa che abbiamo finito di leggere una sequenza.
            if (!current_seq.empty()) {
                reads.push_back(current_seq);  // Aggiungo la sequenza al vettore delle reads.
                current_seq.clear();           // Resetto current_seq per la prossima sequenza.
            }
            // L'header può essere processato se necessario, ma qui viene ignorato.
        } else {
            // Se la riga non è un header, la concateno a current_seq.
            // Esempio: "ATCG" e poi "GTCA" diventeranno "ATCGGTCA".
            current_seq += line;
        }
    }

    // Dopo aver letto tutte le righe, se c'è ancora una sequenza in current_seq, la aggiungo.
    if (!current_seq.empty())
        reads.push_back(current_seq);

    // Ritorno il vettore contenente tutte le sequenze lette.
    return reads;
}
