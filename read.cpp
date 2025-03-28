#include "read.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

vector<string> read_fasta(const string &filename) {
    vector<string> reads;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Errore nell'apertura del file " << filename << endl;
        return reads;
    }
    
    string line;
    size_t readCount = 0;
    // Primo passaggio: contare il numero di reads (righe che iniziano con '>')
    while (getline(infile, line)) {
        if (!line.empty() && line[0] == '>')
            readCount++;
    }
    
    // Riposiziona il file all'inizio
    infile.clear();
    infile.seekg(0);
    
    // Pre-alloca il vettore in base al conteggio ottenuto
    reads.reserve(readCount);
    string current_seq;
    
    // Secondo passaggio: lettura delle sequenze
    while (getline(infile, line)) {
        if (line.empty())
            continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                reads.push_back(current_seq);
                current_seq.clear();
            }
        } else {
            current_seq.append(line);
        }
    }
    if (!current_seq.empty())
        reads.push_back(current_seq);
    
    return reads;
}
