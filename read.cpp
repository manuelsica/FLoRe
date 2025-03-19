#include "read.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

vector<string> read_fasta(const string &filename) {
    vector<string> reads;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Errore nell'apertura del file " << filename << endl;
        return reads;
    }
    string line;
    string current_seq;
    while(getline(infile, line)) {
        if(line.empty())
            continue;
        if(line[0] == '>') {
            if(!current_seq.empty()) {
                reads.push_back(current_seq);
                current_seq.clear();
            }
            // Ignora la riga di header
        } else {
            current_seq += line;
        }
    }
    if(!current_seq.empty())
        reads.push_back(current_seq);
    return reads;
}
