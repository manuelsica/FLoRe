#include <iostream>
#include <string>
#include <vector>
#include "util.hpp"
#include "overlap.hpp"
#include "read.hpp"
#include <tuple>
#include <sstream>
using namespace std;

struct OverlapResult {
    int overlap_len;
    int start1;
    int end1;
    int start2;
    int end2;
    string orientation1;
    string orientation2;
    string combination;
    string r1;
    string r2;
};

struct Combo {
    string label;
    ProcessedRead pr1;
    ProcessedRead pr2;
    string orient1;
    string orient2;
    string orig1;  // read originale (forward o reverse complement) per il confronto
    string orig2;
};

struct ReadData {
    string forward;
    string reverse;
    ProcessedRead pr_fwd;
    ProcessedRead pr_rev;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Utilizzo: " << argv[0] << " <file_fasta>" << endl;
        return 1;
    }
    string fasta_file = argv[1];
    vector<string> reads = read_fasta(fasta_file);
    if (reads.empty()) {
        cerr << "Nessuna read trovata nel file." << endl;
        return 1;
    }
    
    int k = 15;
    
    // Pre-processa tutte le read in parallelo:
    vector<ReadData> all_reads(reads.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < reads.size(); i++) {
        all_reads[i].forward = reads[i];
        all_reads[i].pr_fwd = process_read(reads[i], k);
        string rev = reverse_complement(reads[i]);
        all_reads[i].reverse = rev;
        all_reads[i].pr_rev = process_read(rev, k);
    }
    
    // Per ogni coppia di read, confronta tutte le 4 combinazioni di orientamento
    for (size_t i = 0; i < all_reads.size(); i++) {
        for (size_t j = i + 1; j < all_reads.size(); j++) {
            OverlapResult bestResult;
            bestResult.overlap_len = 0;
            vector<Combo> comparisons = {
                { "ff", all_reads[i].pr_fwd, all_reads[j].pr_fwd, "forward", "forward", all_reads[i].forward, all_reads[j].forward },
                { "fr", all_reads[i].pr_fwd, all_reads[j].pr_rev, "forward", "reverse", all_reads[i].forward, all_reads[j].reverse },
                { "rf", all_reads[i].pr_rev, all_reads[j].pr_fwd, "reverse", "forward", all_reads[i].reverse, all_reads[j].forward },
                { "rr", all_reads[i].pr_rev, all_reads[j].pr_rev, "reverse", "reverse", all_reads[i].reverse, all_reads[j].reverse }
            };
            
            for (auto &comp : comparisons) {
                int overlap_len, start1, end1, start2, end2;
                tie(overlap_len, start1, end1, start2, end2) =
                    graph_overlap_fp(comp.pr1.comp.comp_fp, comp.pr2.comp.comp_fp,
                                     comp.pr1.comp.comp_indices, comp.pr2.comp.comp_indices, k);
                if (overlap_len > bestResult.overlap_len) {
                    bestResult.overlap_len = overlap_len;
                    bestResult.combination = comp.label;
                    bestResult.orientation1 = comp.orient1;
                    bestResult.orientation2 = comp.orient2;
                    bestResult.start1 = start1;
                    bestResult.end1 = end1;
                    bestResult.start2 = start2;
                    bestResult.end2 = end2;
                    bestResult.r1 = comp.orig1;
                    bestResult.r2 = comp.orig2;
                }
            }
            
            // Stampa il risultato per la coppia se viene trovato un overlap
            if (bestResult.overlap_len > 0) {
                cout << "Coppia di read " << i+1 << " e " << j+1 << ":" << endl;
                cout << "La piu' lunga sottosequenza contigua comune (graph-based su fingerprint) ha lunghezza (in match fingerprint): " 
                     << bestResult.overlap_len << endl;
                cout << "Regione di overlap in Read " << i+1 << " (" << bestResult.orientation1 << "): " 
                     << bestResult.r1.substr(bestResult.start1, bestResult.end1 - bestResult.start1) << endl;
                cout << "Regione di overlap in Read " << j+1 << " (" << bestResult.orientation2 << "): " 
                     << bestResult.r2.substr(bestResult.start2, bestResult.end2 - bestResult.start2) << endl;
                cout << "Tabella riassuntiva:" << endl;
                cout << "Read1\tRead2\tOrientation\tStart1\tEnd1\tLen(Read1)\tStart2\tEnd2\tLen(Read2)" << endl;
                int orient_flag = (bestResult.combination == "ff" || bestResult.combination == "rr") ? 0 : 1;
                cout << i+1 << "\t" << j+1 << "\t" << orient_flag << "\t"
                     << bestResult.start1 << "\t" << bestResult.end1 << "\t" << bestResult.r1.size() << "\t"
                     << bestResult.start2 << "\t" << bestResult.end2 << "\t" << bestResult.r2.size() << endl;
                cout << "---------------------------------------" << endl;
            } else {
                cout << "Coppia di read " << i+1 << " e " << j+1 << ": Nessuna sottosequenza comune trovata." << endl;
            }
        }
    }
    
    return 0;
}
