// jsonoutput.cpp
// Implementa la funzione per scrivere i risultati JSON ordinati in un file.

#include "jsonoutput.hpp"
#include <algorithm>    // per std::sort
#include <fstream>      // per std::ofstream
#include <iostream>     // per std::cerr

/*
 * Funzione write_sorted_json:
 * Ordina il vettore di JsonResult in base a read1 e read2 e lo scrive in un file JSON.
 */
void write_sorted_json(std::vector<JsonResult> &json_results,
                       const std::string &filename)
{
    // Ordina i risultati in ordine crescente di read1; se uguali, ordina per read2
    std::sort(json_results.begin(), json_results.end(), [](const JsonResult &a, const JsonResult &b){
        if (a.read1 != b.read1)
            return a.read1 < b.read1;
        return a.read2 < b.read2;
    });
    std::ofstream jf(filename);
    if (!jf.is_open())
    {
        std::cerr << "Errore apertura file " << filename << "\n";
        return;
    }
    // Scrive l'array JSON, con ogni oggetto separato da una virgola
    jf << "[\n";
    for (size_t i = 0; i < json_results.size(); i++)
    {
        jf << json_results[i].json;
        if (i < json_results.size() - 1)
            jf << ",\n";
    }
    jf << "\n]\n";
    jf.close();
}