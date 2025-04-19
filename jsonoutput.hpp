// jsonoutput.hpp
// Header per il modulo di output JSON.
#ifndef JSONOUTPUT_HPP
#define JSONOUTPUT_HPP

#include <vector>
#include <string>

// Struttura per memorizzare il risultato JSON associato a una coppia di read.
struct JsonResult {
    int read1;            // Indice della prima read (1-indexed)
    int read2;            // Indice della seconda read (1-indexed)
    std::string json;     // Stringa JSON formattata
};

/*
 * write_sorted_json: Ordina i risultati per read1 e read2 e li scrive in un file JSON.
 */
void write_sorted_json(std::vector<JsonResult> &json_results,
                       const std::string &filename);

#endif