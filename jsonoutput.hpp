// jsonoutput.hpp
#ifndef JSONOUTPUT_HPP
#define JSONOUTPUT_HPP

#include <vector>
#include <string>

// Struttura per memorizzare il risultato JSON associato a una coppia di read.
struct JsonResult {
    int read1;            // Indice della prima read (1-indexed)
    int read2;            // Indice della seconda read (1-indexed)
    std::string json;     // Stringa JSON formattata

    // Costruttore per abilitare emplace_back(read1, read2, json)
    JsonResult(int r1, int r2, const std::string &j)
      : read1(r1), read2(r2), json(j) {}
};

/// Ordina e scrive i risultati in un file JSON.
void write_sorted_json(std::vector<JsonResult> &json_results,
                       const std::string &filename);

#endif // JSONOUTPUT_HPP
