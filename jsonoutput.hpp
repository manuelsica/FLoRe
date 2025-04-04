// jsonoutput.hpp
// Header per la scrittura dei risultati in formato JSON.

#ifndef JSONOUTPUT_HPP
#define JSONOUTPUT_HPP

#include <vector>
#include <string>

/*
 * Struttura JsonResult:
 * Contiene le informazioni per ogni risultato JSON, inclusi gli indici delle read e la stringa JSON formattata.
 */
struct JsonResult
{
    int read1;
    int read2;
    std::string json;
};

/*
 * Funzione write_sorted_json:
 * Ordina i risultati JSON in base ai campi read1 e read2 e li scrive su un file.
 */
void write_sorted_json(std::vector<JsonResult> &json_results,
                       const std::string &filename);

#endif