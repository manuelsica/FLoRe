// jsonoutput.cpp
#include "jsonoutput.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>

void write_sorted_json(std::vector<JsonResult> &json_results,
                       const std::string &filename)
{
    // Ordina per read1, poi per read2
    std::sort(json_results.begin(),
              json_results.end(),
              [](const JsonResult &a, const JsonResult &b){
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

    jf << "[\n";
    for (size_t i = 0; i < json_results.size(); ++i)
    {
        jf << json_results[i].json;
        if (i + 1 < json_results.size())
            jf << ",\n";
    }
    jf << "\n]\n";
}
