// cfl_icfl.cpp
#include "cfl_icfl.hpp"
#include "cfl.hpp"    // cfl_factors()
#include "icfl.hpp"   // icfl_factors()
#include <string>
#include <vector>
#include <string_view>

std::vector<std::string_view> cfl_icfl(const std::string &word, std::size_t C) {
    // 1) decomposizione CFL
    auto cfl_list = cfl_factors(word);

    // 2) per ogni fattore troppo lungo, esegue ICFL e avvolge con << … >>
    std::vector<std::string_view> out;
    out.reserve(cfl_list.size());
    for (auto &fv : cfl_list) {
        if (fv.size() <= C) {
            out.push_back(fv);
        }
        else {
            // marker inizio sub-decomposizione
            out.emplace_back("<<");

            // icfl_factors prende std::string, quindi creiamo una copia
            std::string tmp(fv);
            auto icfl_list = icfl_factors(tmp);

            for (auto &iv : icfl_list)
                out.push_back(iv);

            // marker fine
            out.emplace_back(">>");
        }
    }
    return out;
}
