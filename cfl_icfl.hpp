// cfl_icfl.hpp
#ifndef CFL_ICFL_HPP
#define CFL_ICFL_HPP

#include <string>
#include <vector>
#include <string_view>

/*
 * CFL : fattorizzazione di Lyndon (Duval)
 * ICFL : factorization inverse di Lyndon + grouping
 *
 * Per ogni fattore di Lyndon > C invoca ICFL, incapsulando la sub-decomposizione
 * tra i marker "<<" e ">>".
 */
std::vector<std::string_view> cfl_icfl(const std::string &word, std::size_t C);

#endif // CFL_ICFL_HPP
