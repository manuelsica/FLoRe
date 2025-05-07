#ifndef ICFL_HPP
#define ICFL_HPP
#include <string>
#include <vector>

/*  Canonical Inverse Lyndon Factorization  (ICFL)
 *  Ritorna i fattori   ι₁, … , ι_h   in O(n)  (Duval su ordine inverso + grouping). */
std::vector<std::string_view> icfl_factors(const std::string &s);

#endif
