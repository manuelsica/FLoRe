#ifndef CFL_HPP
#define CFL_HPP

#include <string>
#include <vector>

/* Ritorna la sequenza dei fattori di Lyndon di s (algoritmo di Duval). */
std::vector<std::string_view> cfl_factors(const std::string &s);

#endif
