// icfl.cpp
#include "icfl.hpp"
#include <algorithm>    // std::equal
#include <vector>
#include <string_view>

/*
 * Duval’s algorithm su **lex-order inverso**:
 * genera la fattorizzazione di Lyndon rispetto
 * all’ordine inverso (>=).
 */
static std::vector<std::string_view> cfl_inv(const std::string &s)
{
    size_t n = s.size(), i = 0;
    std::vector<std::string_view> out;
    while (i < n) {
        size_t j = i + 1, k = i;
        while (j < n && s[k] >= s[j]) {    // “>=” per ordine inverso
            if (s[k] > s[j])
                k = i;                     // restart
            else
                ++k;
            ++j;
        }
        size_t l = j - k;
        while (i <= k) {
            out.emplace_back(s.data() + i, l);
            i += l;
        }
    }
    return out;
}

/*
 * Controlla se `small` è prefisso di `big`.
 * (utility per il grouping successivo)
 */
static inline bool is_prefix(std::string_view small, std::string_view big)
{
    return small.size() <= big.size()
        && std::equal(small.begin(), small.end(), big.begin());
}

/*
 * Raggruppamento: seviene una sequenza di fattori di Lyndon inversi
 * ciascun fattore v[j] che è prefisso di v[i] viene “accorpato”
 * in un unico fattore esteso.
 */
std::vector<std::string_view> icfl_factors(const std::string &s)
{
    auto v = cfl_inv(s);
    std::vector<std::string_view> res;
    for (size_t i = 0; i < v.size(); ) {
        size_t j = i + 1;
        // allarga finché v[j] è prefisso di v[i]
        while (j < v.size() && is_prefix(v[j], v[i]))
            ++j;
        // calcola la lunghezza totale tra inizio di v[i] e fine di v[j-1]
        size_t len = (v[j-1].data() + v[j-1].size()) - v[i].data();
        res.emplace_back(v[i].data(), len);
        i = j;
    }
    return res;
}
