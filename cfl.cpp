#include "cfl.hpp"

// ——— helper “contains” e “find_index” (traduzione diretta dal Python) ———————————————————

static bool contains(const std::string &window,
                     const std::string &preview,
                     char next_c)
{
    // equivale a: ''.join(preview+[next_c]) in ''.join(window+preview)
    std::string dict = window + preview;
    std::string seq = preview + next_c;
    return dict.find(seq) != std::string::npos;
}

static size_t find_index(const std::string &window,
                         const std::string &preview)
{
    // equivale a: len(window) - ''.join(window+preview).find(preview)
    std::string dict = window + preview;
    auto pos = dict.find(preview);
    if (pos == std::string::npos)
        return std::string::npos;
    return window.size() - pos;
}

// ——— implementazione di Duval (CFL) ———————————————————————————————————————————

std::vector<std::string_view> cfl_factors(const std::string &s)
{
    size_t n = s.size(), k = 0;
    std::vector<std::string_view> out;

    while (k < n)
    {
        size_t i = k + 1;
        size_t j = k + 2;
        while (true)
        {
            if (j == n + 1 || s[j - 1] < s[i - 1])
            {
                // emetti i fattori di lunghezza (j-i)
                while (k < i)
                {
                    size_t len = j - i;
                    out.emplace_back(s.data() + k, len);
                    k += len;
                }
                break;
            }
            else
            {
                // restart di i o avanzamento
                if (s[j - 1] > s[i - 1])
                    i = k + 1;
                else
                    ++i;
                ++j;
            }
        }
    }

    return out;
}