// filter.hpp
#ifndef FILTER_HPP
#define FILTER_HPP
#include "util.hpp"  // per ProcessedRead e processReadFingerprint
#include <string>

namespace pseudo_overlap {

/*
 * entropy:
 *   Calcola l’entropia di Shannon di una regione di DNA.
 */
double entropy(const std::string &region);

/*
 * low_complexity:
 *   Ritorna true se l’entropia < min_entropy.
 */
bool low_complexity(const std::string &region, double min_entropy = 1.5);

/*
 * fingerprint_chained_local_align:
 *   FCLA: fingerprint-chaining + local extension sui gap.
 */
bool fingerprint_chained_local_align(const std::string &r1,
                                     const std::string &r2,
                                     int k,
                                     double min_identity,
                                     int max_gap);

/*
 * spectrum_similarity:
 *   Calcola distribuzioni di k-mer (default k=3) usando codifica bit,
 *   misura la Jensen–Shannon divergence e ritorna true se ≤ max_js.
 */
bool spectrum_similarity(const std::string &r1,
                         const std::string &r2,
                         int kmer_len = 3,
                         double max_js = 0.5);

/*
 * block_entropy_consistency:
 *   Divide la regione in finestre (block_size, step block_step),
 *   calcola entropia su ciascuna finestra via string_view
 *   e ritorna true se varianza ≤ max_var.
 */
bool block_entropy_consistency(const std::string &region,
                               int block_size = 10,
                               int block_step = 5,
                               double max_var = 0.25);

} // namespace pseudo_overlap

#endif // FILTER_HPP
