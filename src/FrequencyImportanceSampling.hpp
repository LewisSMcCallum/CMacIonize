#ifndef FREQUENCYIMPORTANCESAMPLING_HPP
#define FREQUENCYIMPORTANCESAMPLING_HPP

#include "Error.hpp"
#include "RandomGenerator.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

/** Sample q=(1-a)p+a*uniform and return the photon-NUMBER weight p/q.
 * p is the piecewise-constant density represented by the existing linear CDF.
 * No spectral range, source luminosity or cross section is changed. */
inline double sample_frequency_importance(
    const std::vector<double> &frequency, const std::vector<double> &cdf,
    RandomGenerator &random, const double uniform_fraction, double &weight) {
  if (!std::isfinite(uniform_fraction) || uniform_fraction <= 0. ||
      uniform_fraction > 1.) {
    cmac_error("Importance sampler requires a uniform fraction in (0,1].");
  }
  if (frequency.size() < 2 || cdf.size() != frequency.size()) {
    cmac_error("Invalid spectrum table for frequency importance sampling.");
  }
  const double width = frequency.back() - frequency.front();
  double nu;
  if (random.get_uniform_random_double() < uniform_fraction) {
    nu = frequency.front() + width * random.get_uniform_random_double();
  } else {
    const double u = random.get_uniform_random_double();
    const size_t i = std::min(size_t(std::upper_bound(cdf.begin(), cdf.end(), u)
                                    - cdf.begin() - 1), cdf.size() - 2);
    nu = frequency[i] + (frequency[i+1] - frequency[i]) *
                           (u - cdf[i]) / (cdf[i+1] - cdf[i]);
  }
  const size_t i = std::min(size_t(std::upper_bound(
      frequency.begin(), frequency.end(), nu) - frequency.begin() - 1),
      frequency.size() - 2);
  const double p = (cdf[i+1] - cdf[i]) / (frequency[i+1] - frequency[i]);
  weight = p / ((1. - uniform_fraction) * p + uniform_fraction / width);
  if (!std::isfinite(weight) || weight < 0.) {
    cmac_error("Invalid spectral packet weight: p=%g, frequency=%g.", p, nu);
  }
  return nu;
}
#endif
