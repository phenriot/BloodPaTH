// Standalone copy of the deterministic RNG helpers added to src/bloodpath_model_0.2.cpp,
// exposed to R via Rcpp::sourceCpp() so their output distributions can be compared against
// R's own rbinom/runif/rnorm/rlnorm/rpert without touching the BloodPaTH package itself
// (no NAMESPACE/RcppExports changes, no package (re)build required).
//
// Usage: see Validation/validate_rng.R in this folder.
//
// IMPORTANT: keep this file's function bodies in sync with src/bloodpath_model_0.2.cpp if the
// latter changes -- this is a deliberate, self-contained copy for validation purposes, not a
// shared header.

#include <Rcpp.h>
#include <cstdint>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

inline uint64_t v_splitmix64_next(uint64_t &state) {
  uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

inline double v_to_unit_double(uint64_t r) {
  return (r >> 11) * (1.0 / 9007199254740992.0);
}

inline int v_seeded_bernoulli(double seed_value, double prob) {
  uint64_t state = static_cast<uint64_t>(seed_value);
  double u = v_to_unit_double(v_splitmix64_next(state));
  return (u < prob) ? 1 : 0;
}

inline double v_seeded_unif(double seed_value, double min_v, double max_v) {
  uint64_t state = static_cast<uint64_t>(seed_value);
  double u = v_to_unit_double(v_splitmix64_next(state));
  return min_v + (max_v - min_v) * u;
}

inline double v_seeded_norm(double seed_value, double mean, double sd) {
  static const double TWO_PI = 6.283185307179586476925286766559;
  uint64_t state = static_cast<uint64_t>(seed_value);
  double u1 = std::max(v_to_unit_double(v_splitmix64_next(state)), 1e-300);
  double u2 = v_to_unit_double(v_splitmix64_next(state));
  double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(TWO_PI * u2);
  return mean + sd * z;
}

inline double v_seeded_lnorm(double seed_value, double meanlog, double sdlog) {
  return std::exp(v_seeded_norm(seed_value, meanlog, sdlog));
}

inline double v_seeded_gamma_from_state(uint64_t &state, double shape) {
  static const double TWO_PI = 6.283185307179586476925286766559;
  double d = shape - 1.0/3.0;
  double c = 1.0/std::sqrt(9.0*d);
  for(;;) {
    double x, v;
    do {
      double u1 = std::max(v_to_unit_double(v_splitmix64_next(state)), 1e-300);
      double u2 = v_to_unit_double(v_splitmix64_next(state));
      x = std::sqrt(-2.0*std::log(u1)) * std::cos(TWO_PI*u2);
      v = 1.0 + c*x;
    } while (v <= 0);
    v = v*v*v;
    double u = v_to_unit_double(v_splitmix64_next(state));
    double x2 = x*x;
    if (u < 1.0 - 0.0331*x2*x2) {return d*v;}
    if (std::log(u) < 0.5*x2 + d*(1.0 - v + std::log(v))) {return d*v;}
  }
}

inline double v_seeded_pert(double seed_value, double min_v, double mode_v, double max_v, double shape = 4.0) {
  if (max_v <= min_v) {return min_v;}
  double alpha1 = 1.0 + shape*(mode_v - min_v)/(max_v - min_v);
  double alpha2 = 1.0 + shape*(max_v - mode_v)/(max_v - min_v);
  uint64_t state = static_cast<uint64_t>(seed_value);
  double g1 = v_seeded_gamma_from_state(state, alpha1);
  double g2 = v_seeded_gamma_from_state(state, alpha2);
  double beta_draw = g1/(g1+g2);
  return min_v + (max_v - min_v)*beta_draw;
}

// [[Rcpp::export]]
IntegerVector val_bernoulli(NumericVector seeds, double prob) {
  IntegerVector out(seeds.size());
  for (int i = 0; i < seeds.size(); i++) {out[i] = v_seeded_bernoulli(seeds[i], prob);}
  return out;
}

// [[Rcpp::export]]
NumericVector val_unif(NumericVector seeds, double min_v, double max_v) {
  NumericVector out(seeds.size());
  for (int i = 0; i < seeds.size(); i++) {out[i] = v_seeded_unif(seeds[i], min_v, max_v);}
  return out;
}

// [[Rcpp::export]]
NumericVector val_norm(NumericVector seeds, double mean, double sd) {
  NumericVector out(seeds.size());
  for (int i = 0; i < seeds.size(); i++) {out[i] = v_seeded_norm(seeds[i], mean, sd);}
  return out;
}

// [[Rcpp::export]]
NumericVector val_lnorm(NumericVector seeds, double meanlog, double sdlog) {
  NumericVector out(seeds.size());
  for (int i = 0; i < seeds.size(); i++) {out[i] = v_seeded_lnorm(seeds[i], meanlog, sdlog);}
  return out;
}

// [[Rcpp::export]]
NumericVector val_pert(NumericVector seeds, double min_v, double mode_v, double max_v, double shape = 4.0) {
  NumericVector out(seeds.size());
  for (int i = 0; i < seeds.size(); i++) {out[i] = v_seeded_pert(seeds[i], min_v, mode_v, max_v, shape);}
  return out;
}
