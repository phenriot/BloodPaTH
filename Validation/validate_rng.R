# Validates that the native C++ RNG helpers added to src/bloodpath_model_0.2.cpp (as part of
# removing the model's per-draw calls into R's set.seed()+sample()/rlnorm()/rnorm()/rpert(),
# a major performance bottleneck) are statistically equivalent to their R counterparts.
#
# This does NOT build the BloodPaTH package -- it compiles a small standalone copy of the RNG
# helpers (Validation/validate_rng.cpp) via Rcpp::sourceCpp() and compares empirical
# distributions against base R / mc2d, using many distinct (fake) seeds like the model would
# produce (one per (time, patient, id_sim, draw-type) combination).
#
# Requires: Rcpp, mc2d. Run from the package root, or adjust the path below.

library(Rcpp)
library(mc2d)

sourceCpp("Validation/validate_rng.cpp")

n <- 200000
seeds <- round(runif(n, 1, 1e7)) # stand-in for the model's round((((time+1)+(p+1))/(p+1))*id_sim*K)

cat("========================================\n")
cat("Bernoulli(p) via val_bernoulli() vs rbinom(1, p)\n")
cat("========================================\n")
for (p in c(0.05, 0.3, 0.7, 0.95)) {
  draws <- val_bernoulli(seeds, p)
  cat(sprintf("p=%.2f | empirical mean=%.4f (expected %.4f)\n", p, mean(draws), p))
}

cat("\n========================================\n")
cat("Uniform(min, max) via val_unif() vs runif()\n")
cat("========================================\n")
draws <- val_unif(seeds, 4, 20)
ref <- runif(n, 4, 20)
print(ks.test(draws, ref))
cat(sprintf("mean: new=%.4f  runif()=%.4f (expected 12)\n", mean(draws), mean(ref)))

cat("\n========================================\n")
cat("Normal(mean, sd) via val_norm() vs rnorm()\n")
cat("========================================\n")
draws <- val_norm(seeds, 0.02, 0.01)
ref <- rnorm(n, 0.02, 0.01)
print(ks.test(draws, ref))
cat(sprintf("mean: new=%.5f  rnorm()=%.5f | sd: new=%.5f  rnorm()=%.5f\n",
            mean(draws), mean(ref), sd(draws), sd(ref)))

cat("\n========================================\n")
cat("Log-Normal(meanlog, sdlog) via val_lnorm() vs rlnorm()\n")
cat("========================================\n")
draws <- val_lnorm(seeds, -4, 0.5)
ref <- rlnorm(n, -4, 0.5)
print(ks.test(draws, ref))
cat(sprintf("mean: new=%.5f  rlnorm()=%.5f\n", mean(draws), mean(ref)))

cat("\n========================================\n")
cat("PERT(min, mode, max) via val_pert() vs mc2d::rpert() [most important check]\n")
cat("========================================\n")
draws <- val_pert(seeds, 0, 20, 100)
ref <- rpert(n, min = 0, mode = 20, max = 100)
print(ks.test(draws, ref))
cat(sprintf("mean: new=%.3f  rpert()=%.3f | sd: new=%.3f  rpert()=%.3f\n",
            mean(draws), mean(ref), sd(draws), sd(ref)))

cat("\n========================================\n")
cat("Determinism check: same seed -> same draw, for all 5 functions\n")
cat("========================================\n")
s <- c(123456, 42, 987654321)
stopifnot(identical(val_bernoulli(s, 0.3), val_bernoulli(s, 0.3)))
stopifnot(identical(val_unif(s, 4, 20), val_unif(s, 4, 20)))
stopifnot(identical(val_norm(s, 0, 1), val_norm(s, 0, 1)))
stopifnot(identical(val_lnorm(s, 0, 1), val_lnorm(s, 0, 1)))
stopifnot(identical(val_pert(s, 0, 20, 100), val_pert(s, 0, 20, 100)))
cat("OK: all five functions are deterministic per seed.\n")

cat("\nDone. For each block above, the KS test p-value should not be small/significant,\n")
cat("and the reported means/sds should be close -- if not, do not trust the corresponding\n")
cat("branch of src/bloodpath_model_0.2.cpp's RNG replacement.\n")
