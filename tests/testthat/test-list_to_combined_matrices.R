test_that("type = 1 row-binds transition matrices and preserves values", {
  m1 <- matrix(1:6, nrow = 2, ncol = 3)
  m2 <- matrix(7:12, nrow = 2, ncol = 3)

  combined <- list_to_combined_matrices(list(m1, m2), type = 1, nb_wards = 2)

  expect_equal(dim(combined), c(4, 3))
  expect_equal(combined, rbind(m1, m2))
})

test_that("type = 2 row-binds procedure probability matrices and preserves values", {
  m1 <- matrix(1:4, nrow = 2, ncol = 2)
  m2 <- matrix(5:8, nrow = 2, ncol = 2)

  combined <- list_to_combined_matrices(list(m1, m2), type = 2, nb_procedures = 2)

  expect_equal(dim(combined), c(4, 2))
  expect_equal(combined, rbind(m1, m2))
})

test_that("type = 1 errors without nb_wards", {
  m1 <- matrix(1:4, nrow = 2, ncol = 2)
  expect_error(list_to_combined_matrices(list(m1), type = 1))
})

test_that("type = 2 errors without nb_procedures", {
  m1 <- matrix(1:4, nrow = 2, ncol = 2)
  expect_error(list_to_combined_matrices(list(m1), type = 2))
})
