context("Kernel functions")

test_that("canonical kernel works", {

  # x is a vector
  x <- 1:2
  res <- kern_canonical(x)
  # > res
  #      [,1] [,2]
  # [1,]    1    2
  # [2,]    2    4

  # x is a matrix
  x <- matrix(1:10, ncol = 5)
  res <- kern_canonical(x)
  # > res
  #      [,1] [,2]
  # [1,]  165  190
  # [2,]  190  220
})

