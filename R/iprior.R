n <- 100
x <- seq(0, 10, length = n)
y <- 1 + 2 * x + rt(n, df = 2)
qplot(x, y)

iprior_noscale_nopsi <- function(y, ..., kernel = "canonical", psi = 1) {
  x_list <- list(...)
  x_list <- lapply(x_list, scale)
  y <- scale(y)

  # Build kernel matrices
  kernel_matrix_list <- lapply(x_list, kernel_matrix, kernel = kernel)
  K <- Reduce("+", lapply(kernel_matrix_list, as.matrix))

  # Calculated fitted values
  K_eigen <- iprior::eigenCpp(K)
  u <- K_eigen$val
  V <- K_eigen$vec
  w <- as.numeric(crossprod(y, V) %*% (t(V) / (psi * u ^ 2 + 1 / psi)))

  y_hat <-  V %*% (t(V) * (psi * u ^ 2)) %*% w
  y_hat <- scale_back(y_hat, get_scale(y), get_centre(y))

  res <- list(x = x_list, y = y, fitted.values = as.numeric(y_hat))
  class(res) <- "iprior_fitted"
  res
}

iprior_noscale <- function(y, ..., kernel = "canonical") {
  n <- length(y)
  x_list <- list(...)
  x_list <- lapply(x_list, scale)
  y <- scale(y)

  # Build kernel matrices
  kernel_matrix_list <- lapply(x_list, kernel_matrix, kernel = kernel)
  K <- Reduce("+", lapply(kernel_matrix_list, as.matrix))

  # Calculated fitted values
  K_eigen <- iprior::eigenCpp(K)
  u <- K_eigen$val
  V <- K_eigen$vec

  # Optimise psi
  optim_res <- optim(0, loglik_iprior_noscale, method = "L-BFGS",
                     u = u, V = V, y = y, control = list(fnscale = -2))
  psi <- exp(optim_res$par)
  w <- as.numeric(crossprod(y, V) %*% (t(V) / (psi * u ^ 2 + 1 / psi)))

  y_hat <- psi * V %*% (t(V) * (u ^ 2)) %*% w
  y_hat <- scale_back(y_hat, get_scale(y), get_centre(y))

  res <- list(x = x_list, y = y, fitted.values = as.numeric(y_hat),
              optim = optim_res)
  class(res) <- "iprior_fitted"
  res
}

loglik_iprior_noscale <- function(theta, u, V, y) {
  psi <- exp(theta)
  n <- length(y)
  w <- as.numeric(crossprod(y, V) %*% (t(V) / (psi * u ^ 2 + 1 / psi)))
  res <- -(n / 2) * log(2 * pi) - sum(log(u ^ 2 + 1 / psi)) / 2 -
    0.5 * crossprod(y, w)
  # if (isTRUE(deviance)) -2 * as.numeric(res)
  # else res
  res
}

iprior_scale <- function(y, ..., kernel = "canonical") {
  n <- length(y)
  x_list <- list(...)
  x_list <- lapply(x_list, function(x) scale(x, center = TRUE, scale = FALSE))
  y <- scale(y, center = TRUE, scale = FALSE)

  # Build kernel matrices
  kernel_matrix_list <- lapply(x_list, kernel_matrix, kernel = kernel)

  # Optimise psi
  optim_res <- optim(rnorm(length(kernel_matrix_list) + 1),
                     loglik_iprior_scale, method = "CG",
                     kernel_matrix_list = kernel_matrix_list, y = y,
                     control = list(fnscale = -2))
  theta <- optim_res$par
  lambda <- exp(theta)[-length(theta)]
  psi <- exp(theta)[length(theta)]
  K <- Reduce("+", mapply("*", lapply(kernel_matrix_list, as.matrix), lambda,
                          SIMPLIFY = FALSE))
  K_eigen <- iprior::eigenCpp(K)
  u <- K_eigen$val
  V <- K_eigen$vec
  z <- psi * u ^ 2 + 1 / psi
  w <- (V * rep(1 / z, each = nrow(V))) %*% crossprod(V, y)
  y_hat <- psi * V %*% (t(V) * (u ^ 2)) %*% w
  y_hat <- y_hat + get_centre(y)

  res <- list(x = x_list, y = y, fitted.values = as.numeric(y_hat),
              optim = optim_res)
  class(res) <- "iprior_fitted"
  res
}

loglik_iprior_scale <- function(theta, kernel_matrix_list, y) {
  lambda <- exp(theta)[-length(theta)]
  psi <- exp(theta)[length(theta)]
  n <- length(y)

  K <- Reduce("+", mapply("*", lapply(kernel_matrix_list, as.matrix), lambda,
                          SIMPLIFY = FALSE))
  K_eigen <- iprior::eigenCpp(K)
  u <- K_eigen$val
  V <- K_eigen$vec
  w <- as.numeric(crossprod(y, V) %*% (t(V) / (psi * u ^ 2 + 1 / psi)))
  res <- -(n / 2) * log(2 * pi) - sum(log(u ^ 2 + 1 / psi)) / 2 -
    0.5 * crossprod(y, w)
  res
}

loglik_test <- function(theta, x, y) {
  x <- scale(x, scale = FALSE)
  y <- scale(y, scale = FALSE)
  lambda <- exp(theta)[-length(theta)]
  psi <- exp(theta)[length(theta)]
  n <- length(y)

  K <- as.matrix(kernel_matrix(x))
  K <- lambda * K
  K_eigen <- iprior::eigenCpp(K)
  u <- K_eigen$val
  V <- K_eigen$vec
  z <- psi * u ^ 2 + 1 / psi
  w <- (V * rep(1 / z, each = nrow(V))) %*% crossprod(V, y)
  res <- -(n / 2) * log(2 * pi) - sum(log(u ^ 2 + 1 / psi)) / 2 -
    0.5 * crossprod(y, w)
  as.numeric(res)

}

get_scale <- function(x) attr(x, "scaled:scale")

get_centre <- function(x) attr(x, "scaled:center")

scale_back <- function(x, scale = NULL, centre = NULL) {
  if (is.null(scale)) scale <- get_scale(x)
  if (is.null(centre)) centre <- get_centre(x)
  if (is.null(scale)) scale <- 1
  if (is.null(centre)) centre <- 0
  res <- sweep(x, 2, scale, "*")
  res <- sweep(res, 2, centre, "+")
  res
}

plot.iprior_fitted <- function(x) {
  plot_df <- data.frame(
    x = scale_back(x$x[[1]]),
    y = scale_back(x$y),
    fitted = x$fitted.values
  )

  ggplot(plot_df, aes(x, y)) +
    geom_point() +
    geom_line(aes(y = fitted)) +
    theme_bw()
}


