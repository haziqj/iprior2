kernel_matrix <- function(x, y = NULL, kernel = "canonical") {
  # Check argument
  kernel <- match.arg(kernel, c("canonical"))

  # Check dimensions
  if (is.matrix(x) & is.vector(y)) {
    y <- matrix(y, nrow = 1)
    if (ncol(x) != ncol(y)) {
      stop("Dimensions of x and y do not match.")
    }
  }
  if (is.vector(x) & is.matrix(y)) {
    stop("y should be a vector.")
  }

  # Pass to correct kernel function
  if (kernel == "canonical") kern_function <- kern_canonical

  res <- list(matrix = kern_function(x, y), kernel = kernel)
  class(res) <- "kernel_matrix"
  res
}

print.kernel_matrix <- function(x) {
  if (x$kernel == "canonical") {
    cat("canonical kernel matrix\n")
  }
  print(x$matrix)
}

as.matrix.kernel_matrix <- function(x) x$matrix

plot.kernel_matrix <- function(x) {
  plot.df <- reshape2::melt(x$matrix)
  ggplot(plot.df, aes(Var1, Var2)) +
    geom_tile(aes(fill = value)) +
    scale_y_reverse(breaks = seq_len(nrow(x$matrix))) +
    scale_x_continuous(breaks = seq_len(ncol(x$matrix))) +
    scale_fill_gradient(low = "grey80", high = "black") +
    labs(x = NULL, y = NULL) +
    theme_minimal()
}

kern_canonical <- function(x, y = NULL) {
  # x and y must both be matrices
  if (is.null(y))  res <- tcrossprod(x)
  else res <- tcrossprod(y, x)
  res
}