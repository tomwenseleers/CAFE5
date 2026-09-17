# Independent R generator-matrix check. Only base R is required.
# Usage: Rscript check_kernel.R /path/to/cafe5
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
binary <- normalizePath(args[1L], mustWork = TRUE)
lambda <- 0.2; nu <- 0.3; t <- 0.7; K <- 80L
Q <- matrix(0, K + 1L, K + 1L)
for (n in 0:K) {
  if (n < K) Q[n + 1L, n + 2L] <- lambda * n + nu
  if (n > 0) Q[n + 1L, n] <- lambda * n
}
diag(Q) <- -rowSums(Q)
# Independent uniformization: exp(tQ)=sum Poisson(m;omega*t)(I+Q/omega)^m.
omega <- max(-diag(Q))
B <- diag(K + 1L) + Q / omega
power <- diag(K + 1L)
reference <- dpois(0, omega*t) * power
for (m in seq_len(qpois(1-1e-14, omega*t))) {
  power <- power %*% B
  reference <- reference + dpois(m, omega*t) * power
}
f <- tempfile(fileext = '.tsv')
status <- system2(binary, c('--innovation', '--lambda', lambda, '--nu', nu,
  '--matrix-time', t, '--max-count', K, '--matrix-output', shQuote(f)))
stopifnot(status == 0L)
actual <- as.matrix(read.table(f, header=FALSE))
unlink(f)
# Compare low states, far from the generator's reflecting upper boundary.
error <- max(abs(actual[1:16,1:16] - reference[1:16,1:16]))
stopifnot(error < 1e-12)
cat(sprintf('R_UNIFORMIZATION_OK max_abs_error=%.17g\n', error))
