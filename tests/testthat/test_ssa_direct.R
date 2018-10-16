#' parms <- c(c = 0.5)
#' initial.state  <- c(X = 10000)
#' a   <- c("c*X")
#' nu  <- matrix(-1)
#' out <- ssa(initial.state,a,nu,parms,final.time = 10,simName = "Irreversible isomerization") # Direct method
#' plot(out$data[,1],out$data[,2]/10000,col = "red",cex = 0.5,pch = 19)
#'
#' a <- runif(100)
#' nu <- matrix(runif(100 * length(a)), ncol = length(a))
#' j <- sample.int(length(a), 1)
#' fastgssa:::ssa.nutiling(a, nu, j)
#' fastgssa:::ssa_nutiling(a, nu, j)
#' ssa_nutiling(a, nu, j)
#'
#'
#' Rcpp::cppFunction("NumericVector ssa_nutiling(NumericVector a, NumericMatrix nu, int j) {
#'   int M = nu.ncol();                  // Number of reaction channels in nu-tile
#'   int N = nu.nrow();                  // Number of states in nu tile
#'   int U = a.length() / M;             // Number of tessallations of nu tile
#'   int f = ceil((j / M) - 1);          // Frameshift factor
#'   int jp = j - f * M;                 // Relative reaction channel index
#'   NumericVector nu_j = rep(0.0, U*N);
#'   // for (int i = 0; i < N; i++) {
#'   //   nu_j(i + f * N) = nu(i, jp);
#'   // }
#'   return nu_j;                        // Return output
#' }")
#'
#'
#'
#' ssa.nutiling <- function(a, nu, j) {
#'   M  <- ncol(nu)              # Number of reaction channels in nu-tile
#'   N  <- nrow(nu)              # Number of states in nu tile
#'   U  <- length(a) / M         # Number of tessallations of nu tile
#'   f  <- ceiling((j / M) - 1)  # Frameshift factor
#'   jp <- j - f * M             # Relative reaction channel index
#'   nu_jp <- nu[, jp]
#'   nu_j <- c(
#'     rep(0, f * N),            # Leading zeros
#'     nu_jp,                    # Relative state-change matrix
#'     rep(0,(U*N-(f*N+N)))      # Lagging zeros
#'   )
#'   return(nu_j)
#' }
#'
#'
#'
#' #' #' ## Setup
#' #' num.as <- 20
#' #' num.runs <- 10000
#' #'
#' #' #' lower bound was approximated by calculating the cor for varying num runs
#' #' lowerbound.cor <- 1 - 10^(2-log10(num.runs))
#' #'
#' #' nu <- diag(runif(num.as))
#' #' a <- runif(num.as)
#' #' a <- a / sum(a)
#' #'
#' #'
#' #'
#' #' #' ## Tests
#' #' context("Direct SSA in R")
#' #' outs <- lapply(seq_len(num.runs), function(i) {
#' #'   fastgssa::ssa.d(a, nu)
#' #' })
#' #'
#' #' js <- sapply(outs, function(o) o$j)
#' #' emperical.freq <- sapply(seq_len(num.as), function(ai) mean(js == ai))
#' #'
#' #' test_that("weighted sampling", {
#' #'   expect_gte(cor(a, emperical.freq), lowerbound.cor)
#' #' })
#' #'
#' #' test_that("good selection of nu_j", {
#' #'   nus <- do.call(rbind, lapply(outs, function(o) o$nu_j))
#' #'   manual.nus <- do.call(rbind, lapply(outs, function(o) nu[,o$j]))
#' #'   expect_equal(nus, manual.nus)
#' #' })
#' #'
#' #'
#' #' context("Direct SSA in Rcpp")
#' #' outs <- lapply(seq_len(num.runs), function(i) {
#' #'   fastgssa::ssa_d(a, nu)
#' #' })
#' #'
#' #' js <- sapply(outs, function(o) o$j)
#' #' emperical.freq <- sapply(seq_len(num.as), function(ai) mean(js == ai))
#' #'
#' #' test_that("weighted sampling", {
#' #'   expect_gte(cor(a, emperical.freq), lowerbound.cor)
#' #' })
#' #'
#' #' test_that("good selection of nu_j", {
#' #'   nus <- do.call(rbind, lapply(outs, function(o) o$nu_j))
#' #'   manual.nus <- do.call(rbind, lapply(outs, function(o) nu[,o$j]))
#' #'   expect_equal(nus, manual.nus)
#' #' })
