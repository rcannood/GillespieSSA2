#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// //[[Rcpp::depends(RcppArmadillo)]]
//
// using namespace Rcpp;
//
// //' Direct method (D)
// //'
// //' Direct method implementation of the \acronym{SSA} as described by Gillespie
// //' (1977). It is usually called from within \code{\link{ssa}}, but can be
// //' invoked directly.
// //'
// //' Performs one time step using the Direct method.
// //'
// //' @param a vector of evaluated propensity functions.
// //' @param nu state-change matrix.
// //' @return A list with two elements, 1) the time leap (\code{tau}) and 2) the
// //' realized state change vector (\code{nu_j}).
// //' @seealso \link{fastgssa-package}, \code{\link{ssa}}
// //' @references Gillespie (1977)
// //' @keywords misc datagen ts
// //'
// //' @import RcppArmadillo
// //' @useDynLib fastgssa
// //' @export
// //' @examples
// //'
// //' ## Logistic growth model
// //' a = function(parms,x){
// //'  b <- parms[1]
// //'  d <- parms[2]
// //'  K <- parms[3]
// //'  N <- x[1]
// //'  return(c(b*N , N*b + (b-d)*N/K))
// //' }
// //' parms <- c(2,1,1000,500)
// //' x <- 500
// //' nu <- matrix(c(+1, -1),ncol=2)
// //' t <- 0
// //' for (i in seq(100)) {
// //'   out <- ssa_direct(NULL, a(parms, x), nu)
// //'   x <- x + out$nu_j
// //'   t <- t + 1
// //'   cat("t:",t,", x:",x,"\n")
// //' }
// //'
// //[[Rcpp::export]]
// List ssa_direct(NumericVector x, NumericVector a, NumericMatrix nu, List method_state) {
//   // generate a random index j weighted by a
//   //CharacterVector ret = RcppArmadillo::sample(x, size, replace, prob) ;
//   //RcppArmadillo::
//   NumericVector csa = cumsum(a);
//   double sum_a = sum(a);
//   double rw = runif(1, 0.0, sum_a)[0];
//   int j = 0;
//   while (csa[j] < rw) {
//     j++;
//   }
//
//   // take column nu[,j]
//   NumericVector nu_j = nu(_, j);
//
//   // generate tau
//   double tau = -log(runif(1)[0]) / sum_a;
//
//   // return output
//   return List::create(
//     Named("tau") = tau,
//     Named("nu_j") = nu_j,
//     Named("j") = j+1,
//     Named("method_state") = method_state
//   );
// }
//
// // [[Rcpp::export]]
// List ssa_direct_diag(NumericVector x, NumericVector a, NumericMatrix nu, List method_state) {
//   // generate a random index j weighted by a
//   NumericVector csa = cumsum(a);
//   double sum_a = sum(a);
//   double rw = runif(1, 0.0, sum_a)[0];
//   int j = 0;
//   while (csa[j] < rw) {
//     j++;
//   }
//
//   // take column nu[,j]
//   NumericVector nu_j = nu(_, j);
//
//   // generate tau
//   double tau = -log(runif(1)[0]) / sum_a;
//
//   // return output
//   return List::create(
//     Named("tau") = tau,
//     Named("nu_j") = nu_j,
//     Named("j") = j+1,
//     Named("method_state") = method_state
//   );
// }
//
//
//
// // ssa.direct.fun <- function(x, a, nu, method.state) {
// //   j    <- sample(seq_along(a), size=1, prob=a)
// //   nu_j <- nu[,j]
// //   tau  <- -log(runif(1))/sum(a)
// //   list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method.state = method.state)
// // }
// // ssa.direct.diag.fun <- function(x, a, nu, method.state) {
// //   j    <- sample(seq_along(a), size=1, prob=a)
// //   nu_j <- ssa.nutiling(a, nu, j)
// //   tau  <- -log(runif(1))/sum(a)
// //   list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method.state = method.state)
// // }
