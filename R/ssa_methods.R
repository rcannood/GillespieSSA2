ssa.method <- function(name, fun, diag_fun, params, initial_method_state) {
  l <- list(name = name, fun = fun, diag_fun = diag_fun, params = params, initial_method_state = initial_method_state)
  class(l) <- "fastgssa::ssamethod"
  l
}

ssa.nutiling <- function(a, nu, j) {
  M  <- ncol(nu)              # Number of reaction channels in nu-tile
  N  <- nrow(nu)              # Number of states in nu tile
  U  <- length(a) / M         # Number of tessallations of nu tile
  f  <- ceiling((j / M) - 1)  # Frameshift factor
  jp <- j - f * M             # Relative reaction channel index
  nu_jp <- nu[, jp]
  nu_j <- c(
    rep(0, f * N),            # Leading zeros
    nu_jp,                    # Relative state-change matrix
    rep(0,(U*N-(f*N+N)))      # Lagging zeros
  )
  return(nu_j)
}
ssa.direct.fun <- function(x, a, nu, method_state) {
  j    <- sample(seq_along(a), size=1, prob=a)
  nu_j <- nu[,j]
  tau  <- -log(runif(1))/sum(a)
  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state)
}
ssa.direct.diag.fun <- function(x, a, nu, method_state) {
  j    <- sample(seq_along(a), size=1, prob=a)
  nu_j <- ssa.nutiling(a, nu, j)
  tau  <- -log(runif(1))/sum(a)
  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state)
}
#' @title Direct method (D)
#'
#' @description Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).
#'
#' @return an object of to be used by \code{\link{ssa}}.
#' @seealso \link{fastgssa-package}, \code{\link{ssa}}
#' @references Gillespie (1977)
#' @keywords misc datagen ts
#' @export
ssa.direct <- function() {
  ssa.method(
    nam = "direct",
    fun = ssa.direct.fun,
    diag_fun = ssa.direct.diag.fun,
    params = c(),
    initial_method_state = list()
  )
}


ssa.etl.fun <- function(x, a, nu, method_state, tau) {
  M <- length(a)
  k <- rpois(M, (a*tau))
  # MPK: Strictly speaking it is not correct to call the realized state-change
  # vector nu_j here since, in leap methods, actually is not reaction specific.
  # In Pineda-Krch (JSS ms) it is refered to as the \bm{\tilde{nu}}
  nu_j <- (nu %*% k)[,1]
  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state)
}
ssa.etl.diag.fun <- function(x, a, nu, method_state, tau) {   # RC: did not optimise this method yet
  MU <- length(a)         # Toto nr of reaction channels
  k  <- rpois(MU,(a*tau)) # Nr of firings per channel
  M  <- ncol(nu)          # Nr of reaction channel per patch (nu)
  U  <- MU/M              # Nr of tilings
  nu_j <- NULL
  for(f in (seq(U)-1))
    nu_j <- c(nu_j, rowSums(matrix(rep(k[1:M+f*M],nrow(nu)),byrow=TRUE,ncol=M)*nu))
  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state)
}
#' @title Explicit tau-leap method (ETL)
#'
#' @description Explicit tau-leap method implementation of the \acronym{SSA} as described by Gillespie (2001).
#'
#' @param tau the step-size (default 0.3)
#' @return an object of to be used by \code{\link{ssa}}.
#' @seealso \link{fastgssa-package}, \code{\link{ssa}}
#' @references Gillespie (2001)
#' @keywords misc datagen ts
#' @export ssa.etl
ssa.etl <- function(tau = 0.3) {
  if (tau <= 0) stop(sQuote("tau"), "needs to be >0")
  ssa.method(
    name = "ETL",
    fun = ssa.etl.fun,
    diag_fun = ssa.etl.diag.fun,
    params = as.list(match.call()),
    initial_method_state = list()
  )
}


ssa.otl.fun <- function(x, a, nu, method_state, hor, nc, epsilon, dtf, nd) { # RC: did not optimise this method yet
  if (method_state$suspensions != 0) {
    method_state$suspensions <- method_state$suspensions - 1
    ssa.direct.fun(x, a, nu, method_state)
  } else {
    # 1. Identify the current critical reactions
    # Calculate the minimum number of firings for reaction before one of it's reactants becomes extinct (L). The 'L' notation is from Eq. 10 in Cao et al. (2006), J. Chem. Phys. 124:044109.
    tmp_nu <- nu
    tmp_nu[nu>=0] <- NA  # We are only interested in negative state changes

    # Turn off warning temporarily because min() throws a warning if it tries to
    # evaluate only 'NA's, which it will for reaction channels that have no
    # negative entries. Despite the warning the end result is correct, i.e. the
    # number of firings for such a channel becomes Inf.
    options(warn = -1)
    L <- apply(floor(x/abs(tmp_nu)),2,min,na.rm=TRUE)
    options(warn = 0)

    Jncr <- L >= nc                     # Indices of the non-critical reactions

    # 2. Compute the first candidate time leap, tau1
    if (sum(Jncr) == 0) {
      tau1 <- Inf                       # It is simple if there are no critical reactions present
    } else {                            # It is complicate if there are critical reactions present
      Irs <- apply((nu != 0),1,any)     # Subset the reactant species
      g <- rep(NA,length(x))
      g[hor==1]  <- 1                   # First-order reactions (S1->...)
      g[hor==2]  <- 2                   # Interspecific second-order reaction, first type (S1+S2->...)
      g[hor==22] <- (2+1/(x[hor==2]-1)) # Intraspecific second-order reaction (S1+S1->...)

      # Define mu ($\hat{\mu$}_i(\matnbf{x}) in Eq. 32a)
      tmp_nu <- matrix(nu[apply(nu, 1, function(x) any(x != 0))],ncol=dim(nu)[2]) # Remove non-reacting species from nu
      tmp <- tmp_nu[,Jncr]*a[Jncr]
      if (is.matrix(tmp)) mu <- rowSums(tmp)
      else mu <- tmp

      # Define sigma ($\hat{\sigma}^2_i(\mathbf{x})$ in Eq. 32b)
      tmp <- tmp_nu[,Jncr]^2*a[Jncr]
      if (is.matrix(tmp)) sigmaSq <- rowSums(tmp)
      else sigmaSq <- tmp

      # Calculate tau1 (Eq. 33). If there are no noncritical reactions (Jncr only
      # has FALSE elements) tau1<-Inf (see #2 in paper)
      leftTerm  <- max(epsilon*x/g,1) / abs(mu)
      rightTerm <- max(epsilon*x/g,1)^2 / abs(sigmaSq)
      tau1      <- min(leftTerm[Irs],rightTerm[Irs])
    }

    # We need to the 'while' loop with it's constructs so that we can recaulate
    # tau if we end up with negative population sizes (see step #6 in paper, page 4)
    calculateTau <- TRUE
    while (calculateTau) {
      # 3. If tau1 is "too small" return to stochRxn() and execute a number of direct method steps.
      if (tau1 < (dtf*1/sum(a))) {
        method_state$suspensions <- nd-1
        return(ssa.direct.fun(x, a, nu, method_state)) # already run a direct
        # return(list(tau = NA, nu_j = NA, j = NA, method_state = method_state))
      }

      # 4. Compute the second candidate time leap from the set of critical reactions, tau2. If there are no critical reactions tau2=Inf
      tau2 <- -log(runif(1))/sum(a[!Jncr])

      # 5. Select the actual tau from the two candidate tau (the smaller of the two)
      # and determine the number of firings each reaction will have
      if (tau1 < tau2) {                                           # Step 5a
        tau <- tau1
        k <- as.numeric(!Jncr)                                     # Sets all critical reactions to one firings and non-critical to zero firings
        k[k==0] <- rpois(sum(Jncr),(a[Jncr]*tau))                  # Sets the number of firings for non-critical reactions
      } else {                                                     # Step 5b
        tau <- tau2
        jc <- sample(seq(dim(nu)[2]),size=1,prob=(a/sum(a[!Jncr]))) # Pick one of the critical reactions that will fire once
        k <- rep(0,dim(nu)[2])                                     # Setting up an empty vector
        k[jc] <- 1                                                 # Add the selected critical reaction that is firing
        k[Jncr %in% TRUE] <- rpois(sum(Jncr),(a*tau))              # The number of firings of non-critical reactions is drawn from a Poisson distribution
      }

      # 6. Update the state vector and check for negative elements. If negative elements are found reduce
      # tau1 by half and return to step 3
      nu_j <- rowSums(matrix(rep(k,nrow(nu)),byrow=TRUE,ncol=length(a))*nu)
      if (any((x+nu_j)<0)) {
        tau1 <- tau1/2
        calculateTau <- TRUE
      } else {
        calculateTau <- FALSE
      }
    }
    method_state$suspensions <- 0
    return(list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state))
  }
}
ssa.otl.diag.fun <- function(x, a, nu, method_state, hor, nc, epsilon, dtf, nd) { # RC: did not optimise this method yet
  if (method_state$suspensions != 0) {
    method_state$suspensions <- method_state$suspensions - 1
    ssa.direct.diag.fun(x, a, nu, method_state)
  } else {
    # 1. Identify the current critical reactions
    # Calculate the minimum number of firings for reaction before one of it's
    # reactants becomes extinct (L). The 'L' notation is from Eq. 10 in Cao
    # et al. (2006). We have to turn off warning messages temporarily because
    # 'min()' throws a warning if it tries to evaluate only 'NA's, which it will
    # for reaction channels that have no negative entries. Despite the warning the
    # end result is correct, i.e. the number of firings for such a channel becomes
    # 'Inf'.
    N <- dim(nu)[1]  # Nr of states per tile
    M <- dim(nu)[2]  # Nr of reaction channels per tile
    U <- length(a)/M      # Nr of tilings
    nu_negative <- nu
    nu_negative[nu>=0] <- NA  # We are only interested in negative state changes
    L <- NULL
    options(warn = -1)
    for(f in (seq(U)-1))
      L <- c(L, apply(floor(x[1:N+f*N]/abs(nu_negative)),2,min,na.rm=TRUE))
    options(warn = 0)
    Jncr <- L >= nc  # Indices of the non-critical reactions

    # 2. Compute the first candidate time leap, tau1
    if (sum(Jncr) == 0) {
      tau1 <- Inf                       # No critical reactions present
    } else {                            # Critical reactions are present
      Irs <- rep(apply((nu != 0),1,any),U) # Subset the reactant species
      g <- rep(NA,length(x))
      g[hor==1]  <- 1                   # First-order reactions (S1->...)
      g[hor==2]  <- 2                   # Interspecific 2nd order reaction, first type (S1+S2->...)
      g[hor==22] <- (2+1/(x[hor==2]-1)) # Intraspecific 2nd order reaction (S1+S1->...)

      # Define mu ($\hat{\mu$}_i(\matnbf{x}) in Eq. 32a)
      # Define sigma ($\hat{\sigma}^2_i(\mathbf{x})$ in Eq. 32b)
      nu_reacting <- matrix(nu[apply(nu,1,any)],ncol=M) # Remove non-reacting species (i.e. rows) from nu
      mu <- NULL
      sigmaSq <- NULL
      for(f in (seq(U)-1)) {
        a_current_frame <- a[1:M+f*M]
        Jncr_current_frame <- Jncr[1:M+f*M]
        mu_tmp <- nu_reacting[,Jncr_current_frame]*a_current_frame[Jncr_current_frame]
        sigmaSq_tmp <- nu_reacting[,Jncr_current_frame]^2*a_current_frame[Jncr_current_frame]
        if (is.matrix(mu_tmp)) mu <- c(mu, rowSums(mu_tmp))
        else mu <- c(mu, mu_tmp)
        if (is.matrix(sigmaSq_tmp)) sigmaSq <- c(sigmaSq, rowSums(sigmaSq_tmp))
        else sigmaSq <- c(sigmaSq, sigmaSq_tmp)
      }

      # Calculate tau1 (Eq. 33). If there are no noncritical reactions (Jncr only
      # has FALSE elements) tau1<-Inf (see #2 in paper)
      leftTerm  <- max(epsilon*x/g,1) / abs(mu)
      rightTerm <- max(epsilon*x/g,1)^2 / abs(sigmaSq)
      tau1      <- min(leftTerm[Irs],rightTerm[Irs])
      if (is.infinite(tau1)) cat("tau1=Inf\n") # DEBUG
      if (is.na(tau1)) browser() # DEBUG
    } # if (sum(Jncr) == 0)

    # We need to the 'while' loop with it's constructs so that we can recaulate
    # tau if we end up with negative population sizes (see #6 in paper, page 4)
    calculateTau <- TRUE
    while (calculateTau) {

      # 3. If tau1 is "too small" return to stochRxn() and execute a number of
      # direct method steps.
      if (tau1 < (dtf*1/sum(a))) {
        method_state$suspensions <- nd-1
        return(ssa.direct.diasg.fun(x, a, nu, method_state)) # already run a direct
      }

      # 4. Compute the second candidate time leap from the set of critical
      # reactions, tau2. If there are no critical reactions tau2=Inf
      tau2 <- -log(runif(1))/sum(a[!Jncr])

      # 5. Select the actual tau from the two candidate tau (the smaller of the
      # two) and determine the number of firings each reaction will have
      if (tau1 < tau2) {                           # Step 5a
        tau <- tau1
        k <- as.numeric(!Jncr)                     # Sets all critical reactions to one firings and non-critical to zero firings
        lambda <- (a[Jncr]*tau) # Fudge for negative probabilities
        lambda[lambda<0] <- 0
        if (any(lambda<0)) {cat("1\n"); browser()}
        k[k==0] <- rpois(sum(Jncr),lambda)  # Sets the number of firings for non-critical reactions
      } else {                                     # Step 5b
        tau <- tau2
        pr <- (a/sum(a[!Jncr])) # Fudge for negative probabilities
        pr[pr<0] <- 0
        if (any(pr<0)) { cat("3\n"); browser()}
        jc <- sample(seq(M*U),size=1,prob=pr) # Pick one of the critical reactions that will fire once
        k <- rep(0,(M*U))                          # Setting up an empty vector
        k[jc] <- 1                                 # Add the selected critical reaction that is firing
        lambda <- (a*tau) # Fudge for negative probabilities
        lambda[lambda<0] <- 0
        if (any(lambda<0)) {cat("2\n"); browser()}
        k[Jncr %in% TRUE] <- rpois(sum(Jncr),(a*tau))  # The number of firings of non-critical reactions is drawn from a Poisson distribution
      }

      # 6. Update the state vector and check for negative elements. If negative
      # elements are found reduce tau1 by half and return to step 3

      # Update the state-change vector by nu-tiling
      nu_j <- NULL
      for(f in (seq(U)-1))
        nu_j <- c(nu_j, rowSums(matrix(rep(k[1:M+f*M],dim(nu)[1]),byrow=TRUE,ncol=M)*nu))

      if (any((x+nu_j)<0)) {
        tau1 <- tau1/2
        calculateTau <- TRUE
      } else {
        calculateTau <- FALSE
      }
    } # while
    method_state$suspensions <- 0
    return(list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state))
  }
}
#' @title Optimized tau-leap method (OTL)
#'
#' @description Optimized tau-leap method implementation of the \acronym{SSA} as described by Cao et al. (2006).
#'
#' @param x0 x0 param
#' @param hor highest order reaction vector (one entry per species in \code{x})
#' @param nc number of critical reactions threshold parameter.
#' @param epsilon error control parameter.
#' @param dtf Direct method threshold factor for temporarily suspending the method.
#' @param nd number of Direct method steps to perform during a suspension.
#' @return an object of to be used by \code{\link{ssa}}.
#' @seealso \link{fastgssa-package}, \code{\link{ssa}}
#' @references Cao et al. (2006)
#' @keywords misc datagen ts
#' @export ssa.otl
ssa.otl <- function(x0, hor = NaN, nc = 10, epsilon = 0.03, dtf = 10, nd = 100) {
  if (any(is.na(hor)))
    hor <- rep(2, length(x0)) # Undefined hor - use default values
  else if (length(hor) != length(x0))
    stop("length of hor vector is different from length of 'x0'")
  else if (any(hor!=1 & hor!=2 & hor!=22))
    stop("wrong value(s) in hor vector (can only be 1, 2, or 22)")

  ssa.method(
    name = "OTL",
    fun = ssa.otl.fun,
    diag_fun = ssa.otl.diag.fun,
    params = list(hor = hor, nc = nc, epsilon = epsilon, dtf = dtf, nd = nd),
    initial_method_state = list(
      suspensions = 0
    )
  )
}



ssa.btl.fun <- function(x, a, nu, method_state, f) { # RC: did not optimise this method yet
  coercing <- FALSE

  # Calculate tau
  tau <- f/sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- length(a)    # Number of reaction channels
  tilde_x <- x
  nu_j <- matrix(rep(0,length(x)))

  # Loop over all reaction channels having propensity fun>0
  for (j in seq(M)[a>0]) {
    if (any(nu[,j]<0)) { # do this if there are limiting reactions
      mask <- nu[,j]<0
      L <- min(floor(tilde_x[mask]/abs(nu[mask,j])))
      if (a[j]*tau>L) {
        p <- 1
        coercing <- TRUE
      } else {
        p <- a[j]*tau/L
      }
      k <- rbinom(1,L,p)
    } else { # do this if there are no limiting reactions
      k <- rpois(1,(a[j]*tau))
    }

    # Update tilde_x for the current reaction j
    tmp_nu_j <- matrix(rep(k,nrow(nu)), byrow=TRUE, ncol=1)*nu[,j]
    tilde_x <- tilde_x + tmp_nu_j

    # Record the current state change (it is returned by ssa.btl)
    nu_j <- nu_j + tmp_nu_j
  }

  if(coercing) warning("coerced p to unity - consider lowering f")

  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state)
}
ssa.btl.diag.fun <- function(x, a, nu, method_state, f) { # RC: did not optimise this method yet
  coercing <- FALSE

  # Calculate tau
  tau <- f/sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- dim(nu)[2] # Number of reaction channels per nu-tile
  N <- dim(nu)[1] # Number of states per nu-tile
  MU <- length(a)      # Toto nr of reaction channels
  U <- MU/M            # Nr of tilings
  tilde_x <- x
  nu_j <- rep(0,(N*U))

  # Loop over all reaction channels having a non-zero (>0) propensity fun
  for (j in seq(U*M)[a>0]) {
    f <- ceiling((j/M)-1)
    jp <- j-f*M  # Intra-patch reaction channel index (j->jp)
    x1 <- 1+f*N
    x2 <- 1+f*N+(N-1)
    if (any(nu[,jp]<0)) {  # Do this if there are limiting reactions
      mask <- nu[,jp]<0    # Which species has the limiting reaction
      tilde_xt <- tilde_x[x1:x2]
      L <- min(floor(tilde_xt[mask]/abs(nu[mask,jp])))
      if (a[j]*tau>L) {
        p <- 1
        coercing <- TRUE
      } else {
        p <- a[j]*tau/L
        k <- rbinom(1,L,p)
      }
    } else { # do this if there are no limiting reactions
      k <- rpois(1,(a[j]*tau))
    }

    # Update tilde_x for the current reaction j
    tmp_nu_j <- rep(k,dim(nu)[1])*nu[,jp]
    tilde_x[x1:x2] <- tilde_x[x1:x2] + tmp_nu_j

    # Record the current state change
    nu_j[x1:x2] <- nu_j[x1:x2] + tmp_nu_j
  }

  if(coercing) warning("coerced p to unity - consider lowering f")

  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state)
}
#' @title Binomial tau-leap method (BTL)
#'
#' @description Binomial tau-leap method implementation of the \acronym{SSA} as described by Chatterjee et al. (2005).
#'
#' @param f coarse-graining factor (see page 4 in Chatterjee et al. 2005).
#' @return an object of to be used by \code{\link{ssa}}.
#' @seealso \link{fastgssa-package}, \code{\link{ssa}}
#' @references Chatterjee et al. (2005)
#' @keywords misc datagen ts
#' @export ssa.btl
ssa.btl <- function(f = 10) {
  if (f <= 1) stop(sQuote("f"), " has to be >1")
  ssa.method(
    nam = "BTL",
    fun = ssa.btl.fun,
    diag_fun = ssa.btl.diag.fun,
    params = c(f = f),
    initial_method_state = list()
  )
}



ssa.em.fun <- function(x, a, nu, method_state, h, noise_strength) {
  tau <- h
  nu_j <- rowSums(t(t(nu) * a)) * h + sqrt(abs(x)) * noise_strength * rnorm(length(x), 0, h)
  # this is the fastest way according to
  # http://stackoverflow.com/questions/3643555/multiply-rows-of-matrix-by-vector

  list(tau = tau, nu_j = nu_j, j = which(nu_j != 0), method_state = method_state)
}
ssa.em.diag.fun <- function(x, a, nu, method_state, h, noise_strength) {
  error("Diag deterministic not supported yet")
}
#' @title Euler-Maruyama method (EM)
#'
#' @description Euler-Maruyama method implementation
#' @param h h parameter
#' @param noise_strength noise_strength parameter
#' @return an object of to be used by \code{\link{ssa}}.
#' @seealso \link{fastgssa-package}, \code{\link{ssa}}
#' @keywords misc datagen ts
#' @export
ssa.em <- function(h = 0.01, noise_strength = 2) {
  ssa.method(
    nam = "em",
    fun = ssa.em.fun,
    diag_fun = ssa.em.diag.fun,
    params = c(h = h, noise_strength = noise_strength),
    initial_method_state = list()
  )
}

