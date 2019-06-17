#ifndef DYNGEN_SSA_ALGORITHM_H
#define DYNGEN_SSA_ALGORITHM_H

#include <Rcpp.h>
#include "ssa.hpp"
#include "ssa_em.hpp"
#include "ssa_direct.hpp"
#include <limits>

using namespace Rcpp;

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, double, NumericVector&);

// [[Rcpp::export]]
List simulate(
    SEXP transition_fun,
    SEXP ssa_alg,
    const NumericVector& initial_state,
    const NumericVector& params,
    const NumericMatrix& nu,
    const double final_time,
    const double census_interval,
    const double max_walltime,
    const bool stop_on_neg_state,
    const bool verbose,
    const double console_interval
) {
  List output(10);

  TR_FUN transition_fun_ = *XPtr<TR_FUN>(transition_fun);
  SSA *ssa_alg_ = XPtr<SSA>(ssa_alg);

  double simtime = 0.0;
  NumericVector state = clone(initial_state);
  NumericVector transition_rates(nu.ncol());

  // calculate initial transition rates
  transition_fun_(state, params, simtime, transition_rates);

  output(0) = List::create(
    _["time"] = simtime,
    _["state"] = clone(state),
    _["transition_rates"] = clone(transition_rates)
  );
  int output_nexti = 1;

  int walltime_start = time(NULL);
  int walltime_nextconsole = walltime_start, walltime_nextinterrupt = walltime_start, walltime_curr = walltime_start;

  double simtime_nextcensus = simtime + census_interval;

  if (verbose) {
    Rcout << "Running SSA " << ssa_alg_->name << " with console output every " << console_interval << " seconds" << std::endl;
    Rcout << "Start time: " << "CURRTIME" << std::endl;
  }
  /*if (verbose) {
   cat("Running ", method$name, " method with console output every ", console.interval, " time step\n", sep = "")
   cat("Start wall time: ", format(time.start), "\n" , sep = "")
   utils::flush.console()
  }*/

  while (simtime < final_time && (walltime_curr - walltime_start) <= max_walltime)  {
    walltime_curr = time(NULL);

    if (walltime_nextinterrupt <= walltime_curr) {
      checkUserInterrupt();
      walltime_nextinterrupt += 1;
    }

    if (verbose && walltime_nextconsole <= walltime_curr) {
      Rcout << walltime_curr << " | time = " << simtime << " : " << "STATE" << std::endl;
      walltime_nextconsole += console_interval;
    }

    double dtime = 0.0;
    NumericVector dstate(state.size());
    // TODO: pass time and state directly to step, instead of dtime and dstate?

    ssa_alg_->step(state, transition_rates, nu, &dtime, dstate);
    // step(state, transition_rates, nu, dtime, dstate);

    state += dstate;
    simtime += dtime;

    /*# Check that no states are negative (can occur in some tau-leaping methods)
     invalid_ix <- is.na(state) | state < 0
     if (any(invalid_ix)) {
     if (stop_on_neg_state) {
     stop("state vector contains negative values\n", paste(names(state)[invalid_ix], collapse = ", "))
     } else {
     state[invalid_ix] <- 0
     }
     }*/

    transition_fun_(state, params, simtime, transition_rates);

    if (simtime_nextcensus <= simtime) {
      simtime_nextcensus += census_interval;
      if (output_nexti == output.size()) {
        output = resize(output, output.size() * 2);
      }

      output(output_nexti) = List::create(
        _["time"] = simtime,
        _["state"] = clone(state),
        _["transition_rates"] = clone(transition_rates)
      );
      output_nexti++;
    }
  }

  // Display the last time step on the console
  if (verbose) {
    Rcout << "time = " << simtime << " : " << "STATE" << std::endl;
  }

  int walltime_end = time(NULL);
  int walltime_elapsed = walltime_end - walltime_start;

  DataFrame stats = DataFrame::create(
    _["method"] = ssa_alg_->name,
    _["final_time_reached"] = simtime >= final_time,
    // _["extinction"] = all(state == 0),
    // _["negative_state"] = any(state < 0),
    // _["zero_prop"] = all(transition_rates == 0),
    _["max_walltime"] = walltime_elapsed >= max_walltime,
    _["walltime_start"] = walltime_start,
    _["walltime_end"] = walltime_end,
    _["walltime_elapsed"] = walltime_elapsed
  );
  // if (verbose) {
  //   //print(stats)
  // }

  return List::create(
    _["output"] = output,
    _["stats"] = stats
  );
  // return(output);
}

#endif
