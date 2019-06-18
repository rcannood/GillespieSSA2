#ifndef DYNGEN_SSA_ALGORITHM_H
#define DYNGEN_SSA_ALGORITHM_H

#include <Rcpp.h>
#include "ssa.hpp"
#include "ssa_em.hpp"
#include "ssa_direct.hpp"
#include "ssa_etl.hpp"
#include "ssa_btl.hpp"

using namespace Rcpp;

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const double, NumericVector&);

// [[Rcpp::export]]
List simulate(
    SEXP transition_fun,
    SEXP ssa_alg,
    const NumericVector& initial_state,
    const NumericVector& params,
    const IntegerMatrix& nu,
    const double final_time,
    const double census_interval,
    const double max_walltime,
    const bool stop_on_neg_state,
    const bool verbose,
    const double console_interval
) {
  List output(10);

  // fetch ssa algorithm from pointer
  TR_FUN transition_fun_ = *XPtr<TR_FUN>(transition_fun);
  SSA *ssa_alg_ = XPtr<SSA>(ssa_alg);

  // initialise data structures
  double simtime = 0.0;
  NumericVector state = clone(initial_state);
  NumericVector transition_rates(nu.ncol());
  double simtime_nextcensus = simtime + census_interval;

  // preallocate data structures
  ssa_alg_->allocate(nu.ncol(), nu.nrow());
  double dtime = 0.0;
  NumericVector dstate(state.size());

  // check whether nu is filled with single values
  IntegerVector nu_row(nu.ncol()), nu_effect(nu.ncol());
  bool nu_single = true;
  for (int j = 0; j < nu.ncol() && nu_single; j++) {
    for (int i = 0; i < nu.nrow(); i++) {
      if (nu(i, j) != 0) {
        if (nu_effect[j] == 0) {
          nu_effect[j] = nu(i, j);
          nu_row[j] = i;
        } else {
          nu_single = false;
        }
      }
    }
  }

  // calculate initial transition rates
  transition_fun_(state, params, simtime, transition_rates);

  output(0) = List::create(
    _["time"] = simtime,
    _["state"] = clone(state),
    _["transition_rates"] = clone(transition_rates)
  );
  int output_nexti = 1;

  // track walltime
  int walltime_start = time(NULL);
  int walltime_nextconsole = walltime_start, walltime_nextinterrupt = walltime_start, walltime_curr = walltime_start;

  // verbose
  if (verbose) {
    Rcout << "Running SSA " << ssa_alg_->name << " with console output every " << console_interval << " seconds" << std::endl;
    Rcout << "Start time: " << "CURRTIME" << std::endl;
    // flush console?
  }

  while (simtime < final_time && (walltime_curr - walltime_start) <= max_walltime)  {
    walltime_curr = time(NULL);

    // check for interrupt
    if (walltime_nextinterrupt <= walltime_curr) {
      checkUserInterrupt();
      walltime_nextinterrupt += 1;
    }

    // print if so desired
    if (verbose && walltime_nextconsole <= walltime_curr) {
      Rcout << "walltime: " << (walltime_curr - walltime_start) << ", simtime: " << simtime << std::endl;
      walltime_nextconsole += console_interval;
    }

    // make a transition step
    if (nu_single) {
      ssa_alg_->step_single(state, transition_rates, nu_row, nu_effect, &dtime, dstate);
    } else {
      ssa_alg_->step(state, transition_rates, nu, &dtime, dstate);
    }

    state += dstate;
    simtime += dtime;

    // Check that no states are negative (can occur in some tau-leaping methods)
    for (int i = 0; i < state.length(); i++) {
      if (state[i] < 0) {
        if (stop_on_neg_state) {
          stop("state vector contains negative value at position " + i);
        } else {
          state[i] = 0;
        }
      }
    }

    // recalculate transition rates
    transition_fun_(state, params, simtime, transition_rates);

    // perform census if so desired
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

  // construct output
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
  //   Rcout << "Stats:" << std::endl;
  //   Rcout << stats << std::endl;
  // }

  return List::create(
    _["output"] = output,
    _["stats"] = stats
  );
}

#endif
