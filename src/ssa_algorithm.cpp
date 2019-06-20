#include <Rcpp.h>

#include "ssa.h"
#include "utils.h"

using namespace Rcpp;

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const double, NumericVector&, NumericVector&);

// [[Rcpp::export]]
List simulate(
    SEXP propensity_funs,
    SEXP ssa_alg,
    const NumericVector& initial_state,
    const NumericVector& params,
    const IntegerMatrix& nu,
    const double final_time,
    const double census_interval,
    const int buffer_size,
    const double max_walltime,
    const bool stop_on_neg_state,
    const bool verbose,
    const double console_interval,
    const bool use_vector_optimisation
) {
  // fetch propensity functions from pointer
  TR_FUN* propensity_funs_ = XPtr<TR_FUN>(propensity_funs);

  // fetch ssa algorithm from pointer
  SSA *ssa_alg_ = XPtr<SSA>(ssa_alg);

  // initialise data structures
  double simtime = 0.0;
  NumericVector state(initial_state.size());
  for (int i = 0; i < initial_state.size(); i++) {
    state[i] = initial_state[i];
  }
  NumericVector propensity(nu.ncol());
  double simtime_nextcensus = simtime + census_interval;

  // preallocate data structures
  ssa_alg_->allocate(nu.ncol(), nu.nrow());
  double dtime = 0.0;
  NumericVector dstate(state.size());
  NumericVector buffer(buffer_size);

  // check whether nu is able to be transformed into a vector
  IntegerVector nu_row(nu.ncol()), nu_effect(nu.ncol());
  bool nu_vector = use_vector_optimisation;
  fill_nu_vectors(nu, nu_row, nu_effect, &nu_vector);

  // calculate initial propensity
  for (int i = 0; i < propensity.size(); i++) {
    propensity_funs_[i](state, params, simtime, propensity, buffer);
  }

  int output_nexti = 0;
  NumericVector output_time(10);
  NumericMatrix output_state(10, state.size());
  NumericMatrix output_propensity(10, propensity.size());
  NumericMatrix output_buffer(10, buffer.size());

  output_time[output_nexti] = simtime;
  for (int i = 0; i < state.size(); i++) {
    output_state(output_nexti, i) = state[i];
  }
  for (int i = 0; i < propensity.size(); i++) {
    output_propensity(output_nexti, i) = propensity[i];
  }
  for (int i = 0; i < buffer.size(); i++) {
    output_buffer(output_nexti, i) = buffer[i];
  }
  output_nexti++;

  // track walltime
  int walltime_start = time(NULL);
  int walltime_nextconsole = walltime_start, walltime_nextinterrupt = walltime_start, walltime_curr = walltime_start;

  // verbose
  if (verbose) {
    Rcout << "Running SSA " << ssa_alg_->name << " with console output every " << console_interval << " seconds" << std::endl;
    Rcout << "Start time: " << "CURRTIME" << std::endl;
    // flush console?
  }

  bool negative_state = false;
  bool zero_prop = false;

  while (
      simtime < final_time &&
        (walltime_curr - walltime_start) < max_walltime &&
        !zero_prop &&
        (!negative_state || !stop_on_neg_state)
    )  {
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

    // make a step
    if (nu_vector) {
      ssa_alg_->step_vector(state, propensity, nu_row, nu_effect, &dtime, dstate);
    } else {
      ssa_alg_->step_matrix(state, propensity, nu, &dtime, dstate);
    }

    state += dstate;
    simtime += dtime;

    // Check that no states are negative (can occur in some tau-leaping methods)
    for (int i = 0; i < state.length(); i++) {
      if (state[i] < 0) {
        negative_state = true;
        if (!stop_on_neg_state) {
          state[i] = 0;
        }
      }
    }

    // recalculate propensity
    for (int i = 0; i < propensity.size(); i++) {
      propensity_funs_[i](state, params, simtime, propensity, buffer);
    }

    // check whether all propensity functions are zero
    zero_prop = true;
    for (int i = 0; i < propensity.size() && zero_prop; i++) {
      if (propensity[i] > 0) {
        zero_prop = false;
      }
    }

    // perform census if so desired
    if (simtime_nextcensus <= simtime) {
      simtime_nextcensus += census_interval;
      if (output_nexti == output_time.size()) {
        // output = resize(output, output.size() * 2);
        output_time = resize(output_time, output_nexti * 2);
        output_state = resize_rows(output_state, output_nexti * 2);
        output_propensity = resize_rows(output_propensity, output_nexti * 2);
        output_buffer = resize_rows(output_buffer, output_nexti * 2);
      }

      output_time[output_nexti] = simtime;
      for (int i = 0; i < state.size(); i++) {
        output_state(output_nexti, i) = state[i];
      }
      for (int i = 0; i < propensity.size(); i++) {
        output_propensity(output_nexti, i) = propensity[i];
      }
      for (int i = 0; i < buffer.size(); i++) {
        output_buffer(output_nexti, i) = buffer[i];
      }
      output_nexti++;
    }
  }

  // TODO: record end state if census_interval is set to inf

  // determine whether extinction has occurred
  bool extinction = true;
  for (int i = 0; i < state.size() && extinction; i++) {
    if (state[i] > 0) {
      extinction = false;
    }
  }

  // construct output
  int walltime_end = time(NULL);
  int walltime_elapsed = walltime_end - walltime_start;

  DataFrame stats = DataFrame::create(
    _["method"] = ssa_alg_->name,
    _["final_time_reached"] = simtime > final_time,
    _["extinction"] = extinction,
    _["negative_state"] = negative_state,
    _["zero_prop"] = zero_prop,
    _["max_walltime"] = walltime_elapsed >= max_walltime,
    _["walltime_start"] = walltime_start,
    _["walltime_end"] = walltime_end,
    _["walltime_elapsed"] = walltime_elapsed
  );

  // remove empty output slots
  output_time = resize(output_time, output_nexti);
  output_state = resize_rows(output_state, output_nexti);
  output_propensity = resize_rows(output_propensity, output_nexti);
  output_buffer = resize_rows(output_buffer, output_nexti);

  if (verbose) {
    Rcout << "SSA finished!" << std::endl;
  }

  return List::create(
    _["time"] = output_time,
    _["state"] = output_state,
    _["propensity"] = output_propensity,
    _["buffer"] = output_buffer,
    _["stats"] = stats
  );
}
