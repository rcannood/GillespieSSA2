#include <Rcpp.h>
#include <math.h>

#include "ssa.h"
#include "utils.h"

using namespace Rcpp;

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const double, NumericVector&, NumericVector&);

// [[Rcpp::export]]
List simulate(
    SEXP propensity_funs,
    const int num_functions,
    SEXP ssa_alg,
    const NumericVector& initial_state,
    const NumericVector& params,
    const IntegerVector& nu_i,
    const IntegerVector& nu_p,
    const IntegerVector& nu_x,
    const double final_time,
    const double census_interval,
    const int buffer_size,
    const double max_walltime,
    const bool stop_on_neg_state,
    const bool verbose,
    const double console_interval,
    const CharacterVector& sim_name
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
  NumericVector propensity(nu_p.size() - 1);
  double simtime_nextcensus = simtime + census_interval;

  // fields for storing statistics
  double dtime_mean = 0;
  double dtime_sd = 0;
  int dtime_steps = 0;

  // preallocate data structures
  ssa_alg_->allocate(propensity.size(), state.size());
  double dtime = 0.0;
  NumericVector dstate(state.size());
  NumericVector buffer(buffer_size);

  // calculate initial propensity
  for (int i = 0; i < num_functions; i++) {
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

    // clear dstate and dtime
    for (int i = 0; i < dstate.size(); i++) {
      dstate[i] = 0;
    }
    dtime = 0;

    // make a step
    ssa_alg_->step(state, propensity, nu_i, nu_p, nu_x, &dtime, dstate);

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
    for (int i = 0; i < num_functions; i++) {
      propensity_funs_[i](state, params, simtime, propensity, buffer);
    }

    // check whether all propensity functions are zero
    zero_prop = true;
    for (int i = 0; i < propensity.size() && zero_prop; i++) {
      if (propensity[i] > 0) {
        zero_prop = false;
      }
    }

    // update statistics
    dtime_steps++;
    if (dtime_steps == 1) {
      dtime_sd = 0;
    } else {
      dtime_sd = sqrt((dtime_steps - 2) / (dtime_steps - 1) * pow(dtime_sd, 2) + pow(dtime - dtime_mean, 2) / dtime_steps);
    }
    dtime_mean = (dtime_mean * (dtime_steps - 1) + dtime) / dtime_steps;

    // perform census if so desired
    if (simtime_nextcensus <= simtime) {
      simtime_nextcensus += census_interval;
      if (output_nexti == output_time.size()) {
        output_time = gillespie::resize_vector(output_time, output_nexti * 2);
        output_state = gillespie::resize_rows(output_state, output_nexti * 2);
        output_propensity = gillespie::resize_rows(output_propensity, output_nexti * 2);
        output_buffer = gillespie::resize_rows(output_buffer, output_nexti * 2);
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

  // end state if census_interval is set to inf
  if (std::isinf(census_interval)) {
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
    _["sim_name"] = sim_name,
    _["stop_simtime"] = simtime > final_time,
    _["stop_extinction"] = extinction,
    _["stop_negative_state"] = negative_state,
    _["stop_zero_prop"] = zero_prop,
    _["stop_walltime"] = walltime_elapsed >= max_walltime,
    _["walltime_start"] = walltime_start,
    _["walltime_end"] = walltime_end,
    _["walltime_elapsed"] = walltime_elapsed,
    _["num_steps"] = dtime_steps,
    _["dtime_mean"] = dtime_mean,
    _["dtime_sd"] = dtime_sd
  );

  // remove empty output slots
  output_time = gillespie::resize_vector(output_time, output_nexti);
  output_state = gillespie::resize_rows(output_state, output_nexti);
  output_propensity = gillespie::resize_rows(output_propensity, output_nexti);
  output_buffer = gillespie::resize_rows(output_buffer, output_nexti);

  if (verbose) {
    Rcout << "SSA finished!" << std::endl;
  }

  return List::create(
    _["time"] = output_time,
    _["state"] = output_state,
    _["propensity"] = output_propensity,
    _["buffer"] = output_buffer,
    _["stats"] = stats,
    _["sim_name"] = sim_name
  );
}
