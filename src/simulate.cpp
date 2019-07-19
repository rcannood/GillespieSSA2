#include <Rcpp.h>
#include <math.h>

#include "ssa_method.h"
#include "utils.h"

using namespace Rcpp;

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const double, NumericVector&, NumericVector&);

class SSA_simulation {
public:
  SSA_simulation(
    const int num_functions,
    SEXP propensity_funs,
    SEXP ssa_method,
    const NumericVector& initial_state,
    const NumericVector& params,
    const IntegerVector& nu_i,
    const IntegerVector& nu_p,
    const IntegerVector& nu_x,
    const double final_time,
    const double census_interval,
    const bool stop_on_neg_state,
    const int buffer_size,
    const std::string sim_name,
    const double max_walltime,
    const bool log_propensity,
    const bool log_firings,
    const bool log_buffer,
    const bool verbose,
    const double console_interval
  ) :
  num_functions(num_functions),
  initial_state(initial_state),
  params(params),
  nu_i(nu_i),
  nu_p(nu_p),
  nu_x(nu_x),
  final_time(final_time),
  census_interval(census_interval),
  stop_on_neg_state(stop_on_neg_state),
  buffer_size(buffer_size),
  sim_name(sim_name),
  max_walltime(max_walltime),
  log_propensity(log_propensity),
  log_firings(log_firings),
  log_buffer(log_buffer),
  verbose(verbose),
  console_interval(console_interval)
  {
    prop_funs = XPtr<TR_FUN>(propensity_funs);
    ssa_alg = XPtr<SSA_method>(ssa_method);

    sim_time_nextcensus = census_interval;

    // initialise state
    state = clone(initial_state);
    dstate = NumericVector(state.size());

    // initialise propensity and buffer
    propensity = NumericVector(nu_p.size() - 1);
    buffer = NumericVector(buffer_size);

    // initialise firings
    firings = NumericVector(propensity.size());
    dfirings = NumericVector(propensity.size());

    // // save time in log
    output_time = NumericVector(10);
    output_state = NumericMatrix(10, state.size());

    // save propensity in log, if desired
    if (log_propensity) {
      output_propensity = NumericMatrix(10, propensity.size());
    } else {
      output_propensity = NumericMatrix(0, 0);
    }
    if (log_buffer) {
      output_buffer = NumericMatrix(10, buffer.size());
    } else {
      output_buffer = NumericMatrix(0, 0);
    }
    if (log_firings) {
      output_firings = NumericMatrix(10, firings.size());
    } else {
      output_firings = NumericMatrix(0, 0);
    }
  }

  ~SSA_simulation() {
    // todo: free all objects again?
  }

  List run() {
    // track walltime
    int walltime_start = time(NULL);
    int walltime_nextconsole = walltime_start, walltime_nextinterrupt = walltime_start, walltime_curr = walltime_start;

    // calculate initial propensities
    calculate_propensity();
    do_census();

    // verbose
    if (verbose) {
      Rcout << "Running SSA " << ssa_alg->name << " with console output every " << console_interval << " seconds" << std::endl;
      Rcout << "Start time: " << "CURRTIME" << std::endl;
    }

    while (
        sim_time < final_time &&
          (walltime_curr - walltime_start) < max_walltime &&
          !zero_prop &&
          (!negative_state || !stop_on_neg_state)
    )  {

      // check for interrupt
      if (walltime_nextinterrupt <= walltime_curr) {
        checkUserInterrupt();
        walltime_nextinterrupt += 1;
      }

      // print if so desired
      if (verbose && walltime_nextconsole <= walltime_curr) {
        Rcout << "walltime: " << (walltime_curr - walltime_start) << ", sim_time: " << sim_time << std::endl;
        walltime_nextconsole += console_interval;
      }

      // make a step
      make_step();

      // Check that no states are negative (can occur in some tau-leaping methods)
      for (auto i = state.begin(); i != state.end(); ++i) {
        if (*i < 0) {
          if (!stop_on_neg_state) {
            *i = 0;
          }
          negative_state = true;
        }
      }

      // recalculate propensity
      calculate_propensity();

      // perform census if so desired
      if (sim_time_nextcensus <= sim_time) {
        sim_time_nextcensus += census_interval;
        do_census();
      }

      walltime_curr = time(NULL);
    }

    // log end state if census_interval is set to inf
    if (std::isinf(census_interval)) {
      do_census();
    }

    // determine whether extinction has occurred
    bool extinction = true;
    for (auto i = state.begin(); i != state.end() && extinction; i++) {
      if (*i > 0) {
        extinction = false;
      }
    }

    // construct output
    int walltime_end = time(NULL);
    int walltime_elapsed = walltime_end - walltime_start;

    DataFrame stats = DataFrame::create(
      _["method"] = ssa_alg->name,
      _["sim_name"] = sim_name,
      _["stop_sim_time"] = sim_time > final_time,
      _["stop_extinction"] = extinction,
      _["stop_negative_state"] = negative_state,
      _["stop_zero_prop"] = zero_prop,
      _["stop_walltime"] = walltime_elapsed >= max_walltime,
      _["walltime_start"] = walltime_start,
      _["walltime_end"] = walltime_end,
      _["walltime_elapsed"] = walltime_elapsed,
      _["num_steps"] = num_steps,
      _["dtime_mean"] = dtime_mean,
      _["dtime_sd"] = dtime_sd,
      _["firings_mean"] = firings_mean,
      _["firings_sd"] = firings_sd
    );

    // remove empty output slots
    output_time = resize_vector(output_time, output_nexti);
    output_state = resize_rows(output_state, output_nexti);
    if (log_propensity) {
      output_propensity = resize_rows(output_propensity, output_nexti);
    }
    if (log_buffer) {
      output_buffer = resize_rows(output_buffer, output_nexti);
    }
    if (log_firings) {
      output_firings = resize_rows(output_firings, output_nexti);
    }

    if (verbose) {
      Rcout << "SSA finished!" << std::endl;
    }

    return List::create(
      _["time"] = output_time,
      _["state"] = output_state,
      _["propensity"] = output_propensity,
      _["firings"] = output_firings,
      _["buffer"] = output_buffer,
      _["stats"] = stats,
      _["sim_name"] = sim_name
    );
  }

  void do_census() {
    if (output_nexti == output_time.size()) {
      output_time = resize_vector(output_time, output_nexti * 2);
      output_state = resize_rows(output_state, output_nexti * 2);
      if (log_propensity) {
        output_propensity = resize_rows(output_propensity, output_nexti * 2);
      }
      if (log_buffer) {
        output_buffer = resize_rows(output_buffer, output_nexti * 2);
      }
      if (log_firings) {
        output_firings = resize_rows(output_firings, output_nexti * 2);
      }
    }
    output_time[output_nexti] = sim_time;
    output_state(output_nexti, _) = state;
    if (log_propensity) {
      output_propensity(output_nexti, _) = propensity;
    }
    if (log_buffer) {
      output_buffer(output_nexti, _) = buffer;
    }
    if (log_firings) {
      output_firings(output_nexti, _) = firings;
      std::fill(firings.begin(), firings.end(), 0);
    }
    output_nexti++;
  }

  void calculate_propensity() {
    for (int i = 0; i < num_functions; i++) {
      prop_funs[i](state, params, sim_time, propensity, buffer);
    }

    // check whether all propensity functions are zero
    zero_prop = true;
    for (int i = 0; i < propensity.size() && zero_prop; i++) {
      if (propensity[i] > 0) {
        zero_prop = false;
      }
    }
  }

  void make_step() {
    // clear dtime, dstate and dfirings
    dtime = 0;
    std::fill(dstate.begin(), dstate.end(), 0);
    std::fill(dfirings.begin(), dfirings.end(), 0);

    // make a step
    ssa_alg->step(state, propensity, nu_i, nu_p, nu_x, &dtime, dstate, dfirings);

    num_steps++;
    sim_time += dtime;
    state += dstate;
    firings += dfirings;

    // update statistics
    int firings_sum = sum(dfirings);
    if (num_steps == 1) {
      dtime_sd = 0;
      firings_sd = 0;
    } else {
      dtime_sd = sqrt((num_steps - 2) / (num_steps - 1) * pow(dtime_sd, 2) + pow(dtime - dtime_mean, 2) / num_steps);
      firings_sd = sqrt((num_steps - 2) / (num_steps - 1) * pow(firings_sd, 2) + pow(firings_sum - firings_mean, 2) / num_steps);
    }
    dtime_mean = (dtime_mean * (num_steps - 1) + dtime) / num_steps;
    firings_mean = (firings_mean * (num_steps - 1) + firings_sum) / num_steps;
  }

  template <typename T>
  T resize_vector(const T& x, int n){
    int oldsize = x.size();
    if (n < oldsize) {
      oldsize = n;
    }
    T y(n);
    for( int i = 0; i < oldsize; i++) {
      y[i] = x[i];
    }
    return y;
  }

  template <typename T>
  T resize_rows(const T& x, int n){
    int oldsize = x.nrow();
    if (n < oldsize) {
      oldsize = n;
    }
    T y(n, x.ncol());
    for( int i = 0; i < oldsize; i++) {
      for (int j = 0; j < x.ncol(); j++) {
        y(i, j) = x(i, j);
      }
    }
    return y;
  }

private:
  const int num_functions;
  TR_FUN *prop_funs;
  SSA_method *ssa_alg;
  const NumericVector& initial_state;
  const NumericVector& params;
  const IntegerVector& nu_i;
  const IntegerVector& nu_p;
  const IntegerVector& nu_x;
  const double final_time;
  const double census_interval;
  const bool stop_on_neg_state;
  const int buffer_size;
  const std::string sim_name;
  const double max_walltime;
  const bool log_propensity;
  const bool log_firings;
  const bool log_buffer;
  const bool verbose;
  const double console_interval;

  // state data structures
  double sim_time = 0;
  double dtime = 0;
  double sim_time_nextcensus;
  NumericVector state;
  NumericVector dstate;
  NumericVector propensity;
  NumericVector buffer;
  NumericVector firings;
  NumericVector dfirings;

  // statistics data structures
  int num_steps = 0;
  double dtime_mean = 0;
  double dtime_sd = 0;
  double firings_mean = 0;
  double firings_sd = 0;

  bool negative_state = false;
  bool zero_prop = false;

  // log data structures
  int output_nexti = 0;
  NumericVector output_time;
  NumericMatrix output_state;
  NumericMatrix output_propensity;
  NumericMatrix output_buffer;
  NumericMatrix output_firings;
} ;

// here comes the boilerplate
// [[Rcpp::export]]
List simulate(
  const int num_functions,
  SEXP propensity_funs,
  SEXP ssa_method,
  const NumericVector& initial_state,
  const NumericVector& params,
  const IntegerVector& nu_i,
  const IntegerVector& nu_p,
  const IntegerVector& nu_x,
  const double final_time,
  const double census_interval,
  const bool stop_on_neg_state,
  const int buffer_size,
  const std::string sim_name,
  const double max_walltime,
  const bool log_propensity,
  const bool log_firings,
  const bool log_buffer,
  const bool verbose,
  const double console_interval
) {
  SSA_simulation *sim = new SSA_simulation(
    num_functions,
    propensity_funs,
    ssa_method,
    initial_state,
    params,
    nu_i,
    nu_p,
    nu_x,
    final_time,
    census_interval,
    stop_on_neg_state,
    buffer_size,
    sim_name,
    max_walltime,
    log_propensity,
    log_firings,
    log_buffer,
    verbose,
    console_interval
  );
  List out = sim->run();
  // delete sim;
  return out;
}


// [[Rcpp::export]]
List test_ssa_step(
    SEXP ssa_alg,
    const NumericVector& state,
    const NumericVector& propensity,
    const IntegerVector& nu_i,
    const IntegerVector& nu_p,
    const IntegerVector& nu_x
) {
  SSA_method *ssa_alg_ = XPtr<SSA_method>(ssa_alg);
  double dtime = 0;
  NumericVector dstate(state.size());
  NumericVector firings(propensity.size());
  ssa_alg_->step(state, propensity, nu_i, nu_p, nu_x, &dtime, dstate, firings);
  return List::create(
    _["dtime"] = dtime,
    _["dstate"] = dstate,
    _["firings"] = firings
  );
}
