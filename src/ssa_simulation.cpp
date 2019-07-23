#include <Rcpp.h>
#include <math.h>

#include "ssa_method.h"
#include "utils.h"

using namespace Rcpp;

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const double, NumericVector&, NumericVector&);

class SSA_simulation {
public:
  SSA_simulation() {}

  /*
  ~SSA_simulation() {}
  */

  void initialise(
    // propensity functions
    int num_functions_,
    SEXP propensity_funs_,
    NumericVector& params_,
    int buffer_size_,
    // algorithm
    SEXP ssa_method_,
    // state
    NumericVector& initial_state_,
    // state change matrix
    IntegerVector& nu_i_,
    IntegerVector& nu_p_,
    IntegerVector& nu_x_,
    // output parameters
    double census_interval_,
    bool log_propensity_,
    bool log_firings_,
    bool log_buffer_,
    // stopping conditions
    bool stop_on_neg_state_,
    double final_time_,
    double max_walltime_,
    // meta information
    std::string sim_name_,
    bool verbose_,
    double console_interval_
  ) {
    // process prop funs
    num_functions = num_functions_;
    prop_funs = XPtr<TR_FUN>(propensity_funs_);
    params = params_;
    buffer = NumericVector(buffer_size_);

    // process algorithm
    ssa_alg = XPtr<SSA_method>(ssa_method_);

    // process state
    initial_state = initial_state_;
    state = NumericVector(initial_state.size());
    dstate = NumericVector(initial_state.size());

    // process state change matrix
    nu_i = nu_i_;
    nu_p = nu_p_;
    nu_x = nu_x_;

    propensity = NumericVector(nu_p.size() - 1);
    firings = NumericVector(propensity.size());
    dfirings = NumericVector(propensity.size());

    // process output parameters
    census_interval = census_interval_;
    log_propensity = log_propensity_;
    log_firings = log_firings_;
    log_buffer = log_buffer_;

    // process stopping conditions
    stop_on_neg_state = stop_on_neg_state_;
    final_time = final_time_;
    max_walltime = max_walltime_;

    // meta information
    sim_name = sim_name_;
    verbose = verbose_;
    console_interval = console_interval_;

    // initialise output structures
    output_nexti = 0;
    output_time = NumericVector(10);
    output_state = NumericMatrix(10, state.size());
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

  List run() {
    // reset output strunctures
    output_nexti = 0;
    resize_outputs(10, true);

    // initialise walltime fields
    int walltime_start = time(NULL);
    int walltime_nextconsole = walltime_start, walltime_nextinterrupt = walltime_start, walltime_curr = walltime_start;

    // initialise state data structures
    sim_time = 0;
    dtime = 0;
    double sim_time_nextcensus = census_interval;
    std::copy(initial_state.begin(), initial_state.end(), state.begin()) ;
    std::fill(dstate.begin(), dstate.end(), 0);
    std::fill(buffer.begin(), buffer.end(), 0);
    std::fill(firings.begin(), firings.end(), 0);
    std::fill(dfirings.begin(), dfirings.end(), 0);
    calculate_propensity();

    // statistics data structures
    num_steps = 0;
    dtime_mean = 0;
    dtime_sd = 0;
    firings_mean = 0;
    firings_sd = 0;

    negative_state = false;
    zero_prop = false;

    // do census of initial state
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
    resize_outputs(output_nexti, false);

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
    // enlarge output objects, if needed
    if (output_nexti == output_time.size()) {
      resize_outputs(output_nexti * 2, false);
    }

    // copy state to output objects
    output_time[output_nexti] = sim_time;
    output_state(output_nexti, _) = state;
    if (output_propensity.nrow() > 0) {
      output_propensity(output_nexti, _) = propensity;
    }
    if (output_buffer.nrow() > 0) {
      output_buffer(output_nexti, _) = buffer;
    }
    if (output_firings.nrow() > 0) {
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

    // Check that no states are negative (can occur in some tau-leaping methods)
    for (auto i = state.begin(); i != state.end(); ++i) {
      if (*i < 0) {
        if (!stop_on_neg_state) {
          *i = 0;
        }
        negative_state = true;
      }
    }
  }

  void resize_outputs(int num_rows, bool fill_zero) {
    output_time = resize_vector(output_time, num_rows, !fill_zero);
    output_state = resize_matrix(output_state, num_rows, output_state.ncol(), !fill_zero);
    if (output_propensity.nrow() > 0) {
      output_propensity = resize_matrix(output_propensity, num_rows, output_propensity.ncol(), !fill_zero);
    }
    if (output_buffer.nrow() > 0) {
      output_buffer = resize_matrix(output_buffer, num_rows, output_buffer.ncol(), !fill_zero);
    }
    if (output_firings.nrow() > 0) {
      output_firings = resize_matrix(output_firings, num_rows, output_firings.ncol(), !fill_zero);
    }

    if (fill_zero) {
      std::fill(output_time.begin(), output_time.end(), 0);
      std::fill(output_state.begin(), output_state.end(), 0);
      std::fill(output_propensity.begin(), output_propensity.end(), 0);
      std::fill(output_buffer.begin(), output_buffer.end(), 0);
      std::fill(output_firings.begin(), output_firings.end(), 0);
    }
  }

  template <typename T>
  T resize_vector(const T& x, int n, bool copy){
    int oldsize = x.size();
    if (n == oldsize) {
      return x;
    } else if (n < oldsize) {
      oldsize = n;
    }
    T y(n); // TODO: no init?
    if (copy) {
      for( int i = 0; i < oldsize; i++) {
        y[i] = x[i];
      }
    }
    return y;
  }

  template <typename T>
  T resize_matrix(const T& x, int nr, int nc, bool copy){
    int oldnr = x.nrow();
    int oldnc = x.ncol();
    if (nr == oldnr && nc == oldnc) {
      return x;
    }
    if (nr < oldnr) {
      oldnr = nr;
    }
    if (nc < oldnc) {
      oldnc = nc;
    }
    T y(nr, nc); // TODO: no init?
    if (copy) {
      for( int i = 0; i < oldnr; i++) {
        for (int j = 0; j < oldnc; j++) {
          y(i, j) = x(i, j);
        }
      }
    }
    return y;
  }

  int num_functions;
  TR_FUN *prop_funs;
  SSA_method *ssa_alg;
  NumericVector initial_state;
  NumericVector params;
  IntegerVector nu_i;
  IntegerVector nu_p;
  IntegerVector nu_x;

  // state data structures
  double sim_time;
  double dtime;
  NumericVector state;
  NumericVector dstate;
  NumericVector propensity;
  NumericVector buffer;
  NumericVector firings;
  NumericVector dfirings;

  // statistics data structures
  int num_steps;
  double dtime_mean;
  double dtime_sd;
  double firings_mean;
  double firings_sd;

  // log data structures
  int output_nexti;
  NumericVector output_time;
  NumericMatrix output_state;
  NumericMatrix output_propensity;
  NumericMatrix output_buffer;
  NumericMatrix output_firings;

  // output parameters
  double census_interval;
  bool log_propensity;
  bool log_firings;
  bool log_buffer;

  // stopping conditions
  bool zero_prop;
  bool negative_state;
  bool stop_on_neg_state;
  double final_time;
  double max_walltime;

  // meta information
  std::string sim_name;
  bool verbose;
  double console_interval;
};


// all of these fields are made publicly accessible to allow unit testing
RCPP_MODULE(gillespie) {
  class_<SSA_simulation>("SSA_simulation")
  .constructor()
  .method("initialise", &SSA_simulation::initialise)
  .method("run", &SSA_simulation::run)
  .method("do_census", &SSA_simulation::do_census)
  .method("calculate_propensity", &SSA_simulation::calculate_propensity)
  .method("make_step", &SSA_simulation::make_step)
  .method("resize_outputs", &SSA_simulation::resize_outputs)
  .field("num_functions", &SSA_simulation::num_functions)
  .field("initial_state", &SSA_simulation::initial_state)
  .field("params", &SSA_simulation::params)
  .field("nu_i", &SSA_simulation::nu_i)
  .field("nu_p", &SSA_simulation::nu_p)
  .field("nu_x", &SSA_simulation::nu_x)
  .field("sim_time", &SSA_simulation::sim_time)
  .field("dtime", &SSA_simulation::dtime)
  .field("state", &SSA_simulation::state)
  .field("dstate", &SSA_simulation::dstate)
  .field("propensity", &SSA_simulation::propensity)
  .field("buffer", &SSA_simulation::buffer)
  .field("firings", &SSA_simulation::firings)
  .field("dfirings", &SSA_simulation::dfirings)
  .field("num_steps", &SSA_simulation::num_steps)
  .field("dtime_mean", &SSA_simulation::dtime_mean)
  .field("dtime_sd", &SSA_simulation::dtime_sd)
  .field("firings_mean", &SSA_simulation::firings_mean)
  .field("firings_sd", &SSA_simulation::firings_sd)
  .field("output_nexti", &SSA_simulation::output_nexti)
  .field("output_time", &SSA_simulation::output_time)
  .field("output_state", &SSA_simulation::output_state)
  .field("output_propensity", &SSA_simulation::output_propensity)
  .field("output_buffer", &SSA_simulation::output_buffer)
  .field("output_firings", &SSA_simulation::output_firings)
  .field("census_interval", &SSA_simulation::census_interval)
  .field("log_propensity", &SSA_simulation::log_propensity)
  .field("log_firings", &SSA_simulation::log_firings)
  .field("log_buffer", &SSA_simulation::log_buffer)
  .field("zero_prop", &SSA_simulation::zero_prop)
  .field("negative_state", &SSA_simulation::negative_state)
  .field("stop_on_neg_state", &SSA_simulation::stop_on_neg_state)
  .field("final_time", &SSA_simulation::final_time)
  .field("max_walltime", &SSA_simulation::max_walltime)
  .field("sim_name", &SSA_simulation::sim_name)
  .field("verbose", &SSA_simulation::verbose)
  .field("console_interval", &SSA_simulation::console_interval)
  ;
}


// [[Rcpp::export]]
List test_ssa_step(
    SEXP ssa_alg,
    NumericVector& state,
    NumericVector& propensity,
    IntegerVector& nu_i,
    IntegerVector& nu_p,
    IntegerVector& nu_x
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
