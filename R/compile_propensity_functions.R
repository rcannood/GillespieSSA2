#' Precompile the SSA propensity functions
#'
#' By precompiling the propensity functions, you can run multiple SSA simulations with the
#' same propensity functions without having to recompile the functions every time.
#'
#' See [ssa()] for more information on the format of each of the parameters.
#'
#' @param propensity_funs \[character\] The propensity functions.
#' @param state_ids \[character\] The names of the states in the correct order.
#' @param params \[named numeric\] Constants that are used in the propensity functions.
#' @param buffer_ids \[character, optional\] The order of any buffer calculations that are made as part of the propensity functions.
#' @param hardcode_params \[logical\] Whether or not to hard code the parameters into the source code.
#'
#' @importFrom stringr str_count str_replace_all str_extract_all str_replace
#' @importFrom Rcpp sourceCpp
#' @importFrom dynutils safe_tempdir
#' @importFrom dplyr last
#' @export
compile_propensity_functions <- function(
  propensity_funs,
  state_ids,
  params,
  buffer_ids = NULL,
  hardcode_params = FALSE
) {
  assert_that(
    is.character(propensity_funs)
  )
  # should this be re-enabled?
  reuse_buffer <- FALSE

  # preprocess propensity functions
  variable_names <- propensity_funs %>% str_extract_all("[A-Za-z_0-9]* *=") %>% map(~ str_replace_all(., "[ =]", ""))
  reaction_ids <- variable_names %>% map_chr(last)

  buffer_usages <- map_int(variable_names, length) - 1
  buffer_size <-
    if (reuse_buffer) {
      max(buffer_usages)
    } else {
      sum(buffer_usages)
    }

  if (is.null(buffer_ids)) {
    buffer_ids <- variable_names %>% unlist() %>% discard(~ . %in% reaction_ids)
  }

  # parse propensity functions
  parse <- map_df(
    seq_along(propensity_funs),
    function(i) {
      prop_fun <- propensity_funs[[i]]

      prop_split <- gsub("([A-Za-z][A-Za-z0-9_]*)", " \\1 ", prop_fun) %>% strsplit(" ") %>% first() %>% discard(~ .== "")

      state_match <- match(prop_split, state_ids)
      state_ix <- which(!is.na(state_match))
      prop_split[state_ix] <- paste0("state[", state_match[state_ix] - 1, "]")

      params_match <- match(prop_split, names(params))
      params_ix <- which(!is.na(params_match))
      if (hardcode_params) {
        prop_split[params_ix] <- paste0("CONST_", prop_split[params_ix])
      } else {
        prop_split[params_ix] <- paste0("params[", params_match[params_ix] - 1, "]")
      }

      buffer_match <- match(prop_split, buffer_ids)
      buffer_ix <- which(!is.na(buffer_match))
      if (length(buffer_ix) > 0) {
        buffer_to_ix <- buffer_match[buffer_ix] - 1
        if (reuse_buffer) {
          buffer_to_ix <- match(buffer_to_ix, unique(buffer_to_ix)) - 1
        }
        prop_split[buffer_ix] <- paste0("buffer[", buffer_to_ix, "]")
      }

      reaction_match <- match(prop_split, reaction_ids)
      reaction_ix <- which(!is.na(reaction_match))
      prop_split[reaction_ix] <- paste0("propensity[", reaction_match[reaction_ix] - 1, "]")

      tibble(
        rcpp = paste(c(prop_split, ";"), collapse = ""),
        states = list(unique(prop_split[state_ix])),
        params = list(unique(prop_split[params_ix])),
        buffers = list(unique(prop_split[buffer_ix])),
        reaction = unique(prop_split[reaction_ix])
      )
    }
  )

  # wrap into separate functions
  rcpp_codes <- paste0(
    "void calculate_propensity_", seq_along(parse$rcpp) - 1, "(\n",
    "  const NumericVector& state,\n",
    "  const NumericVector& params,\n",
    "  const double time,\n",
    "  NumericVector& propensity,\n",
    "  NumericVector& buffer\n",
    ") {\n",
    parse$rcpp %>% str_replace_all("([^;]*;)", "  \\1\n"),
    "}\n"
  )

  # create code to be able to return TR functions
  rcpp_code <- paste0(
    "#include <Rcpp.h>\n",
    "using namespace Rcpp;\n",
    if (hardcode_params) paste0("#define CONST_", names(params), " ", params, "\n", collapse = "") else character(),
    "typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const double, NumericVector&, NumericVector&);\n",
    "\n",
    paste(rcpp_codes, collapse = "\n"),
    "\n",
    "// [[Rcpp::export]]\n",
    "SEXP return_functions() {\n",
    "  TR_FUN *functions = new TR_FUN[", length(rcpp_codes), "];\n",
    paste0("  functions[", seq_along(rcpp_codes) - 1, "] = &calculate_propensity_", seq_along(rcpp_codes) - 1, ";\n", collapse = ""),
    "  XPtr<TR_FUN> ptr(functions);\n",
    "  return ptr;\n",
    "}\n"
  )

  # create temporary dir for compilation
  tmpdir <- dynutils::safe_tempdir("fastgssa")
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))

  # compile code
  return_functions <- NULL # sourceCpp will override this
  Rcpp::sourceCpp(code = rcpp_code, cacheDir = tmpdir, env = environment())

  # return propensity functions as pointer
  pointer <- return_functions()

  # return output
  l <- lst(
    buffer_ids,
    buffer_size,
    pointer,
    reuse_buffer,
    hardcode_params
  )
  class(l) <- "fastgssa::propensity_functions"
  l
}
