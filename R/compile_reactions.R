#' Precompile the reactions
#'
#' By precompiling the reactions, you can run multiple SSA simulations repeatedly
#' without having to recompile the reactions every time.
#'
#' @param reactions '[reaction]' A list of multiple [reaction()] objects.
#' @param state_ids `[character]` The names of the states in the correct order.
#' @param params `[named numeric]` Constants that are used in the propensity functions.
#' @param buffer_ids `[character]` The order of any buffer calculations that are made as part of the propensity functions.
#' @param hardcode_params `[logical]` Whether or not to hardcode the values of `params` in the compilation of the propensity functions.
#'   Setting this to `TRUE` will result in a minor sacrifice in accuracy for a minor increase in performance.
#' @param write_rcpp `[character]` If not `NULL`, then the source code of the propensity functions will be written
#'   to this file location before compilation.
#' @param fun_by `[integer]` Combine this number of propensity functions into one function.
#'
#' @return A list of objects solely to be used by [ssa()].
#'
#' * `x[["state_change"]]`: A sparse matrix of reaction effects.
#' * `x[["reaction_ids"]]`: The names of the reactions.
#' * `x[["buffer_ids"]]`: A set of buffer variables found in the propensity functions.
#' * `x[["buffer_size"]]`: The minimum size of the buffer required.
#' * `x[["functions_pointer"]]`: A pointer to the compiled propensity functions.
#' * `x[["num_functions"]]`: The compiled propensity functions are split up into multiple batches of functions.
#' * `x[["hardcode_params"]]`: Whether the parameters were hard coded into the source code.`
#'
#' @importFrom stringr str_count str_replace_all str_extract_all str_replace str_split
#' @importFrom Rcpp sourceCpp
#' @importFrom dynutils safe_tempdir %all_in%
#' @importFrom dplyr first last
#' @importFrom readr write_lines
#'
#' @export
compile_reactions <- function(
  reactions,
  state_ids,
  params,
  buffer_ids = NULL,
  hardcode_params = FALSE,
  write_rcpp = NA_character_,
  fun_by = 100L
) {
  assert_that(is.list(reactions))
  walk(seq_along(reactions), function(i) {
    assert_that(
      is(reactions[[i]], "SSA_reaction"),
      names(reactions[[i]]$effect) %all_in% state_ids
    )
  })
  reaction_ids <- map_chr(reactions, function(reac) reac$name) %|% paste0("reaction", seq_along(reactions))

  # check ids
  assert_that(
    is.character(state_ids),
    !any(duplicated(state_ids)),
    is.numeric(params),
    length(params) == 0 || !is.null(names(params)),
    !any(duplicated(buffer_ids)),
    is_scalar_logical(hardcode_params),
    is_scalar_character(write_rcpp),
    is_scalar_integer(fun_by),
    !any(duplicated(reaction_ids))
  )

  # create nu sparse matrix
  state_change_df <- map_df(seq_along(reactions), function(j) {
    reac <- reactions[[j]]

    data.frame(
      i = match(names(reac$effect), state_ids),
      j = j,
      x = reac$effect
    )
  })
  state_change <- Matrix::sparseMatrix(
    i = state_change_df$i,
    j = state_change_df$j,
    x = state_change_df$x
  )

  # add reaction ids assignments to propensity functions
  propensity_funs <- map_chr(reactions, "propensity")
  for (i in seq_along(propensity_funs)) {
    propensity_funs[[i]] <- gsub("(.*;)?([^;]*)", paste0("\\1", reaction_ids[[i]], " = \\2"), propensity_funs[[i]])
  }

  # preprocess propensity functions
  variable_names <- propensity_funs %>% str_extract_all("[A-Za-z_0-9]* *=") %>% map(function(x) str_replace_all(x, "[ =]", ""))
  reaction_ids <- variable_names %>% map_chr(last)

  buffer_usages <- map_int(variable_names, length) - 1
  buffer_size <- sum(buffer_usages)

  if (is.null(buffer_ids)) {
    buffer_ids <- variable_names %>% unlist() %>% discard(~ . %in% reaction_ids)
  }

  # split propensity functions into words and non-words
  prop_split <-
    paste0(propensity_funs, ";", collapse = "") %>%
    str_replace_all("([A-Za-z][A-Za-z0-9_]*)", " \\1 ") %>%
    str_split(" ") %>%
    first() %>%
    discard(~ . == "")

  # substitute state variables
  state_match <- match(prop_split, state_ids)
  state_ix <- which(!is.na(state_match))
  prop_split[state_ix] <- paste0("state[", state_match[state_ix] - 1, "]")

  # substitute parameters
  params_match <- match(prop_split, names(params))
  params_ix <- which(!is.na(params_match))
  if (hardcode_params) {
    prop_split[params_ix] <- paste0("CONST_", prop_split[params_ix])
  } else {
    prop_split[params_ix] <- paste0("params[", params_match[params_ix] - 1, "]")
  }

  # substitute buffer variables
  buffer_match <- match(prop_split, buffer_ids)
  buffer_ix <- which(!is.na(buffer_match))
  if (length(buffer_ix) > 0) {
    buffer_to_ix <- buffer_match[buffer_ix] - 1
    prop_split[buffer_ix] <- paste0("buffer[", buffer_to_ix, "]")
  }

  # substitute reaction ids
  reaction_match <- match(prop_split, reaction_ids)
  reaction_ix <- which(!is.na(reaction_match))
  prop_split[reaction_ix] <- paste0("propensity[", reaction_match[reaction_ix] - 1, "]")

  # get indices of semi colons
  prop_lines <-
    paste0(prop_split, collapse = "") %>%
    str_split(";") %>%
    first() %>%
    discard(~ . == "")
  prop_lines <- paste0("  ", prop_lines, ";\n")
  forms_start <- seq(1, length(prop_lines), by = fun_by)
  forms_end <- pmin(forms_start + fun_by - 1, length(prop_lines))
  rcpp_function_code_blocks <-
    map_chr(seq_along(forms_start), function(i) {
      paste0(prop_lines[seq(forms_start[[i]], forms_end[[i]])], collapse = "")
    })

  # wrap into separate functions
  rcpp_codes <- paste0(
    "void calculate_propensity_", seq_along(rcpp_function_code_blocks) - 1, "(\n",
    "  const NumericVector& state,\n",
    "  const NumericVector& params,\n",
    "  const double time,\n",
    "  NumericVector& propensity,\n",
    "  NumericVector& buffer\n",
    ") {\n",
    rcpp_function_code_blocks,
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
    "}\n",
    "\n",
    "// [[Rcpp::export]]\n",
    "int num_functions() {\n",
    "  return ", length(rcpp_codes), ";\n",
    "}\n"
  )

  # create temporary dir for compilation
  tmpdir <- dynutils::safe_tempdir("gillespie")
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))

  # compile code
  if (!is.na(write_rcpp)) {
    readr::write_lines(rcpp_code, write_rcpp)
  }
  return_functions <- NULL # sourceCpp will override this
  Rcpp::sourceCpp(code = rcpp_code, cacheDir = tmpdir, env = environment())

  # return propensity functions as pointer
  functions_pointer <- return_functions()
  num_functions <- num_functions()

  # return output
  l <- list(
    state_change = state_change,
    reaction_ids = reaction_ids,
    buffer_ids = buffer_ids,
    buffer_size = buffer_size,
    functions_pointer = functions_pointer,
    num_functions = num_functions,
    hardcode_params = hardcode_params
  )
  class(l) <- "SSA_reactions"
  l
}
