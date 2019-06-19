#' @importFrom stringr str_count str_replace_all str_extract_all str_replace
#' @importFrom RcppXPtrUtils cppXPtr
#' @importFrom dynutils safe_tempdir
#' @export
compile_propensity_functions <- function(
  propensity_funs,
  reaction_ids,
  state,
  params,
  reuse_buffer = TRUE,
  hardcode_params = FALSE,
  buffer_ids = NULL,
  env = parent.frame()
) {
  variable_names <- propensity_funs %>% str_extract_all("[A-Za-z_0-9]* *=") %>% map(~ str_replace_all(., "[ =]", ""))
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

  rcpp_prop_funs <- map_chr(
    seq_along(propensity_funs),
    function(i) {
      prop_fun <- propensity_funs[[i]]

      prop_split <- gsub("([A-Za-z][A-Za-z0-9_]*)", " \\1 ", prop_fun) %>% strsplit(" ") %>% first() %>% discard(~ .== "")

      state_match <- match(prop_split, names(state))
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

      paste(c(prop_split, ";"), collapse = "")
    }
  )

  rcpp_code <- paste0(
    "void calculate_propensity(\n",
    "  const NumericVector& state,\n",
    "  const NumericVector& params,\n",
    "  const double time,\n",
    "  NumericVector& propensity,\n",
    "  NumericVector& buffer\n",
    ") {\n",
    rcpp_prop_funs %>% str_replace_all("([^;]*;)", "  \\1\n") %>% paste(collapse = ""),
    "}\n"
  )

  includes <-
    if (hardcode_params) {
      paste0(
        "#define CONST_", names(params), " ", params,
        collapse = "\n"
      )
    } else {
      character()
    }

  tmpdir <- dynutils::safe_tempdir("fastgssa")
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))

  # TODO: remove debug code:
  readr::write_lines(paste0(includes, "\n", rcpp_code), "~/fastgssa.cpp")

  pointer <- RcppXPtrUtils::cppXPtr(rcpp_code, cacheDir = tmpdir, includes = includes)

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
