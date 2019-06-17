#' @importFrom stringr str_count str_replace_all str_extract_all str_replace
#' @importFrom RcppXPtrUtils cppXPtr
#' @export
compile_propensity_functions <- function(propensity_funs, state, params, env = parent.frame()) {
  buffer_size <- max(str_count(propensity_funs, "="))
  rcpp_prop_funs <- map_chr(
    propensity_funs,
    function(prop_fun) {
      buffer_names <- prop_fun %>% str_extract_all("[A-Za-z_0-9]* *=") %>% first() %>% str_replace_all("[ =]", "")
      prop_split <- gsub("([A-Za-z][A-Za-z0-9_]*)", " \\1 ", prop_fun) %>% strsplit(" ") %>% first() %>% discard(~ .=="")

      state_match <- match(prop_split, names(state))
      state_ix <- which(!is.na(state_match))
      prop_split[state_ix] <- paste0("state[", state_match[state_ix] - 1, "]")

      params_match <- match(prop_split, names(params))
      params_ix <- which(!is.na(params_match))
      # prop_split[params_ix] <- paste0("params[", params_match[params_ix] - 1, "]")
      prop_split[params_ix] <- params[params_match[params_ix]]

      buffer_match <- match(prop_split, buffer_names)
      buffer_ix <- which(!is.na(buffer_match))
      prop_split[buffer_ix] <- paste0("buffer[", buffer_match[buffer_ix] - 1, "]")

      paste(prop_split, collapse = "")
    }
  )

  buffers <- rcpp_prop_funs %>% str_replace(";[^=]*$", ";") %>% {ifelse(grepl("=", .), ., "")} %>% str_replace_all("([^;]*;)", "  \\1\n")
  calculations <- rcpp_prop_funs %>% str_replace("^.*;", "") %>% {paste0("  transition_rates[", seq_along(.)-1, "] = ", ., ";\n")}

  rcpp_code <- paste0(
    "void calculate_transition_rates(\n",
    "  const NumericVector& state,\n",
    "  const NumericVector& params,\n",
    "  const double time,\n",
    "  NumericVector& transition_rates\n",
    ") {\n",
    ifelse(buffer_size == 0, "", paste0("    NumericVector buffer = no_init(", buffer_size, ");\n")),
    paste(paste0(buffers, calculations), collapse = ""),
    "}\n"
  )
  tmpdir <- dynutils::safe_tempdir("fastgssa")
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))

  transition_functon <- RcppXPtrUtils::cppXPtr(rcpp_code, cacheDir = tmpdir)

  transition_functon
}
