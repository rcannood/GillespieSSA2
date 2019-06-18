#' @importFrom stringr str_count str_replace_all str_extract_all str_replace
#' @importFrom RcppXPtrUtils cppXPtr
#' @importFrom dynutils safe_tempdir
#' @export
compile_propensity_functions <- function(
  propensity_funs,
  state,
  params,
  reuse_buffer = TRUE,
  hardcode_params = FALSE,
  env = parent.frame()
) {
  buffer_usages <- str_count(propensity_funs, "=")
  buffer_size <-
    if (reuse_buffer) {
      max(buffer_usages)
    } else {
      sum(buffer_usages)
    }

  buffer_offset <- 0
  buffer_names <- c()

  rcpp_prop_funs <- map_chr(
    propensity_funs,
    function(prop_fun) {
      buff_nam <- prop_fun %>% str_extract_all("[A-Za-z_0-9]* *=") %>% first() %>% str_replace_all("[ =]", "")
      buffer_names <<- c(buffer_names, buff_nam)
      prop_split <- gsub("([A-Za-z][A-Za-z0-9_]*)", " \\1 ", prop_fun) %>% strsplit(" ") %>% first() %>% discard(~ .=="")

      state_match <- match(prop_split, names(state))
      state_ix <- which(!is.na(state_match))
      prop_split[state_ix] <- paste0("state[", state_match[state_ix] - 1, "]")

      params_match <- match(prop_split, names(params))
      params_ix <- which(!is.na(params_match))
      if (hardcode_params) {
        prop_split[params_ix] <- paste0("CONST_", prop_split[params_ix])
        # prop_split[params_ix] <- params[params_match[params_ix]]
      } else {
        prop_split[params_ix] <- paste0("params[", params_match[params_ix] - 1, "]")
      }

      buffer_match <- match(prop_split, buff_nam)
      buffer_ix <- which(!is.na(buffer_match))
      prop_split[buffer_ix] <- paste0("buffer[", buffer_match[buffer_ix] - 1 + buffer_offset, "]")

      if (!reuse_buffer) {
        buffer_offset <<- buffer_offset + length(buff_nam)
      }

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
    "  NumericVector& transition_rates,\n",
    "  NumericVector& buffer\n",
    ") {\n",
    paste(paste0(buffers, calculations), collapse = ""),
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

  pointer <- RcppXPtrUtils::cppXPtr(rcpp_code, cacheDir = tmpdir, includes = includes)

  # TODO: remove debug code:
  write_lines(paste0(includes, "\n", rcpp_code), "~/fastgssa.cpp")

  l <- lst(
    buffer_names,
    buffer_size,
    pointer,
    reuse_buffer,
    hardcode_params
  )
  class(l) <- "fastgssa::propensity_functions"
  l
}
