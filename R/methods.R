do_this <- FALSE
if (do_this){
  method_registry <- tibble::tribble(
    ~method,                        ~fun_name,                               ~category,    ~pkg,
    "lifetable",                    "sen_e0_mx_lt",                          "opt_ok",     NA,
    "arriaga",                      "arriaga",                               "direct",     NA,
    "arriaga_sym",                  "arriaga_sym",                           "direct",     NA,
    "sen_arriaga",                  "sen_arriaga",                           "direct_sen", NA,
    "sen_arriaga_sym",             "sen_arriaga_sym",                        "direct_sen", NA,
    "sen_arriaga_inst",            "sen_arriaga_instantaneous",             "opt_ok",     NA,
    "sen_arriaga_inst2",           "sen_arriaga_instantaneous2",            "opt_ok",     NA,
    "sen_arriaga_sym_inst",        "sen_arriaga_sym_instantaneous",         "opt_ok",     NA,
    "sen_arriaga_sym_inst2",       "sen_arriaga_sym_instantaneous2",        "opt_ok",     NA,
    "chandrasekaran_ii",           "chandrasekaran_II",                      "direct",     NA,
    "sen_chandrasekaran_ii",       "sen_chandrasekaran_II",                  "direct_sen", NA,
    "sen_chandrasekaran_ii_inst",  "sen_chandrasekaran_II_instantaneous",   "opt_ok",     NA,
    "sen_chandrasekaran_ii_inst2", "sen_chandrasekaran_II_instantaneous2",  "opt_ok",     NA,
    "chandrasekaran_iii",          "chandrasekaran_III",                     "direct",     NA,
    "sen_chandrasekaran_iii",      "sen_chandrasekaran_III",                 "direct_sen", NA,
    "sen_chandrasekaran_iii_inst", "sen_chandrasekaran_III_instantaneous",  "opt_ok",     NA,
    "sen_chandrasekaran_iii_inst2","sen_chandrasekaran_III_instantaneous2", "opt_ok",     NA,
    "lopez_ruzicka",               "lopez_ruzicka",                          "direct",     NA,
    "lopez_ruzicka_sym",           "lopez_ruzicka_sym",                      "direct",     NA,
    "sen_lopez_ruzicka",           "sen_lopez_ruzicka",                      "direct_sen", NA,
    "sen_lopez_ruzicka_sym",       "sen_lopez_ruzicka_sym",                  "direct_sen", NA,
    "sen_lopez_ruzicka_inst",      "sen_lopez_ruzicka_instantaneous",       "opt_ok",     NA,
    "sen_lopez_ruzicka_inst2",     "sen_lopez_ruzicka_instantaneous2",      "opt_ok",     NA,
    "sen_lopez_ruzicka_sym_inst",  "sen_lopez_ruzicka_sym_instantaneous",   "opt_ok",     NA,
    "sen_lopez_ruzicka_sym_inst2", "sen_lopez_ruzicka_sym_instantaneous2",  "opt_ok",     NA,
    "numerical",                   "sen_num",                                "opt_ok",     NA,
    "stepwise",                    "stepwise_replacement",                   "general",    "DemoDecomp",
    "horiuchi",                    "horiuchi",                               "general",    "DemoDecomp"
  ) |>
    dplyr::mutate(pkg = ifelse(is.na(pkg),"LEdecomp","pkg"))

usethis::use_data(method_registry, internal = FALSE, overwrite = TRUE)

}

#' Registry of Available Life Expectancy Decomposition Methods
#'
#' A reference table listing all decomposition methods implemented in the `LEdecomp` package,
#' along with their corresponding function names and classification.
#'
#' This registry helps centralize method metadata and supports internal operations
#' like method matching, class-based routing, and printing.
#'
#' @format A data frame with 28 rows and 3 variables:
#' \describe{
#'   \item{method}{Character. The method name used in the `method` argument of `LEdecomp()`.}
#'   \item{fun_name}{Character. The actual function name (as a string) used to compute the decomposition.}
#'   \item{category}{Character. One of `"direct"`, `"direct_sen"`, `"opt_ok"`, or `"general"`, indicating how the method operates internally:
#'     \itemize{
#'       \item `"direct"`: Classic decomposition methods using two full input vectors (e.g., Arriaga, Lopez-Ruzicka).
#'       \item `"direct_sen"`: Sensitivity-based methods that take \code{mx1} and \code{mx2}, and divide the result by their difference.
#'       \item `"opt_ok"`: Sensitivity methods that accept a single \code{mx} input and optionally allow optimization over the interpolation point between \code{mx1} and \code{mx2}.
#'       \item `"general"`: Generic methods like \code{stepwise} and \code{horiuchi} using externally provided tools (e.g., from DemoDecomp).
#'     }
#'   }
#'   \item{pkg}{Character. The package where the underlying decomp function resides.}
#' }
#'
#' @seealso \code{\link{available_methods}}, \code{\link{LEdecomp}}
#' @keywords datasets internal
"method_registry"



#' List available decomposition methods
#'
#' Returns a table of all implemented methods, their function name, and category.
#'
#' @return A data frame of available decomposition methods.
#' @export
#' @examples
#' available_methods()
available_methods <- function() {
  method_registry
}


get_dec_fun <- function(method) {
  row <- method_registry[method_registry$method == method, ]
  if (!is.na(row$pkg)){
    dec_fun <- getExportedValue( row$pkg, row$fun_name)
  } else {
    dec_fun <- get(row$fun_name)
  }
  dec_fun
}
