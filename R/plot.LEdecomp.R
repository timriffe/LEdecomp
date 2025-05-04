#' Function to plot the results from the LEdecomp
#' @description
#' R function to plot the age-specific life-expectancy decomposition value using one of the different method from the `LEdecomp()` function.
#'
#' @param x `x` developed using function `LEdecomp()` which are objects of the `LEdecomp` class.
#' @param ... additional arguments to show in the plot appearance.
#'
#' @return plot the different age-specific decomposition values for the selected method. This function is valid for all the methods including the sensitivities ones.
#'
#' @seealso \code{\link{LEdecomp}}
#'
#' @references
#'
#' \insertRef{arriaga1984measuring}{coddecomp}
#' \insertREf{Chandrasekaran1986}{coddecomp}
#' \insertRef{preston2000demography}{coddecomp}
#' \insertREf{Ponnapalli2005}{coddecomp}
#'
#' @importFrom ggplot ggplot
#'
#' @examples
#'
#' \donttest{
#'
#' }
#' @export
plot.LEdecomp <- function(x, ...){
  if(!is.null(x)){
    if(!"LEdecomp" %in% class(x))
      stop("The 'x' does not have the 'decomp' structure of R LEdecomp package.")
  }
  if(!is.list(x)){
    stop("The 'x' is not a list. Use 'LEdecomp' function first.")
  }

  if(!is.matrix(x$LEdecomp)){
    if(x$method == "arriaga" | x$method == "arriaga_sym" |
       x$method == "chandrasekaran_ii" | x$method == "chandrasekaran_iii" |
       x$method == "lopez_ruzicka" | x$method == "lopez_ruzicka_sym" |
       x$method == "horiuchi" | x$method == "stepwise" | x$method == "numerical"){

      data <- data.frame(age = x$age,
                         LEdecomp = x$LEdecomp)
      title <- paste(x$method, " LE decomposition method")

      barplot(names = data$age, height = data$LEdecomp, col = "lightgreen",
              xlab = "ages", ylab = "", main = title)

    } else if(x$method == "sen_arriaga" | x$method == "sen_arriaga_sym" |
              x$method == "sen_arriaga_inst" | x$method == "sen_arriaga_inst2" |
              x$method == "sen_chandrasekaran_ii" | x$method == "sen_chandrasekaran_ii_inst" | x$method == "sen_chandrasekaran_ii_inst2" |
              x$method == "sen_chandrasekaran_iii" | x$method == "sen_chandrasekaran_iii_inst" | x$method == "sen_chandrasekaran_iii_inst2" |
              x$method == "sen_lopez_ruzicka" | x$method == "sen_lopez_ruzicka_sym" |
              x$method == "sen_lopez_ruzicka_inst" | x$method == "sen_lopez_ruzicka_sym__inst2"){

      values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst", "sen_arriaga_inst2",
                  "sen_chandrasekaran_ii", "sen_chandrasekaran_ii_inst", "sen_chandrasekaran_ii_inst2",
                  "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst", "sen_chandrasekaran_iii_inst2",
                  "sen_lopez_ruzicka", "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_sym__inst2")

      names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2",
                 "chandrasekaran_ii", "chandrasekaran_ii_inst", "chandrasekaran_ii_inst2",
                 "chandrasekaran_iii", "chandrasekaran_iii_inst", "chandrasekaran_iii_inst2",
                 "lopez_ruzicka", "lopez_ruzicka_sym", "lopez_ruzicka_inst", "lopez_ruzicka_sym__inst2")
      dicc <- names[which(x$method == values)[1]]

      data <- data.frame(age = x$age,
                         LEdecomp = x$LEdecomp)
      title <- paste(dicc, " LE sensitivty analysis")

      barplot(names = data$age, height = data$LEdecomp, col = "lightgreen",
              xlab = "ages", ylab = "", main = title)

    }
  } else if(is.matrix(x$LEdecomp)){

    if(is.null(colnames(x$LEdecomp))){
      colnames(x$LEdecomp) <- paste0("c", seq_len(dim(x$LEdecomp)[2]))
    }

    if(x$method == "arriaga" | x$method == "arriaga_sym" |
       x$method == "chandrasekaran_ii" | x$method == "chandrasekaran_iii" |
       x$method == "lopez_ruzicka" | x$method == "lopez_ruzicka_sym" |
       x$method == "horiuchi" | x$method == "stepwise" | x$method == "numerical"){

      data <- data.frame(age = rep(x$age, dim(x$mx1)[2]),
                         LEdecomp = as.vector(x$LEdecomp),
                         cause = rep(colnames(x$LEdecomp), each = length(x$age)))
      title <- paste(x$method, "cause-of-death LE decomposition method")

      ggplot(data, aes(fill = cause, y = LEdecomp, x = age)) +
        geom_bar(position = "stack", stat = "identity") +
        labs(title = title) +
        theme_classic()


    } else if(x$method == "sen_arriaga" | x$method == "sen_arriaga_sym" |
              x$method == "sen_arriaga_inst" | x$method == "sen_arriaga_inst2" |
              x$method == "sen_chandrasekaran_ii" | x$method == "sen_chandrasekaran_ii_inst" | x$method == "sen_chandrasekaran_ii_inst2" |
              x$method == "sen_chandrasekaran_iii" | x$method == "sen_chandrasekaran_iii_inst" | x$method == "sen_chandrasekaran_iii_inst2" |
              x$method == "sen_lopez_ruzicka" | x$method == "sen_lopez_ruzicka_sym" |
              x$method == "sen_lopez_ruzicka_inst" | x$method == "sen_lopez_ruzicka_sym__inst2"){

      values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst", "sen_arriaga_inst2",
                  "sen_chandrasekaran_ii", "sen_chandrasekaran_ii_inst", "sen_chandrasekaran_ii_inst2",
                  "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst", "sen_chandrasekaran_iii_inst2",
                  "sen_lopez_ruzicka", "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_sym__inst2")

      names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2",
                 "chandrasekaran_ii", "chandrasekaran_ii_inst", "chandrasekaran_ii_inst2",
                 "chandrasekaran_iii", "chandrasekaran_iii_inst", "chandrasekaran_iii_inst2",
                 "lopez_ruzicka", "lopez_ruzicka_sym", "lopez_ruzicka_inst", "lopez_ruzicka_sym__inst2")
      dicc <- names[which(x$method == values)[1]]

      data <- data.frame(age = rep(x$age, dim(x$mx1)[2]),
                         LEdecomp = as.vector(x$LEdecomp),
                         cause = rep(colnames(x$LEdecomp), each = length(x$age)))
      title <- paste(dicc, "cause-of-death sensitivity Life-Expectancy decomposition method.")

      ggplot(data, aes(fill = cause, y = LEdecomp, x = age)) +
        geom_bar(position = "stack", stat = "identity") +
        labs(title = title) +
        theme_classic()

    }

  }

}
