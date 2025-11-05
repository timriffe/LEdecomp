#' Plot Life-Expectancy Decomposition Results (ggplot2)
#'
#' Plot contributions (or sensitivities) to a life expectancy difference by age or by age and cause using ggplot2. This is just for a quick default plot method.
#' @details By default, if `what = "LEdecomp"` we plot using bars (`geom = "bar"`), but you can override this. For bar plots, recall it's the area, not the height that the eye reads; for this reason, if your data is in non-single ages, we divide out the interval width, so that the implied uniform graduation to single ages still sums to the gap.
#'
#' @param x An object of class `"LEdecomp"`.
#' @param what Which series to plot: `"LEdecomp"` (default) or `"sens"`.
#' @param col Optional vector of colors for causes. If `NULL`, a fixed palette
#'   is used and recycled as needed.
#' @param geom Plot geometry: "auto", "line", or "bar". If "auto", "bar" is used
#'   for what = "LEdecomp" and "line" for what = "sens".
#' @param lwd Line width for cause lines (default 1.2).
#' @param xlab,ylab,main Axis labels and main title. If `ylab` is `NULL`,
#'   the default is `"Difference explained (years)"` for `what = "LEdecomp"`
#'   and "`Sensitivity d(e0)/d(mx)"` for `what = "sens"`.
#' @param legend Logical. Show legend (primarily relevant for `layout = "overlay"`).
#' @param legend_pos Legend position passed to `theme(legend.position = ...)`.
#'   Accepts `"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`,
#'   or numeric coordinates `c(x, y)` in \[0,1\] inside the panel.
#' @param abridged_midpoints Logical. If `TRUE` and ages are abridged
#'   `(0, 1, 5, 10, ...)`, plot against bin midpoints instead of lower bounds.
#' @param layout Plot layout for cause-of-death results: `"overlay"` (all causes
#'   in one panel) or `"facet"` (one panel per cause).
#' @param ncol Number of columns to use when `layout = "facet"`. If `NULL`,
#'   a default is chosen based on the number of causes.
#' @param ... Reserved for future use.
#'
#' @return Invisibly returns the ggplot object (after printing).
#' @importFrom stats ave
#' @importFrom utils head
#' @examples
#'
## Example 1: All-cause (use All-causes rows)
## US_data_CoD has: Period, Gender, Age, cause, cause_id, mxt
#' data("US_data_CoD", package = "LEdecomp")
#' allc <- subset(US_data_CoD, Period == 2010 & cause == "All-causes") |>
#'   as.data.frame()
#'
#' # Make Female vs Male all-cause schedules, Age 0:100
#' ac_w <- reshape(allc[, c("Gender","Age","mxt")],
#'                 timevar = "Gender", idvar = "Age", direction = "wide")
#' names(ac_w) <- sub("^mxt\\.", "", names(ac_w))
#' ac_w <- ac_w[order(ac_w$Age), ]
#'
#' dec_ac <- LEdecomp(
#'   mx1 = ac_w$Male,
#'   mx2 = ac_w$Female,
#'   age = 0:100,
#'   method = "sen_arriaga"
#' )
#'
#' # Simple single-line plot
#' \donttest{
#' plot(dec_ac, main = "All-cause Arriaga, 2010 Female vs Male")
#' }
#' ## End(Not run)
#' ## Example 2: Cause of death, one year, Female vs Male
#' cod <- subset(US_data_CoD, Period == 2010 & cause != "All-causes")
#' cod_w <- reshape(cod[, c("Gender","Age","cause","mxt")],
#'                  timevar = "Gender", idvar = c("cause","Age"),
#'                  direction = "wide")|>
#'   as.data.frame()
#' names(cod_w) <- sub("^mxt\\.", "", names(cod_w))
#' cod_w <- cod_w[order(cod_w$cause, cod_w$Age), ]
#'
#' dec_cod <- LEdecomp(
#'   mx1 = cod_w$Male,
#'   mx2 = cod_w$Female,
#'   age = 0:100,
#'   n_causes = length(unique(cod_w$cause)),
#'   cause_names = unique(cod$cause_id),
#'   method = "sen_arriaga"
#' )
#'
#' # Overlay of all causes
#' \donttest{
#' plot(dec_cod, layout = "overlay", main = "Arriaga CoD, 2010 Female vs Male", legend.pos = "top")
#'
#' # Facet by cause (3 columns)
#' plot(dec_cod, layout = "facet", ncol = 3, main = "Arriaga by cause (faceted)")
#' }
#'
#' ## Example 3: How to add an all-cause total line yourself (overlay)
#' \donttest{
#' p <- plot(dec_cod, layout = "overlay", main = "Overlay with manual Total")
#' y_mat <- if (is.matrix(dec_cod$LEdecomp)) dec_cod$LEdecomp else
#'   matrix(dec_cod$LEdecomp, nrow = length(dec_cod$age))
#' total <- rowSums(y_mat)
#' p + ggplot2::geom_line(
#'   data = data.frame(age = dec_cod$age, total = total),
#'   mapping = ggplot2::aes(x = .data$age, y = .data$total),
#'   inherit.aes = FALSE, color = "black", linewidth = 1.1)
#' }
#' @import ggplot2
#' @export
plot.LEdecomp <- function(x,
                          what = c("LEdecomp", "sens"),
                          geom = c("auto", "line", "bar"),
                          col = NULL,
                          lwd = 1.2,
                          xlab = "Age",
                          ylab = NULL,
                          main = NULL,
                          legend = TRUE,
                          legend_pos = "right",
                          abridged_midpoints = FALSE,
                          layout = c("overlay", "facet"),
                          ncol = NULL,
                          ...) {
  if (is.null(x) || !"LEdecomp" %in% class(x)) {
    stop("Object is not of class 'LEdecomp'.")
  }
  what   <- match.arg(what)
  layout <- match.arg(layout)
  geom   <- match.arg(geom)

  if (geom == "auto") {
    geom <- if (what == "LEdecomp") "bar" else "line"
  }

  age <- x$age
  if (is.null(age)) stop("No 'age' vector found in object.")
  nages <- length(age)

  y_raw <- if (what == "LEdecomp") x$LEdecomp else x$sens
  if (is.null(y_raw)) stop("Requested series '", what, "' is not available in object.")

  as_matrix <- function(y, nages) {
    if (is.matrix(y)) return(y)
    ny <- length(y)
    if (ny %% nages == 0L && ny > nages) {
      ncauses <- as.integer(ny / nages)
      matrix(y, nrow = nages, ncol = ncauses)
    } else {
      y
    }
  }
  y <- as_matrix(y_raw, nages)

  if (!is.null(x$nx)) {
    nx <- x$nx
  } else {
    if (length(age) > 1L) {
      dx <- diff(age)
      nx <- c(dx, dx[length(dx)])
    } else {
      nx <- 1
    }
  }

  if (isTRUE(abridged_midpoints) && geom != "bar") {
    age_plot <- age
    if (length(age) >= 3L && age[1] == 0 && age[2] == 1 && all(diff(age[-(1:2)]) %in% c(4, 5))) {
      dx <- diff(age)
      dx <- c(dx, dx[length(dx)])
      age_plot <- age + dx / 2
    }
  } else {
    age_plot <- age
  }

  if (is.null(ylab)) {
    ylab <- if (what == "LEdecomp") "Difference explained (years)" else "Sensitivity d(e0)/d(mx)"
  }
  if (is.null(main)) {
    main <- paste0("Method: ", x$method)
  }

  cause_names <- NULL
  if (!is.null(x$cause_names)) {
    cause_names <- if (is.factor(x$cause_names)) as.character(x$cause_names) else x$cause_names
  } else {
    if (is.matrix(y)) cause_names <- colnames(y)
    if (is.null(cause_names)) {
      mx1mat <- try(as.matrix(x$mx1), silent = TRUE)
      if (!inherits(mx1mat, "try-error")) {
        cause_names <- colnames(mx1mat)
      }
    }
  }
  if (is.matrix(y) && is.null(cause_names)) {
    cause_names <- paste0("cause_", seq_len(ncol(y)))
  }

  # colors
  if (is.matrix(y)) {
    ncauses <- ncol(y)
    if (is.null(col)) {
      base_cols <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                     "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
      col <- rep(base_cols, length.out = ncauses)
    } else if (length(col) < ncauses) {
      col <- rep(col, length.out = ncauses)
    }
    names(col) <- cause_names
  } else {
    if (is.null(col)) col <- "steelblue"
    names(col) <- "All-cause"
  }
  lw_vals <- rep(lwd, length(col)); names(lw_vals) <- names(col)

  ## ------------------------ line branch ------------------------
  if (geom == "line") {
    if (is.matrix(y)) {
      df_long <- data.frame(
        age_plot = rep(age_plot, times = ncol(y)),
        value    = as.vector(y),
        cause    = rep(cause_names, each = length(age_plot)),
        stringsAsFactors = FALSE
      )
    } else {
      df_long <- data.frame(
        age_plot = age_plot,
        value    = as.numeric(y),
        cause    = "All-cause",
        stringsAsFactors = FALSE
      )
    }
    df_long <- df_long[is.finite(df_long$value), , drop = FALSE]

    p <- ggplot(df_long,
                aes(x = .data$age_plot,
                    y = .data$value,
                    color = .data$cause,
                    linewidth = .data$cause)) +
      geom_line() +
      scale_color_manual(values = col, guide = "legend") +
      scale_linewidth_manual(values = lw_vals, guide = "legend") +
      labs(x = xlab, y = ylab, title = main) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank()) +
      geom_hline(yintercept = 0, color = "grey80")

    if (!isTRUE(legend)) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + theme(legend.position = legend_pos)
    }

    if (layout == "facet" && is.matrix(y)) {
      if (is.null(ncol)) {
        n_panels <- length(unique(df_long$cause))
        ncol <- if (n_panels <= 3) n_panels else 3L
      }
      p <- p + facet_wrap(~ .data$cause,
                          ncol = ncol,
                          scales = "fixed",
                          labeller = ggplot2::label_value)
    }

    print(p)
    return(invisible(p))
  }

  ## ------------------------ bar branch ------------------------
  if (is.matrix(y)) {
    # build long in cause order
    df_long <- data.frame(
      age_lower = rep(age, times = ncol(y)),
      nx        = rep(nx,  times = ncol(y)),
      value_raw = as.vector(y),
      cause     = factor(rep(cause_names, each = length(age)), levels = cause_names),
      stringsAsFactors = FALSE
    )
    # height per unit age
    df_long$height <- df_long$value_raw / df_long$nx

    if (layout == "overlay") {
      # stack within age
      df_long <- df_long[order(df_long$age_lower, df_long$cause), ]
      df_long$ymin <- ave(df_long$height, df_long$age_lower,
                          FUN = function(z) c(0, head(cumsum(z), -1)))
      df_long$ymax <- df_long$ymin + df_long$height
    } else {
      # facet: each cause starts at zero
      df_long$ymin <- 0
      df_long$ymax <- df_long$height
    }

    p <- ggplot(df_long,
                aes(fill = .data$cause)) +
      geom_rect(aes(xmin = .data$age_lower,
                    xmax = .data$age_lower + .data$nx,
                    ymin = .data$ymin,
                    ymax = .data$ymax),
                color = NA) +
      scale_fill_manual(values = col, guide = "legend") +
      labs(x = xlab, y = ylab, title = main) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank()) +
      geom_hline(yintercept = 0, color = "grey80")

    if (!isTRUE(legend)) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + theme(legend.position = legend_pos)
    }

    if (layout == "facet") {
      if (is.null(ncol)) {
        n_panels <- length(unique(df_long$cause))
        ncol <- if (n_panels <= 3) n_panels else 3L
      }
      p <- p +
        facet_wrap(~ .data$cause,
                   ncol = ncol,
                   scales = "fixed",
                   labeller = ggplot2::label_value)
    }

    print(p)
    return(invisible(p))

  } else {
    # single series
    df_long <- data.frame(
      age_lower = age,
      nx        = nx,
      height    = as.numeric(y) / nx,
      stringsAsFactors = FALSE
    )

    p <- ggplot(df_long, aes()) +
      geom_rect(aes(xmin = .data$age_lower,
                    xmax = .data$age_lower + .data$nx,
                    ymin = 0,
                    ymax = .data$height),
                fill = col %||% "steelblue", color = NA) +
      labs(x = xlab, y = ylab, title = main) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank()) +
      geom_hline(yintercept = 0, color = "grey80")

    if (!isTRUE(legend)) {
      p <- p + theme(legend.position = "none")
    }

    print(p)
    return(invisible(p))
  }
}
