utils::globalVariables(c("name", "station", "time", "prediction",
                         "unit", "r2", "type", "value", "max_val",
                         "value.scaled", "n", "stations", "metric_name",
                         "rr.label", "breaks", "prop"))

#' Generate Scatter Plot
#'
#' @description Creates scatterplots of associated events for all
#'     metrics and pairs of stations including the fitted regression
#'     line.
#'
#' @param AE A data frame. First list element of output of
#'     \code{\link[=estimate_settings_AE]{estimate_settings_AE()}}.
#' @param axis A logical. If \code{TRUE} (default), the aspect ratio
#'     of the plots are fixed and the y-axis and x-axis are set equal.
#' @param accuracy A number to round to. Use (e.g.) 0.01 to show 2 decimal places
#'     of precision. If NULL, no rounding is applied.
#'
#' @return Returns a `\code{gtable}' object.
#'
#' @keywords internal
#' @noRd
#' @examples
#' # file paths
#' Sx <- system.file("testdata", "Events", "100000_2_2014-01-01_2014-02-28.csv",
#'   package = "hydroroute")
#' Sy <- system.file("testdata", "Events", "200000_2_2014-01-01_2014-02-28.csv",
#'   package = "hydroroute")
#' relation <- system.file("testdata", "relation.csv", package = "hydroroute")
#'
#' # read data
#' Sx <- utils::read.csv(Sx)
#' Sy <- utils::read.csv(Sy)
#'
#' relation <- utils::read.csv(relation)
#' relation <- relation[1:2, ]
#'
#'
#' results <- estimate_AE(Sx, Sy, relation)
#' real_AE <- results$real_AE
#' plt <- hydroroute:::plot_scatter(real_AE)
#' plot.new()
#' grid::grid.draw(plt)
plot_scatter <- function(AE, axis = TRUE, accuracy = NULL) {
  AE <- AE[, !(names(AE) %in% c("diff_metric"))]

  data_long <- reshape(AE)

  metrics <- data.frame(new = c("AMP", "MAFR", "MEFR", "DUR", "RATIO"),
                        unit = c("m\u00b3/s", "(m\u00b3/s)/ts", "(m\u00b3/s)/ts",
                                 "ts", "-"))

  data_long <- data_long |>
    dplyr::mutate(stations = paste0(station.x, "-", station.y),
                  unit = factor(toupper(metric), metrics$new, metrics$unit),
                  metric_name = factor(toupper(metric), metrics$new,
                                       with(metrics,
                                            paste0(new, " (", unit, ")"))))

  plots_scatter <- lapply(levels(data_long$metric_name), function(m) {
    data_sub <- data_long |>
      dplyr::mutate(stations = factor(stations)) |>
      dplyr::filter(metric_name == m) |>
      dplyr::mutate(metric_name = factor(m)) |>
      dplyr::filter(!is.na(x), !is.na(y))

    gg1 <- ggplot2::ggplot(data_sub, ggplot2::aes(x, y)) +
      ggplot2::facet_wrap(metric_name ~ stations, scales = "free", nrow = 1,
                          drop = FALSE) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::xlab(data_sub$unit[1]) +
      ggplot2::ylab(data_sub$unit[1]) +
      ggplot2::geom_smooth(formula = y ~ x, method = "lm") +
      ggpmisc::stat_poly_eq(ggplot2::aes(label = ggplot2::after_stat(rr.label)),
                            formula = y ~ x, colour = "blue",
                            label.x.npc = "right", label.y.npc = "bottom")

    if (axis) {
      limits <- unlist(lapply(seq_along(levels(data_sub$stations)),
                              function(s) {
                                lims <- c(ggplot2::layer_scales(gg1, 1, s)[["x"]]$range$range,
                                          ggplot2::layer_scales(gg1, 1, s)[["y"]]$range$range)
                                suppressWarnings(range(lims))
      }))
      dummy <- data.frame(x = limits, y = limits,
                          metric_name = m,
                          stations = factor(rep(levels(data_sub$stations),
                                                each = 2),
                                            levels(data_sub$stations)),
                          stringsAsFactors = FALSE)
      gg1 <- gg1 + ggplot2::geom_blank(data = dummy) +
        ggplot2::scale_x_continuous(guide = ggplot2::guide_axis(check.overlap = TRUE),
                                    if (!is.null(accuracy)) labels = scales::label_number(accuracy = accuracy)) +
        ggplot2::scale_y_continuous(if (!is.null(accuracy)) labels = scales::label_number(accuracy = accuracy)) +
        ggplot2::theme(aspect.ratio = 1)
    }})

    plots_scatter <- do.call(gridExtra::arrangeGrob,
                             c(plots_scatter,
                               list(nrow = length(plots_scatter))))
    plots_scatter
}

#' Generate Threshold Plot
#'
#' @description Called within \code{\link[=estimate_settings_AE]{estimate_settings_AE()}}.
#'
#' @param freq_dist Table containing relative frequencies.
#' @param pars Estimated parameters returned by
#'     \code{\link[=get_parabola]{get_parabola()}}.
#' @param bounds Suitable cut points based on the estimated parameters
#'     or a strict criterion \code{\link[=get_bounds]{get_bounds()}}.
#' @param Sx Character string containing the station (default:
#'     \code{"S1"}).
#' @param Sy Character string containing the station (default:
#'     \code{"S2"}).
#'
#' @return Returns a `\code{ggplot}' object that contains the
#'     histogram plot with the fitted parabola and added cut points.
#'
#' @keywords internal
#' @noRd
plot_threshold <- function(freq_dist, pars, bounds, Sx = "S1", Sy = "S2") {
  prop_dist <- freq_dist / sum(freq_dist)

  plot_threshold <- ggplot2::ggplot(data.frame(breaks = seq(-0.95, 0.95, by = 0.1),
                                               prop = unname(as.vector(prop_dist))),
                                    ggplot2::aes(x = breaks,
                                                 y = prop)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste0(Sx, " - ", Sy))

  if (!missing(bounds)) {
      plot_threshold <- plot_threshold +
          ggplot2::geom_vline(xintercept = bounds, linetype = 2)
  }

  if (!missing(pars) && length(pars) != 0) {
    x_axis <- seq(-0.95, 0.95, length.out = 100)
    index <- seq(1, length(freq_dist), length.out = 100)
    y_axis <- pars$coefs[1] * (index - pars$coefs[3])^2 + pars$coefs[2]
    df <- data.frame(x = x_axis, y = y_axis) |>
      dplyr::filter(y > 0)

    plot_threshold <- plot_threshold +
      ggplot2::geom_line(ggplot2::aes(x, y), data = df,
                         linetype = 2)
  }

  plot_threshold
}

#' Generate Predictions Plot
#'
#' @description Plots the predicted values of the different metrics in
#'     dependence of the distance to the source (\code{fkm}) as well
#'     as the R-squared values of the fitted regressions and the
#'     number of associated events identified.
#'
#' @param predictions Data frame
#'
#' @return Returns a `\code{gtable}' object.
#'
#' @keywords internal
#' @noRd
plot_predict <- function(predictions) {
  Ns <- stats::na.omit(unique(predictions[, c("station", "n", "fkm")]))

  P <- predictions |>
    dplyr::select(fkm, metric, name, prediction, unit) |>
    dplyr::mutate(type = "P") |>
    dplyr::rename(value = prediction)

  R2 <- predictions |>
    dplyr::select(fkm, metric, r2, unit) |>
    dplyr::mutate(name = NA, type = "R2") |>
    dplyr::rename(value = r2) |> unique() |>
    dplyr::filter(!is.na(value))

  x <- rbind(P, R2[, colnames(P)])

  Max <- x |>
    dplyr::group_by(type, metric) |>
    dplyr::summarise(max_val = max(value, na.rm = TRUE), .groups = "drop")

  Max <- Max |>
    dplyr::filter(type == "P") |>
    dplyr::ungroup() |>
    dplyr::select(max_val, metric)

  x <- x |>
    dplyr::left_join(Max, by = "metric") |>
    dplyr::mutate(value.scaled = ifelse(type == "P", value, value * max_val))

  plots <- c(
    lapply(levels(x$metric), function(m) {
      x_metric <- dplyr::filter(x, metric == m)
      max_val <- max(x_metric$max_val)
      ggplot2::ggplot(x_metric,
                      ggplot2::aes(fkm, value.scaled, linetype = factor(type))) +
        ggplot2::facet_wrap(~ metric) +
        ggplot2::geom_line(ggplot2::aes(group = interaction(factor(name),
                                                            factor(type)))) +
        ggplot2::geom_point(shape = 20) +
        ggplot2::geom_blank(ggplot2::aes(y = 0)) +
        ggplot2::xlab("distance (km)") +
        ggplot2::xlim(range(predictions$fkm, na.rm = TRUE)) +
        ggplot2::scale_y_continuous(x_metric$unit[1],
                                    limits = c(0, max_val),
                                    sec.axis = ggplot2::sec_axis(~ . / max_val,
                                                                 name = "R-squared",
                                                                 breaks = seq(0, 1, length.out = 5))) +
        ggplot2::theme_bw() + ggplot2::guides(linetype = "none")}),
    list(ggplot2::ggplot(Ns,
                         ggplot2::aes(fkm, n)) +
           ggplot2::facet_wrap(~ "Number of events") +
           ggplot2::geom_line() +
           ggplot2::geom_point(shape = 20) +
           ggplot2::xlab("distance (km)") +
           ggplot2::xlim(range(predictions$fkm, na.rm = TRUE)) +
           ggplot2::scale_y_continuous("Count",
                                       limits = c(0, max(Ns$n, na.rm = TRUE)),
                                       sec.axis = ggplot2::sec_axis(~ .)) +
           ggplot2::theme_bw())
    )

  gridExtra::arrangeGrob(grobs = plots, nrow = 1, ncol = length(plots))
}
