utils::globalVariables(c("metric", "station.x", "time.x", "station.y", "time.y",
                         "diff_metric", "y", "x", "id.x", "id.y", "fkm"))

#' Potential Associated Events
#' @description Captures several events from the downstream hydrograph
#'     \eqn{S_{x}}{S_x}. It captures either exact time matches according to the
#'     mean translation time, or matches within a time interval. If several matches
#'     at the downstream hydrograph \eqn{S_{y}}{S_y} are detected that match an
#'     event at the upstream hydrograph, the event at \eqn{S_{y}}{S_y} with the
#'     smallest time difference to the event at \eqn{S_{x}}{S_x} is chosen.
#'     Further, an event detected at \eqn{S_{y}}{S_y} with amplitude ("AMP") most
#'     similar to the upstream hydrograph \eqn{S_{x}}{S_x} is designated as
#'     potential associated event (AE). Therefore, the metric deviation is computed
#'     by \eqn{\frac{(S_y - S_x)}{S_x}} and events where the amplitude at
#'     \eqn{S_{y}}{S_y} is at most two times the amplitude at \eqn{S_{x}}{S_x}
#'     are detected as potential AEs.
#' @param Sx Data frame that consists of flow fluctuation events and computed
#'     metrics (see \code{\link[hydropeak:get_events]{hydropeak::get_events()}})
#'     of an upstream hydrograph \eqn{S_{x}}{S_x}.
#' @param Sy Data frame that consists of flow fluctuation events and computed
#'     metrics (see \code{\link[hydropeak:get_events]{hydropeak::get_events()}})
#'     of a downstream hydrograph \eqn{S_{y}}{S_y}.
#' @param relation Data frame that contains the relation between upstream and
#'     downstream hydrograph. Must only contain two rows (one for each hydrograph)
#'     in order of their location in downstream direction.
#'     See the appended example data \code{relation.csv} or vignette for details on
#'     the structure. See \code{\link[=get_lag]{get_lag()}} for further information
#'     about the relation and the lag between the hydrographs.
#' @param timeLag Numeric vector specifying factors to alter the interval to capture
#'     events from the downstream hydrograph. By default it is
#'     \code{timeLag = c(1, 1, 1)}, this refers to matches within a time slot
#'     \eqn{\pm} the mean translation time from \code{relation}. For exact time
#'     matches, \code{timeLag = c(0, 1, 0)} must be specified.
#' @param metricLag Numeric vector specifying factor to alter the interval to capture
#'     events from the downstream hydrograph. By default it is
#'     \code{metricLag = c(1, 1)}, such that events are filtered where the amplitude
#'     at \eqn{S_{y}}{S_y} is at least 0, i.e., amplitude at
#'     \eqn{S_{x} - 1 \cdot} amplitude at \eqn{S_{x}}{S_x}, and at most two times the
#'     amplitude at \eqn{S_{x}}{S_x}, i.e., \eqn{S_{x} + 1 \cdot} amplitude at \eqn{S_{x}}{S_x}.
#'     For exact matches, \code{metricLag = c(0, 0)} must be specified.
#' @param unique Character string specifying if the potential AEs which
#'     meet the \code{timeLag} and \code{metricLag} condition should
#'     be filtered to contain only unique events using \code{"time"},
#'     i.e., by selecting those where the time difference is smallest
#'     compared to the specified factor of the mean translation time, or using
#'     \code{"metric"}, i.e., by selecting those where the relative
#'     difference in amplitude is smallest (default: \code{"time"}).
#' @param TimeFormat Character string giving the date-time format of the
#'     date-time column in the input data frame  (default: "\%Y-\%m-\%d \%H:\%M").
#' @param tz Character string specifying the time zone to be used for the
#'     conversion (default: "Etc/GMT-1").
#'
#' @return Data frame where each row represents an event at the upstream
#'     hydrograph \eqn{S_{x}}{S_x} combined with an event at the downstream
#'     hydrograph \eqn{S_{y}}{S_y} for each metric (long format). Columns \code{x}
#'     and \code{y} contain the values of the metrics (column \code{Metric}) for
#'     upstream (\code{x}) and downstream (\code{y}) hydrograph. \code{diff_metric}
#'     is the computed metric deviation. If no potential AEs are found, \code{NULL}
#'     is returned.
#' @keywords internal
#' @noRd
potential_AE <- function(Sx, Sy, relation, timeLag = c(1, 1, 1), metricLag = c(1, 1),
                         unique = c("time", "metric"),
                         TimeFormat = "%Y-%m-%d %H:%M", tz = "Etc/GMT-1") {

  unique <- match.arg(unique)

  colnames(relation) <- tolower(colnames(relation))
  colnames(Sx) <- tolower(colnames(Sx))
  colnames(Sy) <- tolower(colnames(Sy))

  # argument Sx, Sy, relation checked in estimate_AE()
  # arguments TimeFormat and tz are checked in shift_time()

  if (is.null(Sx$station)) {
    Sx$station <- rep(relation$station[1], nrow(Sx))
  }

  if (is.null(Sy$station)) {
    Sy$station <- rep(relation$station[2], nrow(Sy))
  }

  # lag_opt
  lag <- ifelse(is.na(relation$lag), "00:00", as.character(relation$lag))
  lag_opt <- diff(hhmm_to_min(lag)) * timeLag[2]

  # get all events with time matches
  events <- merge_time(Sx, Sy, relation, timeLag, TimeFormat, tz)

  potential_AE <- NULL
  if (!is.null(events)) {
    # reshape to long format
    events_long <- reshape(events)

    # filter metric, relative difference in metric
    potential_AE <- events_long |>
      dplyr::filter(metric == "amp") |>
      dplyr::mutate(diff_metric = (y - x) / x)

    # potential AEs: event detected at Sy with metric within metricLag range
    potential_AE <- potential_AE |>
      dplyr::filter(y >= (x - metricLag[1] * x), y <= (x + metricLag[2] * x))

    if (nrow(potential_AE) > 0) {
      unique_AE <- NULL
      # for each event at Sx, filter event at Sy with smallest time / AMP difference
      while (nrow(potential_AE) > 0) {
        if (unique == "time") {
          best <- potential_AE |>
            dplyr::group_by(time.x) |>
            dplyr::slice(which.min(abs(time.x - time.y + lag_opt))) |>
            dplyr::group_by(time.y) |>
            dplyr::slice(which.min(abs(time.x - time.y + lag_opt))) |>
            dplyr::ungroup()
        } else {
          best <- potential_AE |>
            dplyr::group_by(time.x) |>
            dplyr::slice(which.min(abs(y - x))) |>
            dplyr::group_by(time.y) |>
            dplyr::slice(which.min(abs(y - x))) |>
            dplyr::ungroup()
        }

        potential_AE <- potential_AE |>
          dplyr::anti_join(best |> dplyr::select(time.x),
                           by = "time.x") |>
          dplyr::anti_join(best |> dplyr::select(time.y),
                           by = "time.y")
        unique_AE <- rbind(unique_AE, best)
      }

      # all potential AEs incl. all metrics
      potential_AE <- events_long |>
        dplyr::inner_join(unique_AE |> dplyr::select(station.x, time.x, station.y,
                                                     time.y, diff_metric),
                          by = c("station.x", "time.x", "station.y", "time.y"))

      # update diff_metric for all metrics
      potential_AE <-  potential_AE |>
        dplyr::mutate(diff_metric = (y - x) / x)
    } else {
      potential_AE <- NULL
    }
  }

  potential_AE
}

#' Merge Events
#'
#' @description Given two event data frames of neighboring stations
#'     \eqn{S_{x}}{S_x} and \eqn{S_{y}}{S_y} that consist of flow
#'     fluctuation events and computed metrics (see
#'     \code{\link[hydropeak:get_events]{hydropeak::get_events()}}),
#'     the translation time indicated by the relation file as well as
#'     \code{timeLag} between these two stations is subtracted from
#'     \eqn{S_{y}}{S_y} and events are merged where matches according
#'     to differences allowed to \code{timeLag} can be found.
#'
#' @param Sx Data frame that consists of flow fluctuation events and computed
#'     metrics (see \code{\link[hydropeak:get_events]{hydropeak::get_events()}})
#'     of an upstream hydrograph \eqn{S_{x}}{S_x}.
#' @param Sy Data frame that consists of flow fluctuation events and computed
#'     metrics (see \code{\link[hydropeak:get_events]{hydropeak::get_events()}})
#'     of a downstream hydrograph \eqn{S_{y}}{S_y}.
#' @param relation Data frame that contains the relation between upstream and
#'     downstream hydrograph. Must only contain two rows (one for each hydrograph)
#'     in order of their location in downstream direction.
#'     See the appended example data \code{relation.csv} or vignette for details on
#'     the structure. See \code{\link[=get_lag]{get_lag()}} for further information
#'     about the relation and the lag between the hydrographs.
#' @param timeLag Numeric vector specifying factors to alter the interval to capture
#'     events from the downstream hydrograph. By default it is
#'     \code{timeLag = c(1, 1, 1)}, this refers to matches within a time slot
#'     \eqn{\pm} the mean translation time from \code{relation}. For exact time
#'     matches, \code{timeLag = c(0, 1, 0)} must be specified.
#' @param TimeFormat Character string giving the date-time format of the
#'     date-time column in the input data frame  (default: "\%Y-\%m-\%d \%H:\%M").
#' @param tz Character string specifying the time zone to be used for the
#'     conversion (default: "Etc/GMT-1").
#'
#' @return Data frame that has a matched event at \eqn{S_{x}}{S_x} and
#'     \eqn{S_{y}}{S_y} in each row. If no matches are detected,
#'     \code{NULL} is returned.
#' @importFrom lubridate %within%
#' @export
#' @examples
#' Sx <- system.file("testdata", "Events", "100000_2_2014-01-01_2014-02-28.csv",
#'                   package = "hydroroute")
#' Sy <- system.file("testdata", "Events", "200000_2_2014-01-01_2014-02-28.csv",
#'                   package = "hydroroute")
#' relation <- system.file("testdata", "relation.csv", package = "hydroroute")
#' # read data
#' Sx <- utils::read.csv(Sx)
#' Sy <- utils::read.csv(Sy)
#' relation <- utils::read.csv(relation)
#' relation <- relation[1:2, ]
#'
#' # exact matches
#' merged <- merge_time(Sx, Sy, relation, timeLag = c(0, 1, 0))
#' head(merged)
#'
#' # matches within +/- mean translation time
#' merged <- merge_time(Sx, Sy, relation)
#' head(merged)
merge_time <- function(Sx, Sy, relation, timeLag = c(1, 1, 1),
                       TimeFormat = "%Y-%m-%d %H:%M", tz = "Etc/GMT-1") {
  # argument Sx, Sy, relation checked in estimate_AE()
  colnames(relation) <- tolower(colnames(relation))
  colnames(Sx) <- tolower(colnames(Sx))
  colnames(Sy) <- tolower(colnames(Sy))

  if (is.null(Sx$station)) {
    Sx$station <- rep(relation$station[1], nrow(Sx))
  }

  if (is.null(Sy$station)) {
    Sy$station <- rep(relation$station[2], nrow(Sy))
  }

  relation$lag <- ifelse(is.na(relation$lag), "00:00",
                         as.character(relation$lag))
  relation$lag <- hhmm_to_min(relation$lag)

  # RATIO does not make sense for S1 if not type 'Gauge'
  if ("S1" %in% relation$station) {
    Sx$ratio <- ifelse(Sx$station == "S1" &
                       subset(relation, station == "S1")$type != "Gauge",
                       NA_real_, Sx$ratio)
  }

  # make sure Sx and Sy are sorted by Time
  Sx <- dplyr::arrange(Sx, time)
  Sy <- dplyr::arrange(Sy, time)

  # optimal lag: mean translation time between two stations (get_lag())
  lag_opt <- diff(relation$lag) * timeLag[2]

  # subtract optimal lag from downstream hydrograph Sy
  Sy$time <- shift_time(datetime = Sy$time, shift = lag_opt,
                        TimeFormat = TimeFormat, tz = tz)

  Sx$time <- as.POSIXct(Sx$time, format = TimeFormat, tz = tz)

  # indices with matches
  index_list <- vector(mode = "list", length = nrow(Sx))

  # time matches
  for (i in seq_len(length(index_list))) {
    int <- lubridate::interval(Sx[i, ]$time - (lag_opt * timeLag[1]), Sx[i, ]$time + (lag_opt * timeLag[3]))
    index_list[[i]] <- which(Sy$time %within% int)
  }

  index_df <- data.frame(x = rep(seq(length(index_list)), lengths(index_list)),
                         y = unlist(index_list))

  events <- NULL
  if (nrow(index_df) > 0) {
    index_split <- index_df |>
      dplyr::group_by(x)

    index_group_split <- dplyr::group_split(index_split)

    index_df_match <- do.call(rbind, index_group_split)

    Sx_new <- Sx[index_df_match$x, ]
    Sy_new <- Sy[index_df_match$y, ]

    # helper variable for merging
    Sx_new$help <- seq.int(1, nrow(Sx_new))
    Sy_new$help <- seq.int(1, nrow(Sy_new))

    events <- merge(Sx_new, Sy_new, by = "help")

    # drop helper variable
    events <- events[ , !(names(events) %in% "help")]

    # shift back
    events$time.y <- shift_time(datetime = events$time.y, shift = -lag_opt,
                                TimeFormat = TimeFormat, tz = tz)
  }
  events
}


#' Shift Time
#' @description Given a vector of date-time values, the date-time is shifted by
#'     the amount of time provided by argument \code{shift}: \code{datetime - shift}.
#'     It is used to shift the date-time values at the downstream hydrograph
#'     \eqn{S_{y}}{S_y} by the estimated translation time between \eqn{S_{x}}{S_x}
#'     and \eqn{S_{y}}{S_y}.
#' @param datetime Character or `\code{POSIXct}' vector containing date-time values to
#'     be shifted by \code{shift}.
#' @param shift `\code{difftime}' or integer value (minutes) that is subtracted from
#'     date-time values in \code{datetime}.
#' @param TimeFormat Character string giving the date-time format of the
#'     date-time column in the input data frame  (default: "\%Y-\%m-\%d \%H:\%M").
#' @param tz Character string specifying the time zone to be used for the
#'     conversion (default: "Etc/GMT-1").
#'
#' @return Returns the shifted date-time vector of class `\code{POSIXct}'.
#' @keywords internal
#' @noRd
shift_time <- function(datetime, shift, TimeFormat = "%Y-%m-%d %H:%M",
                       tz = "Etc/GMT-1") {
  stopifnot(length(shift) == 1)
  stopifnot(is.character(TimeFormat), length(TimeFormat) == 1)
  stopifnot(is.character(tz), length(tz) == 1)

  if (!lubridate::is.difftime(shift)) { # shift provided in minutes
    shift <- shift * 60
  }

  datetime <- as.POSIXct(datetime, format = TimeFormat, tz = tz)

  # shift date-time
  datetime - shift
}

#' Write Objects to csv
#' @description Writes the provided data frame or matrix to a csv file
#'     with path and file name as specified.
#' @param x Matrix or data frame containing the estimated settings or real AEs
#'     from \code{\link[=estimate_AE]{estimate_AE()}}.
#' @param outdir Character string naming a directory where the generated plot
#'     should be saved to.
#' @param metric Character string containing the related metric. It is part of
#'     the file name.
#' @param event Character string containing the event type. It is part of the
#'     file name.
#' @param station Character vector containing the related stations
#'     (e.g. \code{c("S1", "S2")}). It is part of the file name.
#'     It can be \code{NULL}.
#' @param postifx Character string containing one of
#'     \code{c("settings", "events", "models", "real_AE")}. It is part of the file name.
#' @return No return value, called for side effects.
#'
#' @keywords internal
#' @noRd
writeCSV <- function(x, outdir = tempdir(), metric, event, station = NULL,
                     postfix = c("settings", "events", "models", "real_AE")) {
  stopifnot(is.data.frame(x) | is.matrix(x))
  stopifnot(is.character(outdir) & length(outdir) == 1)

  postfix <- match.arg(postfix)
  stopifnot(!missing(event))
  prefix <- "Q_event_"
  event <- event
  stations <- c(station[1], station[2])

  if (!missing(metric)) {
    stopifnot(is.character(metric) & length(metric) == 1)
    metric <- toupper(metric)
    name <- paste(metric, "LAG", sep = "-")
  } else {
    name <- "LAG"
  }

  if (is.null(station)) {
    path <- file.path(outdir, paste0(prefix, event, "_", name, "_",
                                     postfix, ".csv"))
  } else {
    stopifnot(is.character(station))
    path <- file.path(outdir, paste0(prefix, event, "_", name, "_",
                                     paste(stations, collapse = "_"), "_",
                                     postfix, ".csv"))
  }

  utils::write.table(x,
                     file = path,
                     sep = ",", dec = ".", row.names = FALSE)
}


#' Get cut points
#' @description Rescale the cut points from index scale to correspond
#'     to the interval of relative differences or use a strict
#'     criterion bound.
#' @param pars Estimated parameters returned by
#'     \code{\link[=get_parabola]{get_parabola()}}.
#'
#' @return Numeric vector of length two containing the cut points.
#'
#' @keywords internal
#' @noRd
get_bounds <- function(pars) {
  if (length(pars) == 0) {
    # strict criterion
    bounds <- c(-0.1, 0.1)
  } else {
    trafo <- c(a = 2 * 0.95 / 19,
               b = -0.95 - 2 * 0.95 / 19)
    bounds <- trafo["a"] * pars$xstar + trafo["b"]

    # strict criterion
    if (any(bounds < -1) || any(bounds > 1)) {
      bounds <- c(-0.1, 0.1)
    }
  }
  bounds
}


#' Estimate Associated Events
#' @description For two neighboring stations, potential associated
#'     events (AEs) are determined according to the time lag and
#'     metric (amplitude) difference allowed. For all potential AEs,
#'     parabolas are fitted to the histogram obtained for the relative
#'     difference in amplitude binned into intervals from -1 to 1 of
#'     width 0.1 by fixing the vertex at the inner maximum of the
#'     histogram and the width is determined by minimizing the average
#'     squared distances between the parabola and the histogram data
#'     along arbitrary symmetric ranges from the inner maximum.  Based
#'     on the fitted parabola, cut points with the x-axis are
#'     determined such that only those potential AEs are retained
#'     where the relative difference is within these cut points. If
#'     this automatic scheme does not succeed to determine suitable
#'     cut points, e.g., because the estimated cut points are outside
#'     -1 and 1, then a strict criterion for the relative difference
#'     in amplitude is imposed to identify AEs considering only
#'     deviations of at most 10\%.
#' @param Sx Data frame that consists of flow fluctuation events and
#'     computed metrics (see
#'     \code{\link[hydropeak:get_events]{hydropeak::get_events()}}) of
#'     an upstream hydrograph \eqn{S_{x}}{S_x}.
#' @param Sy Data frame that consists of flow fluctuation events and
#'     computed metrics (see
#'     \code{\link[hydropeak:get_events]{hydropeak::get_events()}}) of
#'     a downstream hydrograph \eqn{S_{y}}{S_y}.
#' @param relation Data frame that contains the relation between
#'     upstream and downstream hydrograph. Must only contain two rows
#'     (one for each hydrograph) in order of their location in
#'     downstream direction.  See the appended example data
#'     \code{relation.csv} or the vignette for details on the
#'     structure. See \code{\link[=get_lag]{get_lag()}} for further
#'     information about the relation and the lag between the
#'     hydrographs.
#' @param timeLag Numeric vector specifying factors to alter the
#'     interval to capture events from the downstream hydrograph. By
#'     default it is \code{timeLag = c(1, 1, 1)}, this refers to
#'     matches within a time slot \eqn{\pm} the mean translation time
#'     from \code{relation}. For exact time matches, \code{timeLag =
#'     c(0, 1, 0)} must be specified.
#' @param metricLag Numeric vector specifying factors to alter the
#'     interval of relative metric deviations to capture events from
#'     the downstream hydrograph. By default. it is \code{metricLag =
#'     c(1, 1)}, such that events are filtered where the amplitude at
#'     \eqn{S_{y}}{S_y} is at least 0, i.e., amplitude at \eqn{S_{x} -
#'     1 \cdot} amplitude at \eqn{S_{x}}{S_x}, and at most two times
#'     the amplitude at \eqn{S_{x}}{S_x}, i.e., \eqn{S_{x} + 1 \cdot}
#'     amplitude at \eqn{S_{x}}{S_x}.  For exact matches,
#'     \code{metricLag = c(0, 0)} must be specified.
#' @param unique Character string specifying if the potential AEs which
#'     meet the \code{timeLag} and \code{metricLag} condition should
#'     be filtered to contain only unique events using \code{"time"},
#'     i.e., by selecting those where the time difference is smallest
#'     compared to the specified factor of the mean translation time, or using
#'     \code{"metric"}, i.e., by selecting those where the relative
#'     difference in amplitude is smallest (default: \code{"time"}).
#' @param TimeFormat Character string giving the date-time format of
#'     the date-time column in the input data frame (default:
#'     "\%Y-\%m-\%d \%H:\%M").
#' @param tz Character string specifying the time zone to be used for
#'     the conversion (default: "Etc/GMT-1").
#' @param settings Data.frame with 3 rows and columns
#'     \code{station.x}, \code{station.y}, \code{bound}, \code{lag},
#'     \code{metric}. \code{lag} needs to correspond to the unique
#'     value specified in argument \code{timeLag} and \code{bound}
#'     needs to contain \code{"lower"}, \code{"inner"},
#'     \code{"upper"}.
#'
#' @return A nested list containing the estimated settings, the
#'     histogram obtained for the relative difference data with
#'     estimated cut points, and the obtained \dQuote{real} AEs.
#' @export
#'
#' @examples
#' # file paths
#' Sx <- system.file("testdata", "Events", "100000_2_2014-01-01_2014-02-28.csv",
#'                   package = "hydroroute")
#' Sy <- system.file("testdata", "Events", "200000_2_2014-01-01_2014-02-28.csv",
#'                   package = "hydroroute")
#' relation <- system.file("testdata", "relation.csv", package = "hydroroute")
#'
#' # read data
#' Sx <- utils::read.csv(Sx)
#' Sy <- utils::read.csv(Sy)
#' relation <- utils::read.csv(relation)
#' relation <- relation[1:2, ]
#'
#' # estimate AE, exact time matches
#' results <- estimate_AE(Sx, Sy, relation, timeLag = c(0, 1, 0))
#' results$settings
#' results$plot_threshold
#' results$real_AE
estimate_AE <- function(Sx, Sy, relation, timeLag = c(1, 1, 1), metricLag = c(1, 1),
                        unique = c("time", "metric"),
                        TimeFormat = "%Y-%m-%d %H:%M", tz = "Etc/GMT-1",
                        settings = NULL) {
  unique <- match.arg(unique)

  colnames(relation) <- tolower(colnames(relation))
  colnames(Sx) <- tolower(colnames(Sx))
  colnames(Sy) <- tolower(colnames(Sy))
  stopifnot(is.data.frame(Sx))
  stopifnot(is.data.frame(Sy))

  stopifnot(is.data.frame(relation))
  stopifnot(nrow(relation) == 2)

  stopifnot(length(unique(Sx$event_type)) == 1)
  stopifnot(length(unique(Sy$event_type)) == 1)
  stopifnot(all.equal(unique(Sx$event_type), unique(Sy$event_type)))

  stopifnot(length(unique(Sx$id)) == 1)
  stopifnot(length(unique(Sy$id)) == 1)

  # only one event per station and time
  stopifnot(length(unique(Sx$time)) == nrow(Sx))
  stopifnot(length(unique(Sy$time)) == nrow(Sy))

  stopifnot(unique(Sx$id) %in% relation$id)
  stopifnot(unique(Sy$id) %in% relation$id)

  # events, exact time matches
  events <- merge_time(Sx = Sx, Sy = Sy, relation = relation, timeLag = timeLag,
                       TimeFormat = TimeFormat, tz = tz)

  # potential AE, AMP within fixed interval
  potential_AE <- potential_AE(Sx = Sx, Sy = Sy, relation = relation,
                               timeLag = timeLag, metricLag = metricLag,
                               unique = unique, TimeFormat = TimeFormat, tz = tz)

  # lag in settings and timeLag must match
  if (!is.null(settings)) {
      stopifnot(length(unique(timeLag)) == 1)
      stopifnot(nrow(settings) == 3)
      stopifnot(all(settings$station.x == potential_AE$station.x[1]),
                all(settings$station.y == potential_AE$station.y[1]),
                all(settings$lag == unique(timeLag)),
                all(settings$bound %in% c("lower", "inner", "upper")),
                !anyDuplicated(settings$bound),
                !anyNA(settings$metric))
      settings <- settings |>
          dplyr::arrange(metric)
  }

  plot_threshold <- NULL
  real_AE <- NULL

  if (!is.null(potential_AE)) {
    # estimate cut points
    potential_AE_subset <- potential_AE |>
      dplyr::filter(metric == "amp")

    # relative difference data binned into intervals from -1 to 1 of width 0.1
    freq_dist <- table(cut(potential_AE_subset$diff_metric, seq(-1, 1, by = 0.1)))

    if (is.null(settings)) {
        prop_dist <- freq_dist / sum(freq_dist)

        pars <- get_parabola(prop_dist)

        # cut points (or strict criterion)
        bounds <- get_bounds(pars)

        # fitted cut points -> relative difference
        metric_est <- bounds + 1

        # estimated settings
        settings <- cbind.data.frame(station.x = rep(potential_AE$station.x[1], 3),
                                     station.y = rep(potential_AE$station.y[1], 3),
                                     bound = c("lower", "inner", "upper"),
                                     lag = timeLag,
                                     metric = c(metric_est[1], mean(metric_est), metric_est[2]))
    } else {
        pars <- NULL
        bounds <- settings$metric[c(1, 3)] - 1
    }

    # real AEs
    # filter events within estimated cut points (settings)
    events <- events |>
      dplyr::inner_join(potential_AE_subset |> dplyr::select(station.x, time.x, station.y,
                                                             time.y, diff_metric),
                        by = c("station.x", "time.x", "station.y", "time.y"))

    real_AE <- events |>
      dplyr::filter(dplyr::between(diff_metric + 1, settings$metric[1],
                                   settings$metric[3]))
    if (nrow(real_AE) == 0) {
        real_AE <- NULL
    }

    # threshold plot
    plot_threshold <- plot_threshold(freq_dist, pars, bounds,
                                     unique(potential_AE$station.x),
                                     unique(potential_AE$station.y))
  }

  return(list(settings = settings, plot_threshold = plot_threshold,
              real_AE = real_AE))
}

#' Trace Longitudinal Hydropeaking Waves Along a River Section
#' @description Estimates all settings based on the `relation' file of a river
#'     section. The function uses a single `relation' file and determines the
#'     settings for all neighboring stations with
#'     \code{\link[=estimate_AE]{estimate_AE()}} for all event types specified in
#'     \code{event_type}. It fits models to describe translation and retention
#'     processes between neighboring hydrographs, and generates plots
#'     (see vignette for details). Given a file with initial values (see vignette),
#'     predictions are made and visualized in a plot.
#'     Optionally, the results can be written to a directory.
#'     All files need to have the same separator (\code{inputsep}) and
#'     character for decimal points (\code{inputdec}).
#' @param relation_path Character string containing the path of the file where
#'     the relation file is to be read from with
#'     \code{\link[utils:read.csv]{utils::read.csv()}}. The file must contain a
#'     column \code{ID} that contains the gauging station ID. ID's in the file have to be
#'     in order of their location in downstream direction.
#' @param events_path Character string containing the path of the directory
#'     where the event files corresponding to the `relation' file are located. Only
#'     relevant files in this directory will be used, i.e., files that are related to
#'     the `relation' file.
#' @param initial_values_path Character string containing the path of the file
#'     which contains initial values for predictions (see vignette).
#' @param settings_path Character string containing the path where the
#'     settings files are to be read from with
#'     \code{\link[utils:read.csv]{utils::read.csv()}} if
#'     available. The settings files must be in the format of the
#'     output of \code{\link[=peaktrace]{peaktrace()}}. If missing or
#'     incomplete, the settings are determined automatically.
#' @param unique Character string specifying if the potential AEs which
#'     meet the \code{timeLag} and \code{metricLag} condition should
#'     be filtered to contain only unique events using \code{"time"},
#'     i.e., by selecting those where the time difference is smallest
#'     compared to the specified factor of the mean translation time, or using
#'     \code{"metric"}, i.e., by selecting those where the relative
#'     difference in amplitude is smallest (default: \code{"time"}).
#' @param inputdec Character string for decimal points in input data.
#' @param inputsep Field separator character string for input data.
#' @param event_type Vector specifying the event type that is used to identify
#'     event files by their file names
#'     (see \code{\link[hydropeak:get_events]{hydropeak::get_events()}}).
#'     Default: \code{c(2, 4)}, i.e., increasing and decreasing events.
#' @param saveResults A logical. If \code{FALSE} (default), the generated plots
#'     and the estimated settings are not saved. Otherwise the settings are written
#'     to a csv file and the plots are saved as png and pdf files.
#' @param outdir Character string naming a directory where the estimated
#'     settings should be saved to.
#' @param TimeFormat Character string giving the date-time format of the
#'     date-time column in the input data frame (default: "\%Y-\%m-\%d \%H:\%M").
#' @param tz Character string specifying the time zone to be used for the
#'     conversion (default: "Etc/GMT-1").
#' @param formula An object of class \code{\link[stats:formula]{stats::formula()}}
#'     to fit models.
#' @param model Function which specifies the method used for fitting models
#'     (default: \code{\link[stats:lm]{stats::lm()}}). The model class must have a
#'     \code{\link[stats:predict]{stats::predict()}} function.
#' @param FKM_MAX Numeric value that specifies the maximum fkm (see `relation'
#'     file) for which predictions seem valid.
#' @param impute_method Function which specifies the method used for imputing
#'     missing values in initial values based on potential AEs
#'     (default: \code{\link[base:max]{base::max()}}).'
#' @param ... Additional arguments to be passed to the function specified in
#'     argument \code{model}.
#'
#' @return A nested list containing an element for each event type in order as
#'     defined in \code{event_type}. Each element contains again six elements,
#'     namely a data frame of estimated settings, a `\code{gtable}' object that specifies
#'     the combined plot of all stations (plot it with
#'     \code{\link[grid:grid.draw]{grid::grid.draw()}}), a data frame containing
#'     \dQuote{real} AEs (i.e., events where the relative difference in amplitude is within
#'     the estimated cut points), a grid of scatterplots (`gtable` object) for
#'     neighboring hydrographs with a regression line for each metric, a data
#'     frame of results of the model fitting where each row contains the
#'     corresponding stations and metric, the model type (default: "lm"), formula,
#'     coefficients, number of observations and \eqn{R^2}, and a plot of predicted
#'     values based on the \dQuote{initial values}.
#' @export
peaktrace  <- function(relation_path, events_path, initial_values_path, settings_path,
                       unique = c("time", "metric"),
                       inputdec = ".", inputsep = ",", event_type = c(2, 4),
                       saveResults = FALSE, outdir = tempdir(),
                       TimeFormat = "%Y-%m-%d %H:%M", tz = "Etc/GMT-1",
                       formula = y ~ x, model = stats::lm,
                       FKM_MAX = 65, impute_method = base::max, ...) {
  old <- options(scipen = 999)
  on.exit(options(old))

  unique <- match.arg(unique)

  # time matches plus/minus lag
  timeLag <- c(1, 1, 1)

  stopifnot(is.character(relation_path) & length(relation_path) == 1)
  stopifnot(is.character(events_path) & length(events_path) == 1)
  stopifnot(is.atomic(event_type))

  # TimeFormat, tz checked in potential_AE()

  stopifnot(inherits(formula, "formula"))
  stopifnot(is.function(model))
  stopifnot(is.numeric(FKM_MAX) & length(FKM_MAX) == 1)

  relations <- utils::read.csv(file = relation_path, sep = inputsep,
                               dec = inputdec, header = TRUE)

  colnames(relations) <- tolower(colnames(relations))

  stopifnot(all(c("id", "lag") %in% colnames(relations)))

  relations <- relations[!is.na(relations$id), ]

  initials <- utils::read.csv(file = initial_values_path, sep = inputsep,
                              dec = inputdec, header = TRUE)

  events_files <- list.files(events_path, recursive = TRUE)

  # all potential combinations
  events_list <- expand.grid(relations$id, event_type, stringsAsFactors = FALSE)
  events_list <- apply(events_list, 1, paste, collapse = "_")

  # filter event files
  events_files <- events_files[grepl(paste0("^(", paste(events_list,
                                                        collapse = "|"), ")"),
                                     events_files)]

  stopifnot(length(events_files) == length(events_list))

  results_list <- vector(mode = "list", length = length(event_type))

  # iterate over event types
  for (e in seq_len(length(event_type))) {
    # empty lists
    plot_list_threshold <- vector(mode = "list", length = (nrow(relations) - 1))
    settings_list <- vector(mode = "list", length = (nrow(relations) - 1))
    real_AE_list <- vector(mode = "list", length = (nrow(relations) - 1))

    if (!missing(settings_path)) {
        file <- file.path(settings_path, paste0("Q_event_", event_type[e], "_AMP-LAG_settings.csv"))
        if (file.exists(file)) {
            settings_list_provided <- utils::read.csv(file = file, sep = ",", dec = ".")
        } else {
            settings_list_provided <- NULL
        }
    }

    # iterate over stations
    for (i in seq_len(nrow(relations) - 1)) {
      STATIONS <- relations$station[i + 0:1]

      # load data
      Sx <- file.path(events_path, events_files[grepl(paste0("^",
                                                             relations$id[i],"_",
                                                             event_type[e]),
                                                      events_files)])
      stopifnot(length(Sx) == 1)

      Sy <- file.path(events_path, events_files[grepl(paste0("^",
                                                             relations$id[i + 1],"_",
                                                             event_type[e]),
                                                      events_files)])
      stopifnot(length(Sy) == 1)

      Sx <- utils::read.csv(Sx)
      Sx$station <- rep(STATIONS[1], nrow(Sx))
      Sx$Time <- as.POSIXct(strptime(Sx$Time, format = TimeFormat, tz = tz))
      Sy <- utils::read.csv(Sy)
      Sy$station <- rep(STATIONS[2], nrow(Sy))
      Sy$Time <- as.POSIXct(strptime(Sy$Time, format = TimeFormat, tz = tz))

      # keep only real AEs previously identified
      if (i > 1) {
          IDs <- real_AE_list[[i-1]][, c("id.y", "time.y")] |>
            dplyr::mutate(time.y = time.y)
          Sx <- dplyr::inner_join(Sx, IDs, c("ID" = "id.y", "Time" = "time.y"))
      }

      # difference between 2 stations (translation time)
      relation <- relations[vapply(STATIONS, grep, relations$station,
                                   FUN.VALUE = numeric(1)), ]

      # save all together afterwards
      settings <- NULL
      if (!missing(settings_path) && !is.null(settings_list_provided)) {
          settings <- settings_list_provided |>
              dplyr::filter(station.x == STATIONS[1],
                            station.y == STATIONS[2],
                            !is.na(metric))
          if (nrow(settings) == 0) {
              settings <- NULL
          }
      }

      results_list[[e]] <- estimate_AE(Sx = Sx, Sy = Sy, relation = relation,
                                       timeLag = timeLag, unique = unique,
                                       TimeFormat = TimeFormat, tz = tz,
                                       settings = settings)

      settings_list[[i]] <- results_list[[e]]$settings
      plot_list_threshold[[i]] <- results_list[[e]]$plot_threshold
      real_AE_list[i] <- list(results_list[[e]]$real_AE)

      if (is.null(real_AE_list[[i]])) {
        break
      } else {
        initials <- impute_initials(real_AE_list[[i]], initials, impute_method)
      }
    }

    results_list[[e]]$settings <- do.call(rbind, settings_list)

    plot_list_threshold <- plot_list_threshold[!sapply(plot_list_threshold, is.null)]
    results_list[[e]]$plot_threshold <- gridExtra::arrangeGrob(grobs = plot_list_threshold,
                                                               nrow = 1)
    results_list[[e]]$real_AE <- do.call(rbind, real_AE_list)

    routing_list <- routing(real_AE = results_list[[e]][[3]], initials = initials,
                            relation = relations, formula = formula,
                            model = model, FKM_MAX = FKM_MAX, ...)

    results_list[[e]]$plot_scatter <- routing_list$plot_scatter
    results_list[[e]]$models <- routing_list$models
    results_list[[e]]$plot_predict <- routing_list$plot_predict

    if (saveResults) {
      # save settings
      writeCSV(results_list[[e]]$settings, outdir = outdir, metric = "AMP",
               event = event_type[e], postfix = "settings")

      # save plot threshold
      ggplot2::ggsave(file.path(outdir, paste0("plot_thresholds_",
                                               event_type[e], ".pdf")),
                      plot = results_list[[e]]$plot_threshold,
                      width = 15, height = 5)

      ggplot2::ggsave(file.path(outdir, paste0("plot_thresholds_",
                                               event_type[e], ".png")),
                      plot = results_list[[e]]$plot_threshold,
                      width = 15, height = 5)

      # save real AE
      writeCSV(results_list[[e]]$real_AE, outdir = outdir, metric = "AMP",
               event = event_type[e], postfix = "real_AE")

      # save plot scatter
      ggplot2::ggsave(file.path(outdir, paste0("plot_scatter_",
                                               event_type[e], ".pdf")),
                      plot = results_list[[e]]$plot_scatter,
                      width = 10, height = 12)

      ggplot2::ggsave(file.path(outdir, paste0("plot_scatter_",
                                               event_type[e], ".png")),
                      plot = results_list[[e]]$plot_scatter,
                      width = 10, height = 12)

      # save plot predict
      ggplot2::ggsave(file.path(outdir, paste0("plot_predict_",
                                               event_type[e], ".pdf")),
                      plot = results_list[[e]]$plot_predict,
                      width = 15, height = 5)

      ggplot2::ggsave(file.path(outdir, paste0("plot_predict_",
                                               event_type[e], ".png")),
                      plot = results_list[[e]]$plot_predict,
                      width = 15, height = 5)
    }
  }
  names(results_list) <- event_type
  invisible(results_list)
}


#' Estimate Models and Make Predictions
#' @description Performs the \dQuote{routing} procedure, i.e., based on associated events,
#'     it uses (linear) models to describe translation and retention processes
#'     between neighboring hydrographs.
#' @param real_AE Data frame that contains real AEs of two neighboring hydrographs
#'     estimated with \code{\link[=estimate_AE]{estimate_AE()}}.
#' @param initials Data frame that contains initial values for predictions
#'     (see vignette).
#' @param relation Data frame that contains the relation between upstream and
#'     downstream hydrograph. Must only contain two rows (one for each hydrograph)
#'     in order of their location in downstream direction.
#'     See the appended example data \code{relation.csv} or vignette for details on
#'     the structure. See \code{\link[=get_lag]{get_lag()}} for further information
#'     about the relation and the lag between the hydrographs.
#' @param formula An object of class \code{\link[stats:formula]{stats::formula()}}
#'     to fit models.
#' @param model Function which specifies the method used for fitting models
#'     (default: \code{\link[stats:lm]{stats::lm()}}). The model class must have a
#'     \code{\link[stats:predict]{stats::predict()}} function.
#' @param FKM_MAX Numeric value that specifies the maximum fkm (see relation
#'     file) for which predictions seem valid.
#' @param ... Additional arguments to be passed to the function specified in
#'     argument \code{model}.
#' @return A nested list containing a grid of scatterplots (`gtable` object) for
#'     neighboring hydrographs with a regression line for each metric, a data
#'     frame of results of the model fitting where each row contains the
#'     corresponding stations and metric, the model type (default: "lm"), formula,
#'     coefficients, number of observations and \eqn{R^2}, and a plot of predicted
#'     values based on the \dQuote{initial values}.
#' @export
routing <- function(real_AE, initials, relation, formula = y ~ x,
                    model = stats::lm, FKM_MAX = 65, ...) {

  colnames(relation) <- tolower(colnames(relation))
  stopifnot(all(c("id", "lag") %in% colnames(relation)))

  colnames(real_AE) <- tolower(colnames(real_AE))
  colnames(initials) <- tolower(colnames(initials))

  routing_list <- vector(mode = "list", length = 3)

  ## keep only real AEs with > 1 events per group
  real_AE_split <- real_AE |>
    dplyr::group_by(id.x, id.y) |>
    dplyr::group_split()

  keep <- which(lapply(real_AE_split, nrow) > 1)
  delete <- which(lapply(real_AE_split, nrow) <= 1)

  if (length(delete) > 0) {
    lapply(real_AE_split[delete], function(x)
               warning(paste("Less than 2 AEs:",
                             paste0("Event type ", unique(x$event_type.x), ","),
                             "ID Sx:", unique(x$id.x), "ID Sy:", unique(x$id.y))))
  }

  real_AE <- do.call("rbind", real_AE_split[keep])

  ## scatter plots
  routing_list[[1]] <- real_AE |>
    dplyr::group_by(id.x, id.y) |>
    plot_scatter()

  ## estimate models
  real_AE_long <- real_AE[, !(names(real_AE) %in% c("diff_metric"))]
  real_AE_long <- reshape(real_AE_long)

  split <- real_AE_long |>
    dplyr::group_by(id.x, id.y) |>
    dplyr::group_split()

  model_list <- lapply(split, function(x) fit_model(formula = formula, data = x, model = model, ...))

  models <- do.call(rbind, lapply(model_list, function(x) x[["models"]]))

  routing_list[[2]] <- do.call(rbind, lapply(model_list, function(x) x[["models.df"]]))

  ## omit models not used for predictions
  routing_list[[2]] <- routing_list[[2]] |>
    dplyr::left_join(initials |>
                     dplyr::transmute(station = station,
                                      metric = tolower(metric)) |>
                     unique(), by = "metric") |>
    dplyr::filter(station.x >= station) |>
    dplyr::select(-station)

  ## predictions
  predictions <- get_predictions(models = models, initials = initials)
  predictions$metric <- tolower(predictions$metric)

  # add fkm
  predictions <- merge(predictions, relation[, c("station", "fkm")],
                       by = "station")

  # add n and r2
  predictions <- merge(predictions,
                       unique(routing_list[[2]][, c("station.y", "n")]),
                       by.x = "station", by.y = "station.y", all = TRUE)

  predictions <- merge(predictions,
                       unique(routing_list[[2]][, c("station.y", "metric", "r2")]),
                       by.x = c("station", "metric"),
                       by.y = c("station.y", "metric"), all = TRUE)

  # add metric labels
  metrics <- data.frame(new = c("AMP", "MAFR", "MEFR", "DUR", "RATIO"),
                        unit = c("m\u00b3/s",
                                 "(m\u00b3/s)/ts",
                                 "(m\u00b3/s)/ts", "ts", "-"))

  predictions <- predictions |>
    dplyr::mutate(unit = factor(toupper(metric), metrics$new, metrics$unit),
                  metric = factor(toupper(metric),
                                  metrics$new,
                                  with(metrics, paste0(new, " (", unit, ")"))))

  # Maximum fkm for which predictions seem valid.
  predictions <- subset(predictions, fkm <= FKM_MAX)

  routing_list[[3]] <- plot_predict(predictions)

  names(routing_list) <- c("plot_scatter", "models", "plot_predict")

  routing_list
}


#' Extract Associated Events
#' @description For given relation and event data return the
#'     associated events which comply with the conditions specified in
#'     the settings.
#' @param relation_path Character string containing the path of the file where
#'     the relation file is to be read from with
#'     \code{\link[utils:read.csv]{utils::read.csv()}}. The file must contain a
#'     column \code{ID} that contains the gauging station ID's in the file have to be
#'     in order of their location in downstream direction.
#' @param events_path Character string containing the path of the directory
#'     where the event files corresponding to the `relation` file are located. Only
#'     relevant files in this directory will be used, i.e., files that are related to
#'     the `relation` file.
#' @param settings_path Character string containing the path of the file where
#'     the settings file is to be read from with
#'     \code{\link[utils:read.csv]{utils::read.csv()}}. The file must be in the format
#'     of the output of \code{\link[=peaktrace]{peaktrace()}}.
#' @param unique Character string specifying if the potential AEs which
#'     meet the \code{timeLag} and \code{metricLag} condition should
#'     be filtered to contain only unique events using \code{"time"},
#'     i.e., by selecting those where the time difference is smallest
#'     compared to the specified factor of the mean translation time, or using
#'     \code{"metric"}, i.e., by selecting those where the relative
#'     difference in amplitude is smallest (default: \code{"time"}).
#' @param inputdec Character string for decimal points in input data.
#' @param inputsep Field separator character string for input data.
#' @param saveResults A logical. If \code{FALSE} (default), the extracted AEs are not saved.
#'     Otherwise the extracted AEs are written to a csv file.
#' @param outdir Character string naming a directory where the extraced AEs
#'     should be saved to.
#' @param TimeFormat Character string giving the date-time format of the
#'     date-time column in the input data frame (default: "\%Y-\%m-\%d \%H:\%M").
#' @param tz Character string specifying the time zone to be used for the
#'     conversion (default: "Etc/GMT-1").
#'
#' @return A data frame containing \dQuote{real} AEs (i.e., events
#'     where the time differences and the relative difference in
#'     amplitude is within the limits and cut points provided by the
#'     file in \code{settings_path}). If no AEs can be found between the first
#'     two neighboring stations, \code{NULL} is returned. Otherwise the function
#'     returns all \dQuote{real} AEs that could be found along the river section
#'     specified in the file from \code{relation_path}. A warning is issued when
#'     the extraction is stopped early and shows the \code{IDs} for which no
#'     AEs are determined.
#' @export
#'
#' @examples
#' relation_path <- system.file("testdata", "relation.csv", package = "hydroroute")
#' events_path <- system.file("testdata", "Events", package = "hydroroute")
#' settings_path <- system.file("testdata", "Q_event_2_AMP-LAG_aut_settings.csv",
#'                                    package = "hydroroute")
#' real_AE <- extract_AE(relation_path, events_path, settings_path)
extract_AE <- function(relation_path, events_path, settings_path,
                       unique = c("time", "metric"),
                       inputdec = ".", inputsep = ",",
                       saveResults = FALSE, outdir = tempdir(),
                       TimeFormat = "%Y-%m-%d %H:%M", tz = "Etc/GMT-1") {

  unique <- match.arg(unique)

  stopifnot(is.character(relation_path) & length(relation_path) == 1)
  stopifnot(is.character(events_path) & length(events_path) == 1)
  stopifnot(is.character(settings_path) & length(settings_path) == 1)

  stopifnot(is.character(TimeFormat), length(TimeFormat) == 1)
  stopifnot(is.character(tz), length(tz) == 1)

  relations <- utils::read.csv(file = relation_path, sep = inputsep,
                               dec = inputdec, header = TRUE)

  colnames(relations) <- tolower(colnames(relations))

  stopifnot(all(c("id", "lag") %in% colnames(relations)))

  relations <- relations[!is.na(relations$id), ]

  settings <- utils::read.csv(file = settings_path, sep = inputsep,
                               dec = inputdec, header = TRUE)

  colnames(settings) <- tolower(colnames(settings))
  stopifnot(all(c("station.x", "station.y", "bound", "lag", "metric") %in% colnames(settings)))

  events_files <- list.files(events_path, recursive = TRUE)

  # all potential combinations
  pattern <- "Q_event_(0|1|2|3|4|5)_"

  event_type <- regmatches(basename(settings_path), regexec(pattern, basename(settings_path)))[[1]][2]
  events_list <- expand.grid(relations$id, event_type, stringsAsFactors = FALSE)
  events_list <- apply(events_list, 1, paste, collapse = "_")

  # filter event files
  events_files <- events_files[grepl(paste0("^(", paste(events_list,
                                                        collapse = "|"), ")"),
                                     events_files)]

  stopifnot(length(events_files) == length(events_list))

  # empty list
  real_AE_list <- vector(mode = "list", length = (nrow(relations) - 1))

  # iterate over stations
  for (i in seq_len(nrow(relations) - 1)) {
    STATIONS <- relations$station[i + 0:1]
    real_AE <- NULL

    # load data
    Sx <- file.path(events_path, events_files[grepl(paste0("^",
                                                           relations$id[i],"_",
                                                           event_type),
                                                    events_files)])
    stopifnot(length(Sx) == 1)

    Sy <- file.path(events_path, events_files[grepl(paste0("^",
                                                           relations$id[i + 1],"_",
                                                           event_type),
                                                    events_files)])
    stopifnot(length(Sy) == 1)

    Sx <- utils::read.csv(Sx)
    Sx$station <- rep(STATIONS[1], nrow(Sx))
    Sy <- utils::read.csv(Sy)
    Sy$station <- rep(STATIONS[2], nrow(Sy))

    # keep only real AEs previously identified
    if (i > 1) {
        IDs <- real_AE_list[[i-1]][, c("id.y", "time.y")] |>
          dplyr::mutate(time.y = as.character(time.y))
        Sx <- dplyr::inner_join(Sx, IDs, c("ID" = "id.y", "Time" = "time.y"))
    }

    # difference between 2 stations (translation time)
    relation <- relations[vapply(STATIONS, grep, relations$station,
                                 FUN.VALUE = numeric(1)), ]

    setting <- settings |>
      dplyr::filter(station.x == STATIONS[1], station.y == STATIONS[2])

    # lag_opt
    lag <- ifelse(is.na(relation$lag), "00:00", as.character(relation$lag))
    lag_opt <- diff(as.difftime(lag, format = "%H:%M", units = "mins")) *
        setting$lag[2]

    events <- merge_time(Sx = Sx, Sy = Sy, relation = relation, timeLag = setting$lag,
                         TimeFormat = TimeFormat, tz = tz)

    potential_AE <- potential_AE(Sx = Sx, Sy = Sy, relation = relation,
                                 timeLag = setting$lag, unique = unique,
                                 TimeFormat = TimeFormat, tz = tz)

    if (!is.null(potential_AE)) {
      potential_AE_subset <- potential_AE |>
      dplyr::filter(metric == "amp")

      # real AEs
      # filter events within estimated cut points (settings)
      events <- events |>
        dplyr::inner_join(potential_AE_subset |>
                          dplyr::select(station.x, time.x, station.y,
                                        time.y, diff_metric),
                          by = c("station.x", "time.x", "station.y", "time.y"))

      real_AE <- events |>
        dplyr::filter(dplyr::between(diff_metric + 1, setting$metric[1],
                                     setting$metric[3]))

      # save all together afterwards
      real_AE_list[[i]] <- real_AE
    } else {
      warning(paste0("No potential AEs between ", unique(Sx$ID), " and ", unique(Sy$ID)))
      break
    }
  }

  real_AE <- do.call(rbind, real_AE_list)

  if (saveResults) {
    writeCSV(real_AE, outdir = outdir, metric = "AMP",
             event = event_type, postfix = "real_AE")
  }

  invisible(real_AE)
}


#' Reshape Q Data
#' @description Reshapes data frame from \dQuote{wide} format to \dQuote{long} format.
#' @param x Merged data frame constructed in \code{\link[=merge_time]{merge_time()}}.
#'
#' @return Data frame in \dQuote{wide} format.
#' @keywords internal
#' @noRd
reshape <- function(x) {
  long <- reshape2::melt(x, id = c("station.x", "time.x", "id.x", "event_type.x",
                                   "station.y", "time.y", "id.y", "event_type.y"))
  stopifnot(is.data.frame(x))
  variables <- strsplit(as.character(long$variable), "\\.")
  long$metric <- vapply(variables, "[", 1, FUN.VALUE = character(1))
  long$location <- vapply(variables, "[", 2, FUN.VALUE = character(1))
  wide <- reshape2::dcast(station.x + id.x + time.x + event_type.x + station.y +
                            id.y + time.y + event_type.y + metric ~ location,
                          data = long)
  wide
}
