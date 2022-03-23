#' Get Lag
#'
#' @description Given a data frame (time series) of stage measurements
#'     and a vector of gauging station ID's in order of their location
#'     in downstream direction, the lag (the amount of passing time
#'     between two gauging stations) is estimated based on the
#'     cross-correlation function (ccf) of the time series of two
#'     adjacent gauging stations
#'     (\code{\link[stats:ccf]{stats::ccf()}}).  To ensure that the
#'     same time period is used for every gauging station,
#'     intersecting time steps are determined. These time steps are
#'     used to estimate the lags. The result of
#'     \code{\link[stats:ccf]{stats::ccf()}} is rounded to four
#'     decimals before selecting the optimal time lag so that minimal
#'     differences are neglected.  If there are multiple time steps
#'     with the highest correlation, the smallest time step is
#'     considered. If the highest correlation corresponds to a zero
#'     lag or positive lag (result should usually be negative as
#'     measurements at the lower gauge are later recorded as
#'     measurements at the upper gauge), a time step of length 1 is
#'     selected and a warning message is generated.
#'
#' @param Q Data frame (time series) of stage measurements which
#'     contains at least a column with the gauging station ID
#'     (default: column index 1), a column with date-time values in
#'     character representation (default: column index 2) and a column
#'     with flow rates (default: column index 3).  If the column
#'     indices differ from \code{c(1, 2, 3)}, they have to be
#'     specified in the \code{cols} argument in the format \code{c(i,
#'     j, k)}.
#' @param relation A character vector containing the gauging station
#'     ID's in order of their location in downstream direction.
#' @param steplength Numeric value that specifies the length between
#'     time steps in minutes (default: \code{15} minutes). As time
#'     steps have to be equispaced, this is used by
#'     \code{\link[hydropeak:flow]{hydropeak::flow()}} to get a
#'     compatible format and fill missing time steps with \code{NA}.
#' @param lag.max Numeric value that specifies the maximum lag at
#'     which to calculate the ccf in
#'     \code{\link[stats:ccf]{stats::ccf()}}.
#' @param na.action Function to be called to handle missing values in
#'     \code{\link[stats:ccf]{stats::ccf()}} (default:
#'     \code{na.pass}).
#' @param mc.cores Number of cores to use with
#'     \code{\link[parallel:mclapply]{parallel::mclapply()}}. On
#'     Windows, this is set to 1.
#' @param tz Character string specifying the time zone to be used for
#'     internal conversion (default: \code{Etc/GMT-1}).
#' @param format Character string giving the date-time format of the
#'     date-time column in the input data frame \code{Q}. This is
#'     passed to \code{\link[hydropeak:flow]{hydropeak::flow()}}, to
#'     get a compatible format (default: \code{YYYY.mm.dd HH:MM}).
#' @param cols Integer vector specifying column indices in \code{Q}.
#'     The default indices are 1 (ID), 2 (date-time) and 3 (flow rate,
#'     Q).  This is passed to
#'     \code{\link[hydropeak:flow]{hydropeak::flow()}}.
#'
#' @return A character vector which contains the estimated cumulative
#'     lag between neighboring gauging stations in the format
#'     \code{HH:MM}.
#' @importFrom stats na.pass
#' @export
#'
#' @examples
#' Q_path <- system.file("testdata", "Q.csv", package = "hydroroute")
#' Q <- utils::read.csv(Q_path)
#'
#' relation_path <- system.file("testdata", "relation.csv",
#'                             package = "hydroroute")
#' relation <- utils::read.csv(relation_path)
#' # from relation data frame
#' get_lag(Q, relation$ID, format = "%Y-%m-%d %H:%M", tz = "Etc/GMT-1")
#'
#' # station ID's in downstream direction as vector
#' relation <- c("100000", "200000", "300000", "400000")
#' get_lag(Q, relation, format = "%Y-%m-%d %H:%M", tz = "Etc/GMT-1")
get_lag <- function(Q, relation, steplength = 15, lag.max = 20,
                    na.action = na.pass, mc.cores = getOption("mc.cores", 2L),
                    tz = "Etc/GMT-1", format = "%Y.%m.%d %H:%M",
                    cols = c(1, 2, 3)) {

  stopifnot(length(mc.cores) == 1,
            is.integer(mc.cores) | is.numeric(mc.cores))
  if (.Platform$OS.type == "windows" && mc.cores > 1L) {
    mc.cores <- 1L
  }

  stopifnot(is.data.frame(Q))
  stopifnot(is.atomic(relation))
  stopifnot(length(steplength) == 1,
            is.numeric(steplength) | is.integer(steplength))
  stopifnot(length(lag.max) == 1,
            is.numeric(lag.max) | is.integer(lag.max))

  relation <- relation[!is.na(relation)]
  relation_ID <- as.character(relation)
  len <- length(relation)

  Q$ID <- as.character(Q$ID)

  # use only the IDs specified in relation
  Q <- Q[Q$ID %in% relation_ID, ]

  # get compatible format and fill missing time steps with NA
  if (!inherits(Q, "flow")) {
    Q <- hydropeak::flow(Q, tz = tz, format = format, cols = cols,
                         steplength = steplength)
  }

  if (!inherits(Q$Time, "POSIXct")) {
    Q$Time <- as.POSIXct(Q$Time, origin = "1970-01-01", tz = tz)
  }

  # preserve relation order
  Q$ID <- factor(Q$ID, levels = relation_ID)

  Q <- split(Q, ~ID)

  # empty lag vector to store result from stats::ccf
  lag <- vector(mode = "numeric", length = len - 1)

  for (i in seq_len(len - 1)) {
    # determine equal time steps in neighboring stations
    intersect_ts <- as.POSIXct(intersect(Q[[i]]$Time, Q[[i + 1]]$Time),
                               origin = "1970-01-01", tz = tz)

    # subset intersecting time steps
    Q1 <- Q[[i]][Q[[i]]$Time %in% intersect_ts, ]
    Q2 <- Q[[i + 1]][Q[[i + 1]]$Time %in% intersect_ts, ]

    # cross correlation between flow values
    ccfvalues <- stats::ccf(x = Q1$Q, y = Q2$Q, lag.max = lag.max, plot = FALSE,
                            na.action = na.action)

    # round to not consider minimal differences
    ccfcorr <- round(ccfvalues$acf[,,1], digits = 4)
    max_cor <- which(ccfcorr == max(ccfcorr))

    if (length(max_cor) == 1 && ccfvalues$lag[max_cor] < 0) {
      lag[i] <- ccfvalues$lag[max_cor]

    } else if (length(max_cor) > 1 && all(ccfvalues$lag[max_cor] < 0)) {
      # multiple time steps with highest correlation, get smallest time step
      lag[i] <- ccfvalues$lag[max_cor[length(max_cor)]]

    } else {
      # warning if highest correlation has 0 or positive time step
      # --> assume time step 1
      if (length(max_cor) >= 1 && ccfvalues$lag[max_cor[1]] == 0) {
        warning("ID ", Q1$ID[1], " and ID ", Q2$ID[1], ": Lag 0")
      } else if (length(max_cor) >= 1 && ccfvalues$lag[max_cor[1]] > 0) {
        warning("ID ", Q1$ID[1], " and ID ", Q2$ID[1], ": Positive lag")
      }
      lag[i] <- 1
    }
  }

  time <- cumsum(abs(lag)) * steplength

  hhmm <- sapply(time, function(x) min_to_hhmm(x))

  # return lag, first lag is always 00:00
  c("00:00", hhmm)
}


#' Get Lag from Input File
#'
#' @description Given a file path it reads a data frame (time series)
#'     of stage measurements which combines several ID's and calls
#'     \code{\link[=get_lag]{get_lag()}}. The relation (ID's) of
#'     gauging stations is read from a file (provide file path). Make
#'     sure that the file with \code{Q} data and the relation file
#'     have the same separator (\code{inputsep}) and character for
#'     decimal points (\code{inputdec}). Gauging station ID's have to
#'     be in order of their location in downstream direction. The
#'     resulting lag is appended to the relation file. This can be
#'     saved to a file.
#'
#' @param Q_file Data frame or character string. If it is a data
#'     frame, it corresponds to the \code{Q} data frame in
#'     \code{\link[=get_lag]{get_lag()}}.  It contains at least a
#'     column with the gauging station ID (default: column index 1), a
#'     column with date-time values in character representation
#'     (default: column index 2) and a column with flow rates
#'     (default: column index 3). If the column indices differ from
#'     \code{c(1, 2, 3)}, they have to be specified as \code{cols}
#'     argument in the format \code{c(i, j, k)}. If it is a character
#'     string, it contains the path to the corresponding file which is
#'     then read within the function with
#'     \code{\link[utils:read.csv]{utils::read.csv()}}.
#' @param relation_file A character string containing the path to the
#'     relation file. It is read within the function with
#'     \code{\link[utils:read.csv]{utils::read.csv()}}.  The file must
#'     contain a column \code{ID} that contains the gauging station
#'     ID's in order of their location in downstream direction. The
#'     lag will then be appended as column to the data frame. For more
#'     details on the relation file, see the vignette.
#' @param steplength Numeric value that specifies the length between
#'     time steps in minutes (default: \code{15} minutes). As time
#'     steps have to be equispaced, this is used by
#'     \code{\link[hydropeak:flow]{hydropeak::flow()}} to get a
#'     compatible format and fill missing time steps with \code{NA}.
#' @param lag.max Maximum lag at which to calculate the ccf in
#'     \code{\link[stats:ccf]{stats::ccf()}}.
#' @param na.action Function to be called to handle missing values in
#'     \code{\link[stats:ccf]{stats::ccf()}} (default:
#'     \code{na.pass}).
#' @param tz Character string specifying the time zone to be used for
#'     internal conversion (default: \code{Etc/GMT-1}).
#' @param format Character string giving the date-time format of the
#'     date-time column in the input data frame \code{Q}. This is
#'     passed to \code{\link[hydropeak:flow]{hydropeak::flow()}}, to
#'     get a compatible format (default: \code{YYYY.mm.dd HH:MM}).
#' @param cols Integer vector specifying column indices in the input
#'     data frame which contain gauging station ID, date-time and flow
#'     rate to be renamed. The default indices are 1 (ID), 2
#'     (date-time) and 3 (flow rate, Q).
#' @param inputsep Character string for the field separator in input
#'     data.
#' @param inputdec Character string for decimal points in input data.
#' @param save A logical. If \code{FALSE} (default) the lag, appended
#'     to the relation file, is not written to a file, otherwise it is
#'     written to \code{outfile}.
#' @param outfile A character string naming a file path where the
#'     output file should be written to.
#' @param mc.cores Number of cores to use with
#'     \code{\link[parallel:mclapply]{parallel::mclapply()}}. On
#'     Windows, this will be set to 1.
#' @param overwrite A logical. If \code{FALSE} (default), it produces
#'     an error if a \code{LAG} column already exists in the
#'     \code{relation} file. Otherwise, it overwrites an existing
#'     column.
#'
#' @return Returns invisible the data frame of the relation data with
#'     the estimated cumulative lag between neighboring gauging
#'     stations' lag in the format \code{HH:MM} appended.
#' @export
#'
#' @examples
#' Q_file <- system.file("testdata", "Q.csv", package = "hydroroute")
#' relation_file <- system.file("testdata", "relation.csv",
#'                              package = "hydroroute")
#' get_lag_file(Q_file, relation_file, inputsep = ",", inputdec = ".",
#'              format = "%Y-%m-%d %H:%M", save = FALSE, overwrite = TRUE)
#'
#' Q_file <- read.csv(Q_file)
#' get_lag_file(Q_file, relation_file, inputsep = ",", inputdec = ".",
#'              format = "%Y-%m-%d %H:%M", save = FALSE, overwrite = TRUE)
get_lag_file <- function(Q_file, relation_file, steplength = 15, lag.max = 20,
                         na.action = na.pass, tz = "Etc/GMT-1",
                         format = "%Y.%m.%d %H:%M", cols = c(1, 2, 3),
                         inputsep = ";", inputdec = ".", save = FALSE,
                         outfile = file.path(tempdir(), "relation.csv"),
                         mc.cores = getOption("mc.cores", 2L),
                         overwrite = FALSE) {


  # Q_file can data frame or path to be read from
  stopifnot(is.data.frame(Q_file) | (is.character(Q_file) & length(Q_file) == 1))

  stopifnot(is.character(relation_file) & length(relation_file) == 1)
  stopifnot(is.logical(save) & length(save) == 1)
  stopifnot(is.logical(overwrite) & length(overwrite) == 1)
  stopifnot(is.character(outfile) & length(outfile) == 1)

  if (!is.data.frame(Q_file)) {
    Q_file <- utils::read.csv(Q_file, sep = inputsep, dec = inputdec)
  }

  relation <- utils::read.csv(relation_file, sep = inputsep, dec = inputdec)

  # existing lag and should not be overwritten
  stopifnot(is.null(relation$LAG) | overwrite)

  # omit empty rows
  relation <- relation[!is.na(relation$ID), ]

  # get lag from Q_file
  lag <- get_lag(Q = Q_file, relation = relation$ID, steplength = steplength,
                 na.action = na.action, mc.cores = mc.cores, tz = tz,
                 format = format, cols = cols, )

  # append lag as column
  relation$LAG <- lag

  if (save) {
    utils::write.csv(x = relation, file = outfile, row.names = FALSE)
  }

  invisible(relation)
}

#' Get Lag from Input Directory
#'
#' @description Given a file path it reads a data frame (time series)
#'     of stage measurements. For each \code{relation} file in the
#'     provided directory path it calls
#'     \code{\link[=get_lag_file]{get_lag_file()}}.  Make sure that
#'     the file with Q data and the relation files have the same
#'     separator (\code{inputsep}) and character for decimal points
#'     (\code{inputdec}). Gauging station ID's in the \code{relation}
#'     files have to be in order of their location in downstream
#'     direction. The resulting lags are appended to the relation
#'     files. The resulting list of relation files can be returned and
#'     each relation file can be saved to its input path.
#'
#' @param Q Data frame or character string. If it is a data frame, it
#'     corresponds to the \code{Q} data frame in
#'     \code{\link[=get_lag]{get_lag()}}.  It contains at least a
#'     column with the gauging station ID (default: column index 1), a
#'     column with date-time values in character representation
#'     (default: column index 2) and a column with flow rates
#'     (default: column index 3). If the column indices differ from
#'     \code{c(1, 2, 3)}, they have to be specified as \code{cols}
#'     argument in the format \code{c(i, j, k)}. If it is a character
#'     string, it contains the path to the corresponding file which is
#'     then read within the function with
#'     \code{\link[utils:read.csv]{utils::read.csv()}}.
#' @param relation A character string containing the path to the
#'     directory where the relation files are located. They are read
#'     within the function with
#'     \code{\link[utils:read.csv]{utils::read.csv()}}.
#' @param steplength Numeric value that specifies the length between
#'     time steps in minutes (default: \code{15} minutes). As time
#'     steps have to be equispaced, this is used by
#'     \code{\link[hydropeak:flow]{hydropeak::flow()}} to get a
#'     compatible format and fill missing time steps with \code{NA}.
#' @param lag.max Maximum lag at which to calculate the ccf in
#'     \code{\link[stats:ccf]{stats::ccf()}}.
#' @param na.action Function to be called to handle missing values in
#'     \code{\link[stats:ccf]{stats::ccf()}} (default:
#'     \code{na.pass}).
#' @param tz Character string specifying the time zone to be used for
#'     internal conversion (default: \code{Etc/GMT-1}).
#' @param format Character string giving the date-time format of the
#'     date-time column in the input data frame \code{Q}. This is
#'     passed to \code{\link[hydropeak:flow]{hydropeak::flow()}}, to
#'     get a compatible format (default: \code{YYYY.mm.dd HH:MM}).
#' @param cols Integer vector specifying column indices in the input
#'     data frame which contain gauging station ID, date-time and flow
#'     rate to be renamed. The default indices are 1 (ID), 2
#'     (date-time) and 3 (flow rate, Q).
#' @param inputsep Field separator character string for input data.
#' @param inputdec Character string for decimal points in input data.
#' @param relation_pattern Character string containing a regular
#'     expression to filter \code{relation} files (default: \code{relation}, to
#'     filter files that contain \code{relation} with no restriction) (see
#'     \code{\link[base:grep]{base::grep()}}).
#' @param save A logical. If \code{FALSE} (default) the lag, appended
#'     to the relation file, overwrites the original \code{relation}
#'     input file.
#' @param mc.cores Number of cores to use with
#'     \code{\link[parallel:mclapply]{parallel::mclapply()}}. On
#'     Windows, this will be set to 1.
#' @param overwrite A logical. If \code{FALSE} (default), it produces
#'     an error if a \code{LAG} column already exists in the
#'     \code{relation} file. Otherwise, it overwrites an existing
#'     column.
#'
#' @return Returns invisible a list of data frames where each list
#'     element represents a \code{relation} file from the input
#'     directory. Optionally, the data frames overwrite the existing
#'     \code{relation} files with the appended \code{LAG} column.
#' @export
#'
#' @examples
#' Q_file <- system.file("testdata", "Q.csv", package = "hydroroute")
#' relations_path <- system.file("testdata", package = "hydroroute")
#' lag_list <- get_lag_dir(Q_file, relations_path, inputsep = ",",
#'                         inputdec = ".", format = "%Y-%m-%d %H:%M",
#'                         overwrite = TRUE)
#' lag_list
get_lag_dir <- function(Q, relation, steplength = 15, lag.max = 20,
                        na.action = na.pass, tz = "Etc/GMT-1",
                        format = "%Y.%m.%d %H:%M", cols = c(1, 2, 3),
                        inputsep = ",", inputdec = ".",
                        relation_pattern = "relation", save = FALSE,
                        mc.cores = getOption("mc.cores", 2L),
                        overwrite = FALSE) {

  if (.Platform$OS.type == "windows" & mc.cores > 1L) {
    mc.cores <- 1L
  }
  # Q_file can data frame or path to be read from
  stopifnot(is.data.frame(Q) | (is.character(Q) & length(Q) == 1))

  stopifnot(is.character(relation_pattern) & length(relation_pattern) == 1)

  # single Q data file
  if (!is.data.frame(Q)) {
    Q_file <- utils::read.csv(Q, sep = inputsep, dec = inputdec)
  }

  file_list <- dir(relation, recursive = TRUE, full.names = TRUE, no.. = TRUE)
  relation_list <- grep(relation_pattern, file_list, value = TRUE)

  result <- parallel::mclapply(relation_list,
                               function(x) get_lag_file(Q_file = Q_file,
                                                        relation_file = x,
                                                        steplength = steplength,
                                                        lag.max = lag.max,
                                                        na.action = na.action,
                                                        mc.cores = mc.cores,
                                                        tz = tz,
                                                        format = format,
                                                        cols = cols,
                                                        inputsep = inputsep,
                                                        inputdec = inputdec,
                                                        save = save,
                                                        outfile = x,
                                                        overwrite = overwrite),
                               mc.cores = mc.cores)

   invisible(result)
}

#' Convert Minutes to HH:MM Format
#'
#' @param minutes A positive integer of minutes that has to be
#'     converted to HH:MM format.
#' @return A character that represents the provided minutes in HH:MM
#'     format.
#' @keywords internal
#' @noRd
min_to_hhmm <- function(minutes) {
  stopifnot(length(minutes) == 1,
            as.numeric(minutes) | as.integer(minutes))

  mm <- minutes %% 60
  hh <- as.integer(minutes / 60)

  if (mm < 10) {
    mm <- paste0("0", mm)
  }

  if (hh < 10) {
    hh <- paste0("0", hh)
  }

  hhmm <- paste(hh, mm, sep = ":")
  hhmm
}
