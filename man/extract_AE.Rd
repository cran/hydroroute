% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AE.R
\name{extract_AE}
\alias{extract_AE}
\title{Extract Associated Events}
\usage{
extract_AE(
  relation_path,
  events_path,
  settings_path,
  unique = c("time", "metric"),
  inputdec = ".",
  inputsep = ",",
  saveResults = FALSE,
  outdir = tempdir(),
  TimeFormat = "\%Y-\%m-\%d \%H:\%M",
  tz = "Etc/GMT-1"
)
}
\arguments{
\item{relation_path}{Character string containing the path of the file where
the relation file is to be read from with
\code{\link[utils:read.csv]{utils::read.csv()}}. The file must contain a
column \code{ID} that contains the gauging station ID's in the file have to be
in order of their location in downstream direction.}

\item{events_path}{Character string containing the path of the directory
where the event files corresponding to the `relation` file are located. Only
relevant files in this directory will be used, i.e., files that are related to
the `relation` file.}

\item{settings_path}{Character string containing the path of the file where
the settings file is to be read from with
\code{\link[utils:read.csv]{utils::read.csv()}}. The file must be in the format
of the output of \code{\link[=peaktrace]{peaktrace()}}.}

\item{unique}{Character string specifying if the potential AEs which
meet the \code{timeLag} and \code{metricLag} condition should
be filtered to contain only unique events using \code{"time"},
i.e., by selecting those where the time difference is smallest
compared to the specified factor of the mean translation time, or using
\code{"metric"}, i.e., by selecting those where the relative
difference in amplitude is smallest (default: \code{"time"}).}

\item{inputdec}{Character string for decimal points in input data.}

\item{inputsep}{Field separator character string for input data.}

\item{saveResults}{A logical. If \code{FALSE} (default), the extracted AEs are not saved.
Otherwise the extracted AEs are written to a csv file.}

\item{outdir}{Character string naming a directory where the extraced AEs
should be saved to.}

\item{TimeFormat}{Character string giving the date-time format of the
date-time column in the input data frame (default: "\%Y-\%m-\%d \%H:\%M").}

\item{tz}{Character string specifying the time zone to be used for the
conversion (default: "Etc/GMT-1").}
}
\value{
A data frame containing \dQuote{real} AEs (i.e., events
    where the time differences and the relative difference in
    amplitude is within the limits and cut points provided by the
    file in \code{settings_path}). If no AEs can be found between the first
    two neighboring stations, \code{NULL} is returned. Otherwise the function
    returns all \dQuote{real} AEs that could be found along the river section
    specified in the file from \code{relation_path}. A warning is issued when
    the extraction is stopped early and shows the \code{IDs} for which no
    AEs are determined.
}
\description{
For given relation and event data return the
    associated events which comply with the conditions specified in
    the settings.
}
\examples{
relation_path <- system.file("testdata", "relation.csv", package = "hydroroute")
events_path <- system.file("testdata", "Events", package = "hydroroute")
settings_path <- system.file("testdata", "Q_event_2_AMP-LAG_aut_settings.csv",
                                   package = "hydroroute")
real_AE <- extract_AE(relation_path, events_path, settings_path)
}
