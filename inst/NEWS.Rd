\name{NEWS}
\title{News for Package 'hydroroute'}
\section{Changes in hydroroute version 0.1.1}{
  \itemize{
    \item Ensure that `get_lag_dir` uses `mc.cores = 1` when calling
    `get_lag_file` in the `mclapply` call.
    \item Argument `unique` added which allows to specify if potential
    AEs should be filtered to contain only unique events by selecting
    associated events based on the best time lag match or the best
    metric match.
    \item Implementation of the automatic procedure corrected to
    corrspond to the method described in Greimel et al. (2022).
    \item Identification of potential AEs corrected to follow events
    across all gauging stations in the relations file and ensure that
    events have unique matches.
    \item Minor changes in documentation.
    \item Minor bug fixes to handle states where no potential AEs or matches are
    found.
    \item First version released on CRAN: 0.1.0
  }
}
