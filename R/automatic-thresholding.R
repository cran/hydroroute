#' Fit a Parabola to a Histogram
#' @description Using the relative frequency values of a equi-spaced
#'     histogram fit parabolas where the mode corresponds to the
#'     vertex considering different range lengths left and right to
#'     the mode.
#' @param x Contingency table (a `\code{table}' object).
#' @param omit Number of boundary points not to consider when
#'     determining the mode.
#' @return A list which contains the estimated coefficients and the
#'     RMSE value of the parabola with the minimal RMSE value.
#' @keywords Internal
#' @noRd
get_parabola <- function(x, omit = 2) {
    x <- unname(x)

    d <- diff(x)
    is <- which(d[-length(d)] > 0 & d[-1] < 0) + 1
    is <- is[is > omit & is < (length(x) - omit)]

    if (length(is) == 0) {
        result <- list()
    } else {
        istar <- is[which.max(x[is])]

        error <- vapply(seq_len(min(istar, length(x) - istar) - 1), function(i) {
            index <- -i:i
            y <- x[istar + index]
            b <- y[i + 1]
            fit <- try(stats::nls(y ~ a * index^2 + b, start = list(a = 0)), silent = TRUE)
            if (is(fit, "try-error")) {
                Inf
            } else {
                coefs <- c(stats::coef(fit), b = b, c = istar)
                xstar <- suppressWarnings(coefs["c"] + c(-1, 1) * sqrt(-coefs["b"] / coefs["a"]))
                index <- which(seq_along(x) > xstar[1] & seq_along(x) < xstar[2]) - istar
                ystar <- stats::predict(fit, data.frame(index = index))
                mean((x[istar + index] - ystar)^2)
                }
            }, FUN.VALUE = numeric(1))
        i <- which.min(error)
        if (length(i) == 0 | is.infinite(error[i])) {
            result <- list()
        } else {
            index <- -i:i
            y <- x[istar + index]
            b <- y[i + 1]
            fit <- stats::nls(y ~ a * index^2 + b, start = list(a = 0))
            coefs <- c(stats::coef(fit), b = b, c = istar)
            xstar <- suppressWarnings(coefs["c"] + c(-1, 1) * sqrt(-coefs["b"] / coefs["a"]))

            summary_fit <- summary(fit)

            result <- list(coefs = coefs, xstar = xstar, rmse = summary_fit$sigma)
        }
    }
    return(result)
}
