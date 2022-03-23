#' Fitting Models
#' @description Fit the specified models for all metrics available in the data set
#'     provided. 
#' @param formula An object of class \code{\link[stats:formula]{stats::formula()}}
#'     to fit models.
#' @param data A data frame in a format as available in the first list element from
#'     \code{\link[=potential_AE]{potential_AE()}}.
#' @param model A function which specifies the method used for fitting models
#'     (default: \code{\link[stats:lm]{stats::lm()}}). The model class must have
#'     a \code{\link[stats:predict]{stats::predict()}} function.
#' @param ... Additional arguments to be passed to the function specified in
#'     argument \code{model}.
#'
#' @return A list where the first element contains a list of model
#'     objects and the second element contains the models in a data
#'     frame.
#' @keywords internal
#' @noRd
fit_model <- function(formula = y ~ x, data, model = stats::lm, ...) {
  stopifnot(all(c("station.x", "station.y", "metric") %in% colnames(data)))

  if (is.factor(data$metric)) {
    metrics <- levels(data$metric)
  } else {
    metrics <- sort(unique(data$metric))
  }

  data <- stats::na.omit(data)
  r2 <- vector(mode = "numeric", length = length(metrics))
  models <- vector(mode = "list", length(metrics))
  n <- vector(mode = "numeric", length(metrics))

  for (i in seq_along(metrics)) {
    data_subset <- data[data$metric == metrics[i], ]
    n[i] <- nrow(data_subset)

    if (nrow(data_subset) > 0) {
      models[[i]] <- model(formula = formula, data = data_subset, ...)
      r2[i] <- get_r2(models[[i]], data_subset)
    }
  }
  basicinfo <- unique(data[, c("station.x", "station.y")])
  rownames(basicinfo) <- NULL
  stopifnot(nrow(basicinfo) == 1)

  type <- unname(do.call(rbind, lapply(models, class)))
  type[type == "NULL"] <- NA_character_

  formula <- as.character(lapply(models, formula))
  formula[formula == "list()"] <- NA_character_

  coefs <- lapply(models, stats::coef)
  # replace NULL with NA_real_
  coefs[sapply(coefs, is.null)] <- NA_real_
  coefs <- do.call("rbind", coefs)

  models_df <- cbind(basicinfo, metric = metrics, type, formula, coefs, n, r2)
  rownames(models_df) <- NULL

  models <- append(c(unique(data$station.x), unique(data$station.y)), models)

  dim(models) <- c(1, length(models))
  colnames(models) <- c("station.x", "station.y", metrics)

  list(models = models, models.df = models_df)
}

#' R-Squared
#' @description Extract the R-squared from a model (\code{$r.squared}) or calculate it.
#' @param model Object of a model class that returns an
#'     \code{$r.squared} object or has a \code{summary} and a
#'     \code{predict} function to compute the R-squared value.
#' @param data Data frame that was used to fit the \code{model}. Can be omitted,
#'     if the \code{model} returns an \code{$r.squared} object.
#' @return A numeric value giving the R-squared of the model.
#' @keywords internal
#' @noRd
get_r2 <- function(model, data) {
  if (!is.null(summary(model)$r.squared)) {
    r2 <- summary(model)$r.squared
  } else {
    yi_hat <- stats::predict(model)
    response <- all.vars(stats::formula(model))[1]
    yi <- data[, response]
    y_mean <- mean(yi)

    r2 <- sum((yi_hat - y_mean)^2) / sum((yi - y_mean)^2)
  }
  return(r2)
}

#' Model Predictions
#' @description Given fitted regression models obtain the predicted
#'     values given some initial values.
#' @param models A list matrix containing in each row the fitted
#'     models for a pair of stations with the columns containing the
#'     information on the involved stations as well as the fitted
#'     models for the different metrics. Such an object can be
#'     obtained as the first list element of the object returned by
#'     \code{\link[=fit_model]{fit_model()}}.
#' @param initials A data frame which contains initial values for
#'     predictions (see vignette).
#' @return A data frame containing the predicted values.
#' @keywords internal
#' @noRd
get_predictions <- function(models, initials) {
  # no models
  stopifnot(!is.null(models))

  # missing values
  stopifnot(!is.na(initials))

  # remove empty rows
  for (i in seq_len(nrow(models))) {
    if (is.null(models[i, ])) {
      models <- models[-i, ]
    }
  }

  colnames(initials) <- tolower(colnames(initials))
  stopifnot(!any(is.na(initials)))

  if (!is.factor(initials$metric)) {
    initials$metric <- factor(initials$metric,
                              levels = c("AMP", "MAFR", "MEFR", "DUR", "RATIO"))
  }

  metrics <- levels(initials$metric)

  results <- vector(mode = "list", length = length(metrics))

  # metric
  for (i in seq_along(metrics)) {
    initials_subset <- subset(initials, subset = metric == metrics[i])

    # name / scenario
    results_name <- vector(mode = "list", length = length(initials_subset$name))
    for (j in seq_along(initials_subset$name)) {

      # for each scenario: get predictions from each model, plus initial values
      results_model <- vector(mode = "list", length = (length(models[,1]) + 1))
      results_model[[1]] <- data.frame(station = initials_subset$station[j],
                                       metric = initials_subset$metric[j],
                                       prediction = initials_subset$value[j],
                                       name = initials_subset$name[j])
      k <- 1        
      while (k <= nrow(models)) {
        if (models[[k, 2]] != initials_subset$station[j]) {
          prev <- unname(subset(initials_subset, subset = name == initials_subset$name[j],
                                select = "value"))
          prev <- stats::predict(models[[k, (i + 2)]], newdata = data.frame(x = prev))
          results_model[[k + 1]] <- data.frame(station = models[[k, 2]], # station.y
                                               metric = initials_subset$metric[j],
                                               prediction = prev,
                                               name = initials_subset$name[j])
        }
        k <- k + 1
      }
      results_name[[j]] <- do.call("rbind", results_model)
    }

    results[[i]] <- do.call("rbind", results_name)
  }

  results <- do.call("rbind", results)
  rownames(results) <- NULL
  results
}

#' Impute Missing Initial Values
#' @description Determine initial values in a data-driven way if missing. 
#' @param real_AE Data frame containing identified \dQuote{real} AEs. Output from
#'     \code{\link[=estimate_AE]{estimate_AE()}}.
#' @param initials A data frame which contains initial values for predictions
#'     (see vignette).
#' @param method A function that is used to impute missing values based on data
#'     \code{real_AE}.
#' @return A data frame that corresponds to the \code{initials} input data frame
#'     with no missing values.
#' @keywords internal
#' @noRd
impute_initials <- function(real_AE, initials, method = base::max) {
  colnames(initials) <- tolower(colnames(initials))
  colnames(real_AE) <- tolower(colnames(real_AE))
  initials$metric <- tolower(initials$metric)

  real_AE_long <- real_AE[, !(names(real_AE) %in% "diff_metric")]
  real_AE_long <- reshape(real_AE_long)

  for (i in which(is.na(initials$value))) {
    STATION <- as.character(initials[i, "station"])
    METRIC <- as.character(initials[i, "metric"])

    # impute from Sx or Sy
    if (unique(real_AE_long$station.x == STATION)) {
      initials[i, "value"] <- method(subset(real_AE_long, metric == METRIC)$x)

    } else if (unique(real_AE_long$station.y == STATION)) {
      initials[i, "value"] <- method(subset(real_AE_long, metric == METRIC)$y)
    }
  }
  initials$metric <- toupper(initials$metric)

  initials
}
