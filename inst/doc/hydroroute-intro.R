## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = 'center'
)
old <- options(useFancyQuotes = FALSE)
library("hydroroute")

## -----------------------------------------------------------------------------
Q_path <- system.file("testdata", "Q.csv", package = "hydroroute")
Q <- read.csv(Q_path)
head(Q)

## -----------------------------------------------------------------------------
relation_path <- system.file("testdata", "relation.csv", package = "hydroroute")
relation <- read.csv(relation_path)
relation

## -----------------------------------------------------------------------------
Sx <- system.file("testdata", "Events", "100000_2_2014-01-01_2014-02-28.csv",
                  package = "hydroroute")
Sx <- read.csv(Sx)
head(Sx)

## -----------------------------------------------------------------------------
Q_file <- system.file("testdata", "Q.csv", package = "hydroroute")
relation <- system.file("testdata", "relation.csv", package = "hydroroute")
get_lag_file(Q_file, relation, inputsep = ",", format = "%Y-%m-%d %H:%M",
             save = FALSE, overwrite = TRUE)

## -----------------------------------------------------------------------------
Q_file <- system.file("testdata", "Q.csv", package = "hydroroute")
relations_path <- system.file("testdata", package = "hydroroute")
get_lag_dir(Q_file, relations_path, inputsep = ",",
            inputdec = ".", format = "%Y-%m-%d %H:%M", overwrite = TRUE)

## -----------------------------------------------------------------------------
# file paths
Sx <- system.file("testdata", "Events", "100000_2_2014-01-01_2014-02-28.csv",
                  package = "hydroroute")
Sy <- system.file("testdata", "Events", "200000_2_2014-01-01_2014-02-28.csv",
                  package = "hydroroute")
relation <- system.file("testdata", "relation.csv", package = "hydroroute")

# read data
Sx <- utils::read.csv(Sx)
Sy <- utils::read.csv(Sy)
relation <- utils::read.csv(relation)
relation <- relation[1:2, ]

# estimate AE, exact time matches
results <- estimate_AE(Sx, Sy, relation, timeLag = c(0, 1, 0))

## -----------------------------------------------------------------------------
results$settings

## ---- fig.width = 6, fig.height = 3-------------------------------------------
results$plot_threshold

## -----------------------------------------------------------------------------
head(results$real_AE)

## -----------------------------------------------------------------------------
initial_values_path <- system.file("testdata", "initial_value_routing.csv",
                                   package = "hydroroute")
initials <- read.csv(initial_values_path)
initials

## -----------------------------------------------------------------------------
relation_path <- system.file("testdata", "relation.csv", package = "hydroroute")
events_path <- system.file("testdata", "Events", package = "hydroroute")
initial_values_path <- system.file("testdata", "initial_value_routing.csv",
                                   package = "hydroroute")
res <- peaktrace(relation_path, events_path, initial_values_path)

## -----------------------------------------------------------------------------
res$`2`$settings

## ---- fig.width = 8, fig.height = 3-------------------------------------------
grid::grid.draw(res$`2`$plot_threshold)

## -----------------------------------------------------------------------------
head(res$`2`$real_AE)

## ---- fig.width = 8, fig.height = 15------------------------------------------
grid::grid.draw(res$`2`$plot_scatter)

## -----------------------------------------------------------------------------
res$`2`$models

## ---- fig.width = 8, fig.height = 10------------------------------------------
gridExtra::grid.arrange(grobs = res$`2`$plot_predict$grobs, nrow = 3, ncol = 2)

## -----------------------------------------------------------------------------
initials[1, ]
res$`2`$models[1, ]

## -----------------------------------------------------------------------------
res$`2`$models[6, ]

## -----------------------------------------------------------------------------
res$`2`$models[11, ]

## -----------------------------------------------------------------------------
res$`4`$settings

## -----------------------------------------------------------------------------
head(res$`4`$real_AE)

## ---- fig.width = 8, fig.height = 3-------------------------------------------
grid::grid.draw(res$`4`$plot_threshold)

## ---- fig.width = 8, fig.height = 15------------------------------------------
grid::grid.draw(res$`4`$plot_scatter)

## -----------------------------------------------------------------------------
res$`4`$models

## ---- fig.width = 8, fig.height = 10------------------------------------------
gridExtra::grid.arrange(grobs = res$`4`$plot_predict$grobs, nrow = 3, ncol = 2)

## -----------------------------------------------------------------------------
relation_path <- system.file("testdata", "relation.csv", package = "hydroroute")
events_path <- system.file("testdata", "Events", package = "hydroroute")
settings_path <- system.file("testdata", "Q_event_2_AMP-LAG_settings.csv",
	                     package = "hydroroute")
initials_path <- system.file("testdata", "initial_value_routing.csv",
                             package = "hydroroute")
real_AE <- extract_AE(relation_path, events_path, settings_path)
head(real_AE)

## -----------------------------------------------------------------------------
relation <- utils::read.csv(relation_path)
initials <- utils::read.csv(initials_path)
res <- routing(real_AE, initials, relation)

## ---- fig.width = 8, fig.height = 15------------------------------------------
grid::grid.draw(res$plot_scatter)

## -----------------------------------------------------------------------------
res$models

## ---- fig.width = 8, fig.height = 10------------------------------------------
gridExtra::grid.arrange(grobs = res$plot_predict$grobs, nrow = 3, ncol = 2)

## ---- include = FALSE---------------------------------------------------------
options(old)

