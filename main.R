# Analysis of flow scenarios for three fish species:
#    Murray cod (Maccullochella peelii)
#    River blackfish (Gadopsis marmoratus)
#    Murray-Darling rainbowfish (Melanotaenia fluviatilis)
# 
# Flow scenarios consider observed daily flows with and without
#   environmental water allocations from 2000-2023 and projected
#   daily flows under a scenario of climate change rescaled to 
#   match catchment-specific predicted runoff under the RCP 8.5 
#   emissions scenario projected to 2040.
#
# Author: Jian Yen (jdl.yen [at] gmail.com)
# 
# Date created: 5 July 2023
# Date modified: 12 October 2023

# TODO: update BF template with recruitment flow effects
#       update all templates with info from trends analysis (not many effects on adults, so be careful with this)

# load some packages
library(qs)
library(dplyr)
library(tidyr)
library(lubridate)
library(aae.db)
library(aae.hydro)
library(aae.pop.templates)
library(sf)
library(ggplot2)
library(ragg)

# and some load helpers
source("R/utils.R")
source("R/fish.R")
source("R/flow.R")

# load data
cpue <- fetch_fish(recompile = FALSE)

# load stocking info
stocking <- read.csv("data/stocking.csv")

# fetch hydrological data
flow <- fetch_flow(start = 2009, end = 2023, recompile = FALSE)

# specify futures (individual years and events; can chop and change these to 
#    create specific scenarios)
flow_futures <- specify_flow_future(flow)

# and add some extra information on hypoxia risk and K by species (expands
#   over species)
flow_futures <- flow_futures |>
  left_join(
    fetch_hypoxia_risk(),
    by = "waterbody"
  ) |>
  left_join(
    fetch_carrying_capacity(), 
    by = "waterbody",
    relationship = "many-to-many"
  )

# specify climate change scenarios
# flow_gcm <- apply_gcm(flow)  # add a lookup for catchments, need to add a few still

## TODO: redo metrics calculation fn to handle tibble rather htan lists
flow_futures_tmp |>
  group_by(waterbody, future, scenario) |>
  calculate(value = stream_discharge_mld, date = date_formatted, resolution = survey(9:11))

metrics <- vector("list", length = length(carrying_capacity))
names(metrics) <- names(carrying_capacity)
for (i in seq_along(carrying_capacity)) {
  flow_sub <- flow_futures[names(carrying_capacity[[i]])]
  metrics[[i]] <- vector("list", length = length(flow_sub))
  names(metrics[[i]]) <- names(flow_sub)
  for (j in seq_along(flow_sub)) {
    metrics[[i]][[j]] <- vector("list", length = length(flow_sub[[j]]))
    names(metrics[[i]][[j]]) <- names(flow_sub[[j]])
    for (k in seq_along(flow_sub[[j]])) {
      metrics[[i]][[j]][[k]] <- mapply(
        calculate_metrics, 
        flow_sub,
        names(flow_sub),
        carrying_capacity[[names(carrying_capacity)[i]]],
        MoreArgs = list(species = names(carrying_capacity)[i]),
        SIMPLIFY = FALSE
      )
    }
  }
}

metrics_futures <- vector("list", length = length(carrying_capacity))
names(metrics_futures) <- names(carrying_capacity)
for (i in seq_along(carrying_capacity)) {
  flow_sub <- flow_futures[names(carrying_capacity[[i]])]
  metrics_futures[[i]] <- mapply(
    calculate_metrics, 
    lapply(flow_sub, \(x) x$ave$none),
    names(flow_sub),
    carrying_capacity[[names(carrying_capacity)[i]]],
    MoreArgs = list(species = names(carrying_capacity)[i]),
    SIMPLIFY = FALSE
  )
}


## SEE IF CAN CALCULATE FOR FUTURES (NEED TO WORK OUT NESTED LISTS)

# load data wihtout e-flows (check reaches and filter to those?)
# Ovens, Upper Loddon, and Broken R Reach 3 have no e-water (very small amounts)
## NEED mackenzie R1/R2/R3/SOMETHING
counterfactual_lu <- c(
  "lowerbrokenc_r4" = "broken_creek_r4",
  "campaspe_r4" = "campaspe_river_r4",
  "glenelg_r1b" = "glenelg_river_r1",
  "glenelg_r2" = "glenelg_river_r2",
  "glenelg_r3" = "glenelg_river_r3",
  "goulburn_r4" = "goulburn_river_r4",
  "loddon_r4" = "loddon_river_r4",
  "macalister_r1" = "macalister_river_r1",
  "moorabool_r3b" = "moorabool_river_r3",
  "thomson_r3" = "thomson_river_r3"
)

## NEED TO USE OTHER METRICS ABOVE FOR OBSERVED, SET OBSERVED = COUNTERFACTUAL
##   FOR SYSTEMS WITH MINIMAL E-WATER

# load counterfactual data and calculate flow metrics
counterfactual <- read.csv("data/vic-ewater-collated-long.csv")
counterfactual <- lapply(
  counterfactual |> pull(system) |> unique(),
  \(x, y) y |> filter(system == x),
  y = counterfactual
)

# fetch hydrological data
flow_list <- fetch_flow(start = 1999, end = 2023, recompile = FALSE)

# calculate flow metrics
names(counterfactual) <- sapply(counterfactual, \(x) x$system[1])
counterfactual <- counterfactual[names(counterfactual) %in% names(counterfactual_lu)]
names(counterfactual) <- counterfactual_lu[names(counterfactual)]
counterfactual$wimmera_river_r3 <- counterfactual$wimmera_river_r2
flow_list <- flow_list[names(flow_list) %in% counterfactual_lu]
flow_list <- flow_list[match(names(counterfactual), names(flow_list))]
counterfactual <- mapply(
  \(x, y) x |>
    mutate(
      date_formatted = parse_date_time(date_formatted, orders = c("ymd_HMS", "ymd"))
    ) |> 
    left_join(
      y,
      by = "date_formatted"
    ) |>
    select(date_formatted, system, scenario, discharge_mld, stream_discharge_mld, water_temperature_c) |>
    mutate(
      scenario = ifelse(scenario == "No e-water", "counterfactual", scenario),
      stream_discharge_mld = ifelse(scenario == "counterfactual", discharge_mld, stream_discharge_mld)
    ) |>
    pivot_wider(
      id_cols = c(date_formatted, system, water_temperature_c),
      names_from = scenario,
      values_from = c(stream_discharge_mld)
    ) |>
    mutate(counterfactual = ifelse(is.na(counterfactual), Observed, counterfactual)),
  x = counterfactual,
  y = flow_list,
  SIMPLIFY = FALSE
)
metrics_obs <- mapply(
  calculate_metrics,
  lapply(counterfactual, \(x) x |> rename(stream_discharge_mld = Observed)),
  names(counterfactual), 
  SIMPLIFY = FALSE
)
metrics_cf <- mapply(
  calculate_metrics, 
  lapply(counterfactual, \(x) x |> rename(stream_discharge_mld = counterfactual)), 
  names(counterfactual), 
  SIMPLIFY = FALSE
)
metrics_obs <- bind_rows(metrics_obs)
metrics_cf <- bind_rows(metrics_cf)


# calculate percentage overhang to use as a predictor in blackfish model
# steps: list all ISC reaches in each VEWH reach
#        calculate area of streambed_width_r (total area of the river segment)
#        work out area of overhang clipped to the streambed width
#        divide overhang by total area to give proportion overhang

# add scenarios for futures (2 years out?)

# calculate flow metrics and other predictors

# initialise populations

# simulate

# outputs:
#  0. Validation steps? Generic plots? (Is the model any good? Use hindcasts and forecasts for 2023 against data, plus provide forecasts for 2024)
#  1. scenarios with and without e-flows (Comparison of scenarios)
#  1a. (optional) add long-term climate change scenarios (Comparison of scenarios, use rescaling for convenience)
#  2. near-term forecasts with flow options and climates
#  3. map the climate outputs and effects of e-flows?

