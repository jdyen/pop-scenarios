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
# Date modified: 11 October 2023

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
source("R/fish.R")
source("R/flow.R")

# specify species and locations
sciname_list <- c(
  "Gadopsis marmoratus",
  "Maccullochella peelii",
  "Melanotaenia fluviatilis"
)
waterbody_list <- list(
  gadopsis_marmoratus = c(
    "Glenelg River",
    "Loddon River",
    "Mackenzie River",
    "Macalister River",
    "Moorabool River",
    "Thomson River"
  ),
  maccullochella_peelii = c(
    "Broken Creek",
    "Broken River",
    "Campaspe River",
    "Goulburn River", 
    "Loddon River",
    "Ovens River"
  ),
  melanotaenia_fluviatilis = c(
    "Broken Creek",
    "Broken River",
    "Campaspe River",
    "Goulburn River",
    "Loddon River",
    "Ovens River"
  )
)
reach_list <- c(
  "Broken Creek_r4",
  "Broken River_r3",  # use 404224, very little e-water delivered here
  "Campaspe River_r4",
  "Glenelg River_r1", # use R1b gauge flows
  "Glenelg River_r2",
  "Glenelg River_r3",
  "Goulburn River_r4",
  "Loddon River_r2",  # no compliance assessed here, so no e-water recorded
  "Loddon River_r4",
  "Macalister River_r1",
  "Mackenzie River_r3",  # @ McKenzie Ck, too far down but no gauge in R1
  "Moorabool River_r3",
  "Ovens River_r5",  # very little e-water here
  "Thomson River_r3"
)

carrying_capacity <- list(
  gadopsis_marmoratus = c(
    "glenelg_river_r1" = 50000,
    "glenelg_river_r2" = 100000,
    "glenelg_river_r3" = 200000,
    "loddon_river_r2" = 100000,
    "macalister_river_r1" = 100000,
    "mackenzie_river_r3" = 50000,
    "moorabool_river_r3" = 100000,
    "thomson_river_r3" = 200000
  ),
  maccullochella_peelii = c(
    "broken_creek_r4" = 50000,
    "broken_river_r3" = 200000,
    "campaspe_river_r4" = 50000,
    "goulburn_river_r4" = 100000,
    "loddon_river_r4" = 200000,
    "ovens_river_r5" = 100000
  ),
  melanotaenia_fluviatilis = c(
    "broken_creek_r4" = 1000,
    "broken_river_r3" = 5000,
    "campaspe_river_r4" = 5000,
    "goulburn_river_r4" = 10000,
    "loddon_river_r4" = 3000,
    "ovens_river_r5" = 10000
  )
)

# load data
cpue <- fetch_fish(recompile = FALSE)

# load stocking info
stocking <- read.csv("data/stocking.csv")

# fetch hydrological data
flow <- fetch_flow(start = 2009, end = 2023, recompile = FALSE)

# specify futures (individual years and events; can chop and change these to 
#    create specific scenarios)
flow_futures <- specify_flow_future(flow)

# TODO: do chopping and changing here...
# TODO: see if metrics can be calculated for the futures only (which would
#   make it easier to append these to the main flow scenarios)

# specify climate change scenarios
# flow_gcm <- apply_gcm(flow)  # add a lookup for catchments, need to add a few still

# calculate flow metrics
calculate_metrics <- function(x, xname, k, species) {

  # specify hypoxia risk by system
  hypoxia_risk <- NULL
  if (xname == "broken_creek_r4")
    hypoxia_risk <- c(350, 25)

  # calculate and return k-scaled flow metrics  
  calculate_flow_metrics(
    x, 
    carrying_capacity = k,
    hypoxia_risk = hypoxia_risk
  )
  
}
metrics <- vector("list", length = length(carrying_capacity))
names(metrics) <- names(carrying_capacity)
for (i in seq_along(carrying_capacity)) {
  flow_sub <- flow[names(carrying_capacity[[i]])]
  metrics[[i]] <- mapply(
    calculate_metrics, 
    flow_sub,
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

