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

# calculate flow metrics
metrics <- calculate_metrics(flow_futures, recompile = FALSE)

# list all possible future scenarios 
#   (81 combinations per species and waterbody)
scenario_options <- flow_futures |>
  distinct(species, waterbody, future, scenario) |>
  mutate(future_next = future, scenario_next = scenario) |>
  complete(nesting(species, waterbody), future, future_next, scenario, scenario_next)

# for each future, grab the correct 2024 metrics and then grab the
#   2025 metrics based on the *_next settings, but replace antecedent
#   for that scenario with the appropriate 2024 value
metrics_2024 <- metrics |>
  filter(water_year == max(water_year))
metrics_2025 <- metrics_2024 |> mutate(water_year = water_year + 1)
metrics_future <- scenario_options |>
  left_join(
    metrics_2024,
    by = c("species", "waterbody", "future", "scenario")
  ) |>
  left_join(
    metrics_2025,
    by = c("species", "waterbody", "future_next" = "future", "scenario_next" = "scenario"),
    suffix = c("_2024", "_2025")
  ) |>
  mutate(
    proportional_antecedent_flow_2025 = proportional_annual_flow_2024,
    proportional_max_antecedent_2025 = proportional_max_annual_2024
  )
metrics_2024 <- metrics_future |>
  select(species, waterbody, future, future_next, scenario, scenario_next, contains("2024")) |>
  rename(water_year = water_year_2024) |>
  rename_with(\(x) gsub("_2024", "", x), contains("2024"))
metrics_2025 <- metrics_future |>
  select(species, waterbody, future, future_next, scenario, scenario_next, contains("2025")) |>
  rename_with(\(x) gsub("_2025", "", x), contains("2025"))
metrics_future <- bind_rows(metrics_2024, metrics_2025)

# load data wihtout e-flows (check reaches and filter to those?)
# Ovens, Upper Loddon, and Broken R Reach 3 have no e-water (very small amounts)
counterfactual_lu <- c(
  "lowerbrokenc_r4" = "broken_creek_r4",
  "campaspe_r4" = "campaspe_river_r4",
  "glenelg_r1b" = "glenelg_river_r1",
  "glenelg_r2" = "glenelg_river_r2",
  "glenelg_r3" = "glenelg_river_r3",
  "goulburn_r4" = "goulburn_river_r4",
  "loddon_r4" = "loddon_river_r4",
  "macalister_r1" = "macalister_river_r1",
  "mackenzie_r3" = "mackenzie_river_r3",
  "moorabool_r3b" = "moorabool_river_r3",
  "thomson_r3" = "thomson_river_r3"
)

# load counterfactual data and filter to target systems
counterfactual <- read.csv("data/vic-ewater-collated-long.csv")

# reformat counterfactuals to match flow_futures and remove negative
#   counterfactual flow values
counterfactual <- counterfactual |>
  mutate(
    date_formatted = parse_date_time(date_formatted, orders = c("ymd_HMS", "ymd")),
    waterbody = counterfactual_lu[system],
    scenario = ifelse(scenario == "No e-water", "counterfactual", scenario),
    water_temperature_c = NA,
    future = NA,
    stream_discharge_mld = discharge_mld,
    stream_discharge_mld = ifelse(stream_discharge_mld < 0, 0, stream_discharge_mld)
  ) |>
  filter(scenario == "counterfactual", waterbody %in% counterfactual_lu) |>
  as_tibble() |>
  left_join(
    fetch_hypoxia_risk(),
    by = "waterbody"
  ) |>
  left_join(
    fetch_carrying_capacity(), 
    by = "waterbody",
    relationship = "many-to-many"
  ) |>
  select(all_of(colnames(flow_futures)))

# re-calculate flow metrics for counterfactual scenarios
metrics_counterfactual <- calculate_metrics(
  counterfactual,
  recompile = FALSE,
  suffix = "counterfactual"
)

# then left_join this to the pre-2024 metrics set and fill NAs 
#    with observed values (i.e., no counterfactual implies no e-water)
metrics_observed <- metrics |>
  filter(future == "ave", scenario == "none", water_year < 2024) |>
  select(-future, -scenario) |>
  left_join(
    metrics_counterfactual |> select(-future, -scenario),
    by = c("species", "waterbody", "water_year"),
    suffix = c("", "_counterfactual")
  ) |>
  mutate(
    across(
      contains("_counterfactual"),
      .fns = \(x) ifelse(is.na(x), get(gsub("_counterfactual", "", cur_column())), x)
    )
  )
metrics_counterfactual <- metrics_observed |>
  select(species, waterbody, water_year, contains("_counterfactual")) |>
  rename_with(\(x) gsub("_counterfactual", "", x), contains("_counterfactual"))
metrics_observed <- metrics_observed |>
  select(!contains("_counterfactual"))

# calculate percentage overhang to use as a predictor in blackfish model
# steps: list all ISC reaches in each VEWH reach
#        calculate area of streambed_width_r (total area of the river segment)
#        work out area of overhang clipped to the streambed width
#        divide overhang by total area to give proportion overhang

# initialise populations
specify_initial_conditions <- function(species, waterbody, cpue, start, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    !species %in% c("gadopsis_marmoratus", "maccullochella_peelii", "melanotaenia_fluviatilis")
  )
  
  # rename waterbody to avoid conflicts below
  wb <- waterbody
  
  # pull out CPUE for the target species
  cpue_sub <- cpue |> 
    mutate(
      sciname = tolower(gsub(" ", "_", scientific_name)),
      waterbody = paste0(
        tolower(gsub(" ", "_", waterbody)),
        reach_no,
        sep = "_r"
      )
    ) |>
    filter(
      sciname == species,
      waterbody == wb,
      survey_year == start
    )
  
  # calculate total CPUE per reach
  cpue_sub <- cpue_sub |>
    group_by(survey_year) |>
    summarise(
      catch = sum(catch),
      effort_h = sum(effort_h)
    ) |>
    mutate(
      cpue = catch / effort_h,
      cpue_max = max(cpue)
    )

  # pull out reach length for target waterbody
  reach_length <- .vewh_reach_lengths |> 
    filter(waterbody == wb) |>
    pull(reach_length) |>
    as.numeric()
  
  # extract carrying capacity for species and reach
  k <- flow_futures |> 
    rename(sciname = species) |>
    filter(sciname == species, waterbody == wb) |>
    pull(carrying_capacity) |>
    unique()
  
  # k is 20% of capacity for river blackfish (why?)
  if (species == "gadopis_marmoratus") 
    k <- 0.2 * k
  
  # work out the rescaling rate for catch to abundance based on a capture
  #    rate of 40%
  adult_rescale <- (1 / 0.4) * (reach_length / 100)
  init <- set_initial(
    species = species,
    cpue = cpue$cpue,
    cpue_max = cpue$cpue_max,
    effort_h = cpue$effort_h,
    n = k,
    nsim = nsim,
    rescale = adult_rescale
  )

  # return
  init
  
}


# set up all three population models (in a list?)
specify_pop_model <- function(species, waterbody, ntime, nstocked, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    !species %in% c("gadopsis_marmoratus", "maccullochella_peelii", "melanotaenia_fluviatilis")
  )

  # extract carrying capacity for species and reach
  k <- flow_futures |> 
    rename(sciname = species, wb = waterbody) |>
    filter(sciname == species, wb == waterbody) |>
    pull(carrying_capacity) |>
    unique()
  
  # river blackfish model
  if (species == "gadopsis_marmoratus") {
    mod <- river_blackfish(
      k = k,
      ntime = ntime
    )
  }
  
  # murray cod model
  if (species == "maccullochella_peelii") {
    
    # use a lookup to define system from waterbody
    system <- switch(
      waterbody,
      "broken_creek_r4" = "largetrib",
      "largetrib"
    )
    
    mod <- murray_cod(
      k = k,
      system = system,
      n = list(
        # number stocked, accounting for fingerling mortality and 50:50 sex ratio
        nstocked[seq_len(ntime)],
        rep(0, ntime),
        rep(0, ntime)
      ),
      ntime = ntime, 
      start = rep(1, 3), 
      end = rep(ntime, 3), 
      add = c(TRUE, TRUE, TRUE),
      p_capture = 0.1,   # 10% capture probability for any fish in slot
      slot = c(550, 750) # slot in mm
    )
  }
  
  # rainbowfish model
  if (species == "melanotaenia_fluviatilis") {
    mod <- murray_rainbowfish(
      k = k,
      ntime = ntime
    )
  }
  
  # return
  mod
  
}

simulate_scenario <- function(species, x, nsim, init, metrics, coefs, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    !species %in% c("gadopsis_marmoratus", "maccullochella_peelii", "melanotaenia_fluviatilis")
  )
  
  # river blackfish model
  if (species == "gadopsis_marmoratus") {
    sims_init <- simulate(
      x,
      nsim = nsim,
      init = init,
      args = list(
        covariates = c(
          format_covariates(metrics),
          list(
            coefs = coefs[1:5],
            temperature_coefficient = coefs[6]
          )
        )
      ),
      options = list(
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
  }
  
  # murray cod model
  if (species == "maccullochella_peelii") {
    sims_init <- simulate(
      x,
      nsim = nsim,
      init = init,
      args = list(
        covariates = c(
          format_covariates(metrics),
          list(threshold = 0.05),
          list(coefs = coefs)
        ),
        density_dependence = list(
          kdyn = lapply(
            seq_len(ntime), 
            function(.x) metrics$kdyn[.x]
          )
        )
      ),
      options = list(
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
  }
  
  # rainbowfish model
  # TODO: consider updating to age-based model
  if (species == "melanotaenia_fluviatilis") {
    sims <- simulate(
      x,
      nsim = nsim,
      init = init,
      args = list(
        covariates = c(
          format_covariates(metrics),
          list(coefs = coefs)
        )
      ),
      options = list(
        update = update_multinomial,
        tidy_abundances = floor
      )
    )
  }
  
  # return
  sims
  
}

# function to return coefficients for each species
get_coefs <- function(species) {
  coefs <- list(
    "gadopsis_marmoratus" = list(
      "glenelg_river_r1" = c(-5, 50, 50, 100, 75, 0.1),
      "glenelg_river_r2" = c(-5, 50, 50, 100, 50, 0.1),
      "glenelg_river_r3" = c(-5, 50, 50, 100, 50, 0.1),
      "loddon_river_r2" = c(-5, 50, 25, 70, 80, 0.1),
      "macalister_river_r1" = c(-5, 50, 15, 10, 20, 0.1),
      "mackenzie_river_r3" = c(-5, 50, 25, 70, 80, 0.1),
      "moorabool_river_r3" = c(-5, 50, 40, 40, 30, 0.1),
      "thomson_river_r3" = c(-5, 50, 20, 10, 15, 0.1)
    ),
    "maccullochella_peelii" = list(
      "broken_creek_r4" = c(0, -100, 6, -30, 100, 100),
      "broken_river_r3" = c(-15, 35, 45, -30, 80, 25),
      "campaspe_river_r4" = c(-60, 10, 20, -15, 45, 10),
      "goulburn_river_r4" = c(-5, 20, 6, -10, 30, 10),
      "loddon_river_r4" = c(-30, 20, 10, -10, 30, 10),
      "ovens_river_r5" = c(-10, 60, 20, -30, 60, 25)
    ),
    "melanotaenia_fluviatilis" = list(
      "broken_creek_r4" = c(-130, 40, 20),
      "broken_river_r3" = c(-180, 45, 30),
      "campaspe_river_r4" = c(-100, 10, 5),
      "goulburn_river_r4" = c(30, -10, 20),
      "loddon_river_r4" = c(-120, 10, 30),
      "ovens_river_r5" = c(-20, 40, 20)
    )
  )
  coefs[[species]]
}


# simulate for each species in turn
species_list <- metrics_observed |> pull(species) |> unique()
for (i in seq_along(species_list)) {
  
  # pull out metrics for relevant species
  get_metric_names <- function(species) {
    metric_list <- list(
      "gadopsis_marmoratus" = c("proportional_spring_flow"),
      "maccullochella_peelii" = c(""),
      "melanotaenia_fluviatilis" = c("")
    )
    metric_list[[species]]
  }
  metrics_observed_sp <- metrics_observed |>
    filter(species == all_of(species_list[i])) |>
    select(waterbody, water_year, all_of(get_metric_names(species_list[i])))
  metrics_counterfactual_sp <- metrics_counterfactual |>
    filter(species == all_of(species_list[i])) |>
    select(waterbody, water_year, all_of(get_metric_names(species_list[i])))
  
  # TODO: add extra metrics for rainbowfish (redfin, gambusia, instream cover)
  #   and for blackfish (overhanging veg cover)
  
  # simulate for each waterbody in turn
  waterbodies <- metrics_observed_sp |> pull(waterbody) |> unique()  
  for (j in seq_along(waterbodies)) {
    
    # specify initial conditions
    initial <- specify_initial_conditions(species_list[i], waterbodies[j], cpue)
    
    # filter to each waterbody in turn
    metrics_observed_wb <- metrics_observed_sp |>
      filter(waterbody == all_of(waterbodies[j]))
    metrics_counterfactual_wb <- metrics_counterfactual_sp |>
      filter(waterbody == all_of(waterbodies[j]))
    
    # grab stocking info if required
    if (species_list[i] == "maccullochella_peelii") {
      stocking_rates <- stocking |>
        filter(
          Species == "Murray Cod",
          System == system_lu[waterbodies[j]]
        ) |>
        select(Year, Number) |>
        mutate(Year = Year + 1) |>
        rename(water_year = Year, number_stocked = Number)
      metrics_observed_wb <- metrics_observed_wb |>
        left_join(stocking_rates, by = c("water_year"))
      metrics_counterfactual_wb <- metrics_counterfactual_wb |>
        left_join(stocking_rates, by = c("water_year"))
    }
    
    # initialise population model for this species and waterbody
    # TODO: make this work for all three species
    pop <- specify_pop_model(species_list[i], waterbodies[j], metrics)
    
    # simulate model
    sims_observed <- simulate_scenario(species_list[i], pop_models[[species_list[i]]])
    sims_counterfactual <- simulate_scenario(species_list[i], pop_models[[species_list[i]]])
    
    # save output
    qsave(
      sims_observed, 
      file = paste0("outputs/simulated/observed-", species_list[i], "-", waterbodies[j], ".qs")
    )
    qsave(
      sims_counterfactual, 
      file = paste0("outputs/simulated/counterfactual-", species_list[i], "-", waterbodies[j], ".qs")
    )
    
    # simulate futures
    future_sub <- metrics_future |>
      filter(
        species == species_list[i],
        waterbody == waterbodies[j]
      ) |>
      distinct(future, future_next, scenario, scenario_next)
    for (k in seq_len(nrow(future_sub))) {
      
      # extract initial conditions from sims_observed
      initial_future <- sims_observed[, , dim(sims_observed)[3]]
      
      # pull out metrics for a given scenario
      metrics_future_sub <- metrics_future |>
        filter(
          species == species_list[i],
          waterbody == waterbodies[j],
          future == future_sub$future[k],
          future_next == future_sub$future_next[k],
          scenario == future_sub$scenario[k],
          scenario_next == future_sub$scenario_next[k]
        )
      
      # simulate under the specific scenario
      sims_future <- simulate_scenario(species_list[i], pop_models[[species_list[i]]])
      
      # and save output
      future_name <- paste(
        future_sub$future[k],
        future_sub$future_next[k],
        future_sub$scenario[k],
        future_sub$scenario_next[k],
        sep = "_"
      )
      qsave(
        sims_future, 
        file = paste0("outputs/simulated/future-", species_list[i], "-", waterbodies[j], "-", future_name, ".qs")
      )
      
    }
    
  }
  
}

# summarise all outputs
# model validation (using observed- outputs)

# comparison of flows with and without e-water (using observed- and counterfactual- outputs)

# not considered: long-term scenarios under RCP rescaling

# futures (using future- outputs)

