# Analysis of flow scenarios for three fish species:
#    Murray cod (Maccullochella peelii)
#    River blackfish (Gadopsis marmoratus)
#    Murray-Darling rainbowfish (Melanotaenia fluviatilis)
# 
# Flow scenarios consider observed daily flows with and without
#   environmental water allocations from 2009-2023 and near-term (to 2025) 
#   forecasts of flows under different plausible climates and
#   flow management strategies
#
# Author: Jian Yen (jdl.yen [at] gmail.com)
# 
# Date created: 5 July 2023
# Date modified: 4 January 2024

# TODO: 
#   - update RBF and BF templates to increase average pop growth rates (flat-lining in most cases)
#     - BF done (coldwater impacts too strong)
#     - TODO: RBF
#   - update BF template with recruitment flow effects (??) and with overhang info
#      - template doens't need changing? Check
#   - update RBF model to be age-based
#   - check new metrics (voerhang and iwh) to make sure they're incorporated in all templates

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
library(rstanarm)
library(bayesplot)

# and some load helpers
source("R/utils.R")
source("R/fish.R")
source("R/flow.R")

# settings
set.seed(2023-10-17)
nsim <- 200
nburnin <- 10
simulate_again <- TRUE

# load data
cpue <- fetch_fish(recompile = FALSE)
cpue_exotic <- cpue |> 
  filter(grepl("perca|gambusia", scientific_name, ignore.case = TRUE))
cpue_exotic <- cpue_exotic |>
  mutate(waterbody = paste(
    tolower(gsub(" ", "_", waterbody)), reach_no, sep = "_r")
  ) |>
  group_by(scientific_name, waterbody, survey_year) |>
  summarise(presence = ifelse(sum(cpue) > 0, 1, 0))
cpue <- cpue |> 
  filter(!grepl("perca|gambusia", scientific_name, ignore.case = TRUE))

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

# add extra metrics for rainbowfish (redfin, gambusia, instream_cover)
#   and for blackfish (overhanging veg cover, instream_cover)
veg_overhang <- read.csv("data/eflow_veg_overhang_231113.csv")
veg_overhang <- veg_overhang |>
  as_tibble() |>
  mutate(
    waterbody = gsub(" ", "_", tolower(waterbody)),
    waterbody = paste(waterbody, vewh_reach, sep = "_r")
  ) |>
  filter(waterbody %in% unique(metrics_observed$waterbody)) |>
  mutate(proportion_overhang = overhang_area_m2 / buff5m_area_m2)
instream_cover <- read.csv("data/vefmap_habitat_site_metrics_231208.csv")
instream_cover <- instream_cover |>
  left_join(
    fetch_site_info(instream_cover) |> 
      select(id_site, reach_no) |> 
      collect(),
    by = "id_site"
  )
instream_cover <- instream_cover |>
  group_by(waterbody, reach_no) |>
  summarise(
    iwh = median(cur_pred_m3_m2, na.rm = TRUE),
    proportion_overhang = median(mean_vo_percentage, na.rm = TRUE) / 100
  ) |>
  ungroup() |>
  mutate(
    waterbody = gsub(" ", "_", tolower(waterbody)),
    waterbody = paste(waterbody, reach_no, sep = "_r")
  )
veg_overhang <- veg_overhang |>
  select(waterbody, proportion_overhang) |>
  full_join(
    instream_cover |> select(-reach_no),
    by = "waterbody",
    suffix = c("", "_hab_metrics")
  ) |>
  mutate(
    proportion_overhang_hab_metrics = ifelse(
      is.na(proportion_overhang_hab_metrics),
      1.2 * proportion_overhang, 
      proportion_overhang_hab_metrics
    ),
    iwh = ifelse(waterbody == "loddon_river_r2", 0.0199, iwh),
    iwh = ifelse(waterbody == "macalister_river_r1", 0.00716, iwh),
    iwh = ifelse(waterbody == "mackenzie_river_r3", 0.00886, iwh),
    iwh = ifelse(waterbody == "ovens_river_r5", max(iwh, na.rm = TRUE), iwh),
    iwh = ifelse(waterbody == "moorabool_river_r0", 0.0145, iwh),
  ) |>
  bind_rows(
    tibble(
      waterbody = c("broken_creek_r4", "loddon_river_r4"),
      proportion_overhang = rep(NA, 2),
      iwh = c(0.0128, 0.00471),
      proportion_overhang_hab_metrics = c(0.797, 0.641)
    )
  ) |>
  mutate(
    iwh = iwh / max(iwh),
    proportion_overhang_hab_metrics = proportion_overhang_hab_metrics / 
      max(proportion_overhang_hab_metrics)
  )
metrics_observed <- metrics_observed |>
  left_join(
    cpue_exotic |>
      mutate(
        species = "redfin",
        species = ifelse(
          grepl("gambusia", scientific_name, ignore.case = TRUE),
          "gambusia", 
          species
        )
      ) |>
      pivot_wider(
        id_cols = c(waterbody, survey_year),
        names_from = species,
        values_from = presence,
        names_prefix = "presence_"
      ),
    by = c("waterbody", "water_year" = "survey_year")
  ) |>
  mutate(
    presence_redfin = ifelse(is.na(presence_redfin), 0, presence_redfin),
    presence_gambusia = ifelse(is.na(presence_gambusia), 0, presence_gambusia)
  ) |>
  left_join(
    veg_overhang |>
      select(waterbody, iwh, proportion_overhang_hab_metrics) |>
      rename(
        instream_cover = iwh, 
        veg_overhang = proportion_overhang_hab_metrics
      ),
    by = "waterbody"
  )
metrics_counterfactual <- metrics_counterfactual |>
  left_join(
    metrics_observed |> 
      select(species, waterbody, water_year, contains("presence_"), instream_cover, veg_overhang),
    by = c("species", "waterbody", "water_year")
  )

# repeat for future metrics, but just use most recent observations to project
#   future occurrences
max_year <- cpue_exotic |>
  group_by(waterbody) |>
  summarise(max_year = max(survey_year))
metrics_future <- metrics_future |>
  left_join(
    cpue_exotic |>
      mutate(
        species = "redfin",
        species = ifelse(
          grepl("gambusia", scientific_name, ignore.case = TRUE),
          "gambusia", 
          species
        )
      ) |>
      pivot_wider(
        id_cols = c(waterbody, survey_year),
        names_from = species,
        values_from = presence,
        names_prefix = "presence_"
      ) |>
      left_join(max_year, by = "waterbody") |>
      filter(survey_year == max_year) |>
      select(contains("presence_"), waterbody),
    by = c("waterbody")
  ) |>
  mutate(
    presence_redfin = ifelse(is.na(presence_redfin), 0, presence_redfin),
    presence_gambusia = ifelse(is.na(presence_gambusia), 0, presence_gambusia)
  ) |>
  left_join(
    veg_overhang |>
      select(waterbody, iwh, proportion_overhang_hab_metrics) |>
      rename(
        instream_cover = iwh, 
        veg_overhang = proportion_overhang_hab_metrics
      ),
    by = "waterbody"
  )

# simulate for each species in turn (if required)
if (simulate_again) {
  
  species_list <- metrics_observed |> pull(species) |> unique()
  for (i in seq_along(species_list)) {
    
    # pull out flow/covariate metrics for species
    metrics_observed_sp <- metrics_observed |>
      filter(species == species_list[i]) |>
      select(waterbody, water_year, all_of(get_metric_names(species_list[i])))
    metrics_counterfactual_sp <- metrics_counterfactual |>
      filter(species == species_list[i]) |>
      select(waterbody, water_year, all_of(get_metric_names(species_list[i])))
    
    # rename a few metrics for some species
    if (species_list[i] == "gadopsis_marmoratus") {
      metrics_observed_sp <- metrics_observed_sp |>
        rename(antecedent_flow = proportional_antecedent_flow)
      metrics_counterfactual_sp <- metrics_counterfactual_sp |>
        rename(antecedent_flow = proportional_antecedent_flow)
    }
    if (species_list[i] == "maccullochella_peelii") {
      metrics_observed_sp <- metrics_observed_sp |>
        rename(blackwater_risk = hypoxia_risk)
      metrics_counterfactual_sp <- metrics_counterfactual_sp |>
        rename(blackwater_risk = hypoxia_risk)
    }
    if (species_list[i] == "melanotaenia_fluviatilis") {
      metrics_observed_sp <- metrics_observed_sp |>
        rename(redfin = presence_redfin, gambusia = presence_gambusia)
      metrics_counterfactual_sp <- metrics_counterfactual_sp |>
        rename(redfin = presence_redfin, gambusia = presence_gambusia)
    }
    
    # simulate for each waterbody in turn
    waterbodies <- metrics_observed_sp |> pull(waterbody) |> unique()  
    for (j in seq_along(waterbodies)) {
      
      # extract carrying capacity for species and reach
      k <- flow_futures |> 
        filter(species == species_list[i], waterbody == waterbodies[j]) |>
        pull(carrying_capacity) |>
        unique()
      
      # filter to each waterbody in turn
      metrics_observed_wb <- metrics_observed_sp |>
        filter(waterbody == waterbodies[j])
      metrics_counterfactual_wb <- metrics_counterfactual_sp |>
        filter(waterbody == waterbodies[j])
      
      # specify initial conditions
      initial <- specify_initial_conditions(
        species = species_list[i],
        waterbody = waterbodies[j],
        cpue = cpue,
        start = min(metrics_observed_wb$water_year),
        nsim = nsim,
        k = k
      )
      
      # grab stocking info if required (default to zero, otherwise)
      n_stocked <- rep(0, nrow(metrics_observed_wb))
      system_lu <- c(
        "broken_creek_r4" = "Broken Creek",
        "broken_river_r3" = "Broken River",
        "campaspe_river_r4" = "Campaspe River",
        "goulburn_river_r4" = "Goulburn River",
        "loddon_river_r4" = "Loddon River",
        "ovens_river_r5" = "Ovens River"
      )
      if (species_list[i] == "maccullochella_peelii") {
        stocking_rates <- stocking |>
          filter(
            Species == "Murray Cod",
            System == system_lu[waterbodies[j]]
          ) |>
          select(Year, Number) |>
          mutate(Year = Year + 1) |>
          rename(water_year = Year, number_stocked = Number)
        n_stocked <- metrics_observed_wb |>
          select(water_year) |>
          left_join(stocking_rates, by = c("water_year")) |>
          mutate(number_stocked = ifelse(is.na(number_stocked), 0, number_stocked)) |>
          pull(number_stocked)
      }
      
      # initialise population model for this species and waterbody
      pop <- specify_pop_model(
        species = species_list[i],
        waterbody = waterbodies[j],
        ntime = nrow(metrics_observed_wb), 
        nstocked = n_stocked,
        k = k
      )
      
      # and simulate from this model
      sims_observed <- simulate_scenario(
        species = species_list[i],
        x = pop, 
        nsim = nsim, 
        init = initial,
        metrics = metrics_observed_wb,
        coefs = get_coefs(species_list[i], waterbodies[j]),
        nburnin = nburnin
      )
      sims_counterfactual <- simulate_scenario(
        species = species_list[i],
        x = pop, 
        nsim = nsim, 
        init = initial,
        metrics = metrics_counterfactual_wb,
        coefs = get_coefs(species_list[i], waterbodies[j]),
        nburnin = nburnin
      )
      
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
      for (ff in seq_len(nrow(future_sub))) {
        
        # extract initial conditions from sims_observed
        initial_future <- sims_observed[, , dim(sims_observed)[3]]
        
        # pull out metrics for a given scenario
        metrics_future_sub <- metrics_future |>
          filter(
            species == species_list[i],
            waterbody == waterbodies[j],
            future == future_sub$future[ff],
            future_next == future_sub$future_next[ff],
            scenario == future_sub$scenario[ff],
            scenario_next == future_sub$scenario_next[ff]
          )
        
        # rename a few metrics for some species
        if (species_list[i] == "gadopsis_marmoratus") {
          metrics_future_sub <- metrics_future_sub |>
            rename(antecedent_flow = proportional_antecedent_flow)
        }
        if (species_list[i] == "maccullochella_peelii") {
          metrics_future_sub <- metrics_future_sub |>
            rename(blackwater_risk = hypoxia_risk)
        }
        if (species_list[i] == "melanotaenia_fluviatilis") {
          metrics_future_sub <- metrics_future_sub |>
            rename(redfin = presence_redfin, gambusia = presence_gambusia)
        }
        
        # simulate under the specific scenario
        sims_future <- simulate_scenario(
          species = species_list[i],
          x = pop, 
          nsim = nsim, 
          init = initial_future,
          metrics = metrics_future_sub,
          coefs = get_coefs(species_list[i], waterbodies[j]),
          nburnin = nburnin
        )
        
        # and save output
        future_name <- paste(
          future_sub$future[ff],
          future_sub$future_next[ff],
          future_sub$scenario[ff],
          future_sub$scenario_next[ff],
          sep = "_"
        )
        qsave(
          sims_future, 
          file = paste0("outputs/simulated/future-", species_list[i], "-", waterbodies[j], "-", future_name, ".qs")
        )
        
      }
      
    }
    
  }
  
}

# load all simulated models
mc_sim_obs <- load_simulated(type = "observed", species = "maccullochella")
mc_sim_cf <- load_simulated(type = "counterfactual", species = "maccullochella")
mc_sim_future <- load_simulated(type = "future", species = "maccullochella")
bf_sim_obs <- load_simulated(type = "observed", species = "gadopsis")
bf_sim_cf <- load_simulated(type = "counterfactual", species = "gadopsis")
bf_sim_future <- load_simulated(type = "future", species = "gadopsis")
rb_sim_obs <- load_simulated(type = "observed", species = "melanotaenia")
rb_sim_cf <- load_simulated(type = "counterfactual", species = "melanotaenia")
rb_sim_future <- load_simulated(type = "future", species = "melanotaenia")

# model CPUE using an AR1 model to estimate values (simple AR1 model with
#   random terms to soak up variation)
iter <- 2000
warmup <- 1000
chains <- 4
cores <- 4
use_cached <- FALSE
cpue_mc <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_bf <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Gadopsis marmoratus",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_rb <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Melanotaenia fluviatilis",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       

# posterior checks (saved to figures)
pp_mc <- pp_check(cpue_mc) + scale_x_log10()
pp_bf <- pp_check(cpue_bf) + scale_x_log10()
pp_rb <- pp_check(cpue_rb) + scale_x_log10()
ggsave(
  filename = "outputs/figures/pp-check-mc-cpue.png",
  plot = pp_mc,
  device = ragg::agg_png,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  filename = "outputs/figures/pp-check-bf-cpue.png",
  plot = pp_bf,
  device = ragg::agg_png,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  filename = "outputs/figures/pp-check-rb-cpue.png",
  plot = pp_rb,
  device = ragg::agg_png,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# model validation (using observed minus outputs)
mc_sim_metrics <- calculate_metrics(
  x = mc_sim_obs,
  cpue_mod = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
bf_sim_metrics <- calculate_metrics(
  x = bf_sim_obs,
  cpue = cpue_bf,
  subset = 1:11, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
rb_sim_metrics <- calculate_metrics(
  x = rb_sim_obs,
  cpue = cpue_rb,
  subset = 1:5, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
plot_metric(mc_sim_metrics)
plot_metric(bf_sim_metrics)
plot_metric(rb_sim_metrics)
plot_hindcasts(
  x = mc_sim_obs,
  cpue = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
plot_hindcasts(
  x = bf_sim_obs,
  cpue = cpue_bf,
  subset = 1:11,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
plot_hindcasts(
  x = rb_sim_obs,
  cpue = cpue_rb,
  subset = 1:5, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)

# comparison of flows with and without e-water (using observed- and counterfactual- outputs)
mc_adults <- plot_trajectories(
  x = mc_sim_obs, 
  y = mc_sim_cf,
  subset = 5:50,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9),
  scenarios = c("Observed", "Counterfactual")
) 
mc_recruits <- plot_trajectories(
  x = mc_sim_obs, 
  y = mc_sim_cf,
  subset = 1:2,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9),
  scenarios = c("Observed", "Counterfactual")
) 
bf_adults <- plot_trajectories(
  x = bf_sim_obs, 
  y = bf_sim_cf,
  subset = 1:11,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9),
  scenarios = c("Observed", "Counterfactual")
) 
rb_adults <- plot_trajectories(
  x = rb_sim_obs, 
  y = rb_sim_cf,
  subset = 1:5,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9),
  scenarios = c("Observed", "Counterfactual")
) 

# futures (using future- outputs)
# TODO: automate over all systems
mc_adult_futures_glb <- plot_forecasts(
  x = mc_sim_future, 
  subset = 5:10,
  probs = c(0.1, 0.9),
  system = "goulburn_river_r4"
) 
mc_adult_futures_cmp <- plot_forecasts(
  x = mc_sim_future, 
  subset = 5:10,
  probs = c(0.1, 0.9),
  system = "campaspe_river_r4"
) 

