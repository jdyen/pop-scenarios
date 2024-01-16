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
# Date modified: 16 January 2024

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
library(patchwork)

# and some load helpers
source("R/utils.R")
source("R/fish.R")
source("R/flow.R")
source("R/validation.R")

# settings
set.seed(2023-10-17)
nsim <- 200
nburnin <- 10
simulate_again <- TRUE

# work out reach lengths
vewh_reaches <- fetch_table("eflow_reaches_20171214", "projects") |>
  collect()
st_geometry(vewh_reaches) <- st_as_sfc(vewh_reaches$geom, crs = 4283)
vewh_reach_length <- vewh_reaches |>
  select(eflowriver, shape_leng, vewh_reach) |>
  mutate(segment_length = st_length(vewh_reaches)) |>
  filter(
    grepl(
      "broken|glenelg|goulburn|mackenzie|campaspe|loddon|ovens|macalister|moorabool|thomson",
      eflowriver,
      ignore.case = TRUE
    )
  ) |>
  group_by(eflowriver, vewh_reach) |>
  summarise(surveyed_length_m = sum(segment_length)) |>
  ungroup() |>
  mutate(surveyed_length_km = surveyed_length_m / 1000)

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

# download CPUE of recruits
cpue_recruits <- fetch_fish(recruit = TRUE, recompile = FALSE)

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
  reference = flow_futures,
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
    ),
    hypoxia_risk_temp_counterfactual = hypoxia_risk_temp,
    hypoxia_risk_counterfactual = ifelse(
      hypoxia_risk_discharge_counterfactual & hypoxia_risk_temp_counterfactual,
      1,
      0
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
    iwh = rescale_fn(iwh, base = 0.75, na.rm = TRUE),
    proportion_overhang_hab_metrics = rescale_fn(
      proportion_overhang_hab_metrics,
      base = 0.85,
      na.rm = TRUE
    )
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
      initial <- simulate_scenario(
        species = species_list[i],
        x = pop, 
        nsim = nsim, 
        init = initial,
        metrics = metrics_observed_wb[1, ],
        coefs = get_coefs(species_list[i], waterbodies[j]),
        nburnin = nburnin - 1
      )
      sims_observed <- simulate_scenario(
        species = species_list[i],
        x = pop, 
        nsim = nsim, 
        init = initial[, , dim(initial)[3]],
        metrics = metrics_observed_wb,
        coefs = get_coefs(species_list[i], waterbodies[j]),
        nburnin = 0
      )
      sims_counterfactual <- simulate_scenario(
        species = species_list[i],
        x = pop, 
        nsim = nsim, 
        init = initial[, , dim(initial)[3]],
        metrics = metrics_counterfactual_wb,
        coefs = get_coefs(species_list[i], waterbodies[j]),
        nburnin = 0
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
      
      # extract initial conditions for forecasts from sims_observed
      initial_future <- sims_observed[, , dim(sims_observed)[3]]
      
      # simulate futures
      future_sub <- metrics_future |>
        filter(
          species == species_list[i],
          waterbody == waterbodies[j]
        ) |>
        distinct(future, future_next, scenario, scenario_next)
      for (ff in seq_len(nrow(future_sub))) {
        
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
          nburnin = 0
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
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
use_cached <- TRUE
cpue_mc <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_recruit_mc <- estimate_cpue(
  x = cpue_recruits, 
  recruit = TRUE,
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
cpue_recruit_bf <- estimate_cpue(
  x = cpue_recruits, 
  recruit = TRUE,
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
pp_mc <- pp_check(cpue_mc) + scale_x_log10() + xlab("Catch (total)") + ylab("Density") + theme(legend.position = "none")
pp_bf <- pp_check(cpue_bf) + scale_x_log10() + xlab("Catch (total)") + ylab("Density") + theme(legend.position = "none")
pp_mc_recruit <- pp_check(cpue_recruit_mc) + scale_x_log10() + xlab("Catch (young of year)") + ylab("Density") + theme(legend.position = "none")
pp_bf_recruit <- pp_check(cpue_recruit_bf) + scale_x_log10() + xlab("Catch (young of year)") + ylab("Density") + theme(legend.position = "none")
pp_rb <- pp_check(cpue_rb) + scale_x_log10() + xlab("Catch (total)") + ylab("Density") + theme(legend.position = "none")
pp_all <- (pp_mc | pp_mc_recruit) /
  (pp_bf | pp_bf_recruit) /
  (pp_rb | plot_spacer()) +
  plot_annotation(tag_levels = "a")
ggsave(
  filename = "outputs/figures/pp-checks.png",
  plot = pp_all,
  device = ragg::agg_png,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# model validation (using observed minus outputs)
mc_sim_metrics <- calculate_val_metrics(
  x = mc_sim_obs,
  cpue_mod = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
mc_sim_metrics_adult <- calculate_val_metrics(
  x = mc_sim_obs,
  cpue_mod = cpue_mc,
  subset = 5:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
mc_sim_metrics_recruit <- calculate_val_metrics(
  x = mc_sim_obs,
  cpue_mod = cpue_recruit_mc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year)
)
bf_sim_metrics <- calculate_val_metrics(
  x = bf_sim_obs,
  cpue = cpue_bf,
  subset = 1:11, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
bf_sim_metrics_recruit <- calculate_val_metrics(
  x = bf_sim_obs,
  cpue_mod = cpue_recruit_bf,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year)
)
rb_sim_metrics <- calculate_val_metrics(
  x = rb_sim_obs,
  cpue = cpue_rb,
  subset = 1:7, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
sim_metrics <- bind_rows(
  mc_sim_metrics |> mutate(species = "Murray Cod"),
  mc_sim_metrics_recruit |> mutate(species = "Murray Cod (young of year)"),
  mc_sim_metrics_adult |> mutate(species = "Murray Cod (adults)"),
  rb_sim_metrics |> mutate(species = "Murray-Darling Rainbowfish")
)
metrics_plot_mdb <- plot_metric(sim_metrics)
metrics_plot_bf <- plot_metric(
  bind_rows(
    bf_sim_metrics |> mutate(species = "River Blackfish"),
    bf_sim_metrics_recruit |> mutate(species = "River Blackfish (young of year)")
  )
)
ggsave(
  filename = "outputs/figures/metrics-mdb.png",
  plot = metrics_plot_mdb,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/metrics-bf.png",
  plot = metrics_plot_bf,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# plot all hindcast combinations
mc_hindcast <- plot_hindcasts(
  x = mc_sim_obs,
  cpue = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
mc_hindcast_adult <- plot_hindcasts(
  x = mc_sim_obs,
  cpue = cpue_mc,
  subset = 5:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
mc_hindcast_recruit <- plot_hindcasts(
  x = mc_sim_obs,
  cpue = cpue_recruit_mc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year)
)
bf_hindcast <- plot_hindcasts(
  x = bf_sim_obs,
  cpue = cpue_bf,
  subset = 1:11,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
bf_hindcast_recruit <- plot_hindcasts(
  x = bf_sim_obs,
  cpue = cpue_recruit_bf,
  recruit = TRUE,
  subset = 1,
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year)
)
rb_hindcast <- plot_hindcasts(
  x = rb_sim_obs,
  cpue = cpue_rb,
  subset = 1:7, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
ggsave(
  filename = "outputs/figures/hindcast-mc.png",
  plot = mc_hindcast,
  device = ragg::agg_png,
  width = 6.5,
  height = 5,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcast-mc-recruits.png",
  plot = mc_hindcast_recruit,
  device = ragg::agg_png,
  width = 6.5,
  height = 5,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcast-mc-adults.png",
  plot = mc_hindcast_adult,
  device = ragg::agg_png,
  width = 6.5,
  height = 5,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcast-bf.png",
  plot = bf_hindcast,
  device = ragg::agg_png,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcast-bf-recruit.png",
  plot = bf_hindcast_recruit,
  device = ragg::agg_png,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcast-rb.png",
  plot = rb_hindcast,
  device = ragg::agg_png,
  width = 6.5,
  height = 5,
  units = "in",
  dpi = 600
)


# comparison of flows with and without e-water (using observed- and counterfactual- outputs)
mc_all <- plot_trajectories(
  x = mc_sim_obs, 
  y = mc_sim_cf,
  subset = 1:50,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9)
) 
mc_adults <- plot_trajectories(
  x = mc_sim_obs, 
  y = mc_sim_cf,
  subset = 5:50,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9)
) 
mc_recruits <- plot_trajectories(
  x = mc_sim_obs, 
  y = mc_sim_cf,
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9)
) 
bf_adults <- plot_trajectories(
  x = bf_sim_obs, 
  y = bf_sim_cf,
  subset = 2:11,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9)
) 
bf_all <- plot_trajectories(
  x = bf_sim_obs, 
  y = bf_sim_cf,
  subset = 1:11,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9)
) 
bf_recruits <- plot_trajectories(
  x = bf_sim_obs, 
  y = bf_sim_cf,
  subset = 1,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9)
) 
rb_all <- plot_trajectories(
  x = rb_sim_obs, 
  y = rb_sim_cf,
  subset = 1:7,
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year),
  probs = c(0.1, 0.9)
) 

metrics_counterfactual |>
  filter(
    waterbody == "broken_creek_r4",
    species == "maccullochella_peelii"
  ) |>
  mutate(type = "counterfactual") |>
  bind_rows(
    metrics_observed |>
      filter(
        waterbody == "broken_creek_r4",
        species == "maccullochella_peelii"
      ) |>
      mutate(type = "actual")
  ) |>
  select(
    water_year,
    type,
    proportional_spring_flow,
    proportional_summer_flow,
    proportional_winter_flow,
    proportional_max_antecedent,
    spawning_flow_variability,
    spawning_temperature,
    minimum_daily_flow,
    lagged_median_flow,
    hypoxia_risk
  ) |>
  pivot_longer(
    cols = c(contains("flow"), contains("max"), contains("hyp"), contains("temp"))
  ) |>
  ggplot(aes(x = water_year, y = value, col = type)) +
  geom_point() +
  facet_wrap( ~ name, scales = "free_y")

# futures (using future- outputs)
mc_adult_futures_glb <- plot_forecasts(
  x = mc_sim_future, 
  subset = 5:50,
  probs = c(0.1, 0.9),
  system = "goulburn_river_r4"
) 

# Options: 2024 forecast only
mc_recruit_futures_glb <- plot_forecasts(
  x = mc_sim_future, 
  subset = 1,
  probs = c(0.1, 0.9),
  system = "goulburn_river_r4",
  target = 2024
) 

# 2025 forecast but for a single 2024 forecast
mc_recruit_futures_glb <- plot_forecasts(
  x = mc_sim_future, 
  subset = 1:3,
  probs = c(0.1, 0.9),
  system = "goulburn_river_r4",
  climate = "Wet (2023/2024)"
) 

mc_adult_futures_cmp <- plot_forecasts(
  x = mc_sim_future, 
  subset = 5:50,
  probs = c(0.1, 0.9),
  system = "campaspe_river_r4"
) 
mc_recruit_futures_cmp <- plot_forecasts(
  x = mc_sim_future, 
  subset = 1:3,
  probs = c(0.1, 0.9),
  system = "campaspe_river_r4"
) 

plot_forecasts(
  x = mc_sim_future, 
  subset = 1,
  probs = c(0.1, 0.9),
  system = "campaspe_river_r4",
  climate = "Wet (2023/2024)"
) 
plot_forecasts(
  x = mc_sim_future, 
  subset = 5:50,
  probs = c(0.1, 0.9),
  system = "campaspe_river_r4",
  target = 2024
) 

