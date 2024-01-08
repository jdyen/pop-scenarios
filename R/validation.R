# helper functions for validation summaries
# function to estimate CPUE from catch data
estimate_cpue <- function(
    x,
    species, 
    iter = 2000,
    warmup = 1000, 
    chains = 2, 
    cores = 1, 
    use_cached = TRUE
) {
  
  # check species are in correct list
  if (!species %in% unique(x$scientific_name))
    stop("species must be represented in data set", call. = FALSE)
  
  # use saved version if available, otherwise refit the model
  species_clean <- gsub(" ", "_", tolower(species))
  if (use_cached & any(grepl(species_clean, dir("outputs/fitted/")))) {
    
    # if exists, load fitted version
    mod <- qread(paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
    
  } else {
    
    # subset data to target species and collate over all gear types and surveys
    #   within a year
    x <- x |>
      filter(scientific_name == species) |>
      group_by(id_site, waterbody, reach_no, survey_year) |>
      summarise(
        catch = sum(catch),
        effort_h = sum(effort_h)
      ) |>
      ungroup()
    
    # fit model
    mod <- stan_glmer(
      catch ~ log_cpue_ym1 +
        (1 | waterbody / reach_no) +
        (1 | id_site) +
        (1 | survey_year) +
        (1 | waterbody:survey_year) +
        offset(effort_h),
      family = poisson,
      data = x |> 
        left_join(
          x |>
            mutate(
              survey_year = survey_year + 1,
              log_cpue_ym1 = log(catch + 1) - log(effort_h)
            ) |>
            select(id_site, survey_year, log_cpue_ym1),
          by = c("id_site", "survey_year")
        ) |>
        filter(!is.na(log_cpue_ym1)),
      iter = iter,
      warmup = warmup,
      chains = chains,
      cores = cores
    )
    
    qsave(mod, file = paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
    
  }
  
  # return
  mod
  
}  

# function to calculate mid, lower, and upper bounds from simulated
#   trajectories
summarise_sim <- function(x, y, subset, probs, growth_rate = TRUE) {
  
  # pull out abundances for the subset of ages/stages
  abund <- apply(x[, subset, , drop = FALSE], c(1, 3), sum)
  
  # return growth rates if required (raw abundances returned otherwise)
  if (growth_rate) {
    
    # skip if all values are zero
    if (!all(abund == 0)) {
      
      # calculate the pop growth rate first, then summarise this if
      #   z-scaling
      # fill zeros with half of the min observed value
      abund[abund == 0] <- min(abund[abund > 0], na.rm = TRUE) / 2.0
      
      # calculate pop growth rate
      abund <- abund / abund[, c(1L, seq_len(ncol(abund) - 1L))]
      
      # calculate zscores
      abund <- abund - mean(abund)
      abund <- abund / sd(abund)
      
    }
    
    # remove first column (duplicated in previous for unknwn reason)
    abund <- abund[, -1]
    
  }
  
  # collate raw predicted values, dropping first column (which was included
  #   above for some unknown reason but is just year 1 divided by itself)
  out <- tibble(
    y,
    mid = apply(abund, 2, median),
    lower = apply(abund, 2, quantile, probs = probs[1]),
    upper = apply(abund, 2, quantile, probs = probs[2])
  )
  
  # return
  out
  
}

# function to calculate pop growth rates from obsered data and compare
#   to simulated values
add_cpue <- function(sim, cpue_mod, sim_years = 2010:2023, probs = c(0.1, 0.9)) {
  
  # generate new samples from the fitted posterior for each year/waterbody,
  #    setting previous cpue to 0 to estimate growth rate directly
  #    (no need to divide by catch_ym1)
  newdata <- cpue_mod$data |> 
    distinct(waterbody, reach_no, survey_year) |>
    mutate(
      log_cpue_ym1 = 0,
      effort_h = 1,
      id_site = "abc"
    )
  cpue_pred <- posterior_epred(
    cpue_mod, 
    newdata = newdata,
    re.form = ~ (1 | waterbody / reach_no) +
      (1 | survey_year)
  )
  cpue_ar1 <- tibble(
    newdata,
    cpue = apply(cpue_pred, 2, median),
    lower = apply(cpue_pred, 2, quantile, probs = probs[1]),
    upper = apply(cpue_pred, 2, quantile, probs = probs[2])
  )
  cpue_ar1 <- cpue_ar1 |>
    mutate(
      waterbody = paste0(
        tolower(gsub(" ", "_", waterbody)),
        "_r",
        reach_no
      )
    ) |>
    select(-reach_no) |>
    rename(growth_rate = cpue)
  
  # add in survey year info to simulated values
  nwaterbody <- sim |> pull(waterbody) |> unique() |> length()
  sim <- sim |> mutate(survey_year = rep(sim_years, nwaterbody))
  
  # z-scale it all
  cpue_std <- cpue_ar1 |>
    group_by(waterbody) |>
    summarise(
      center = mean(growth_rate, na.rm = TRUE),
      width = sd(growth_rate, na.rm = TRUE)
    )
  cpue_ar1 <- cpue_ar1 |>
    left_join(cpue_std, by = "waterbody") |>
    mutate(
      growth_rate_z = (growth_rate - center) / width,
      lower_z = (lower - center) / width,
      upper_z = (upper - center) / width
    ) |>
    select(waterbody, survey_year, growth_rate_z, lower_z, upper_z)
  
  # return this value joined to simulated pop growth rates
  sim |> 
    left_join(cpue_ar1, by = c("waterbody", "survey_year")) |>
    pivot_longer(
      cols = c(mid, growth_rate_z, lower, lower_z, upper, upper_z),
      values_to = "value",
      names_to = "type"
    ) |>
    mutate(
      category = ifelse(grepl("_z", type), "Observed", "Simulated"),
      type = gsub("_z", "", type),
      type = gsub("growth_rate", "mid", type)
    ) |>
    pivot_wider(
      id_cols = c(waterbody, survey_year, category),
      names_from = type,
      values_from = value
    )
  
}

# calculate summary metrics
calculate_val_metrics <- function(x, cpue_mod, subset, sim_years, probs = c(0.1, 0.9)) {
  
  # use functions above to summarise the simulated population trajectories
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nrow(x$scenario)),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add estimated CPUE
  x <- add_cpue(
    sim = x,
    cpue_mod = cpue_mod,
    sim_years = sim_years,
    probs = probs
  )
  
  # split out the modelled and observed values and calculate residuals
  x <- x |>
    select(waterbody, survey_year, category, mid) |>
    pivot_wider(
      id_cols = c(waterbody, survey_year),
      names_from = category,
      values_from = mid
    ) |>
    mutate(eps = Simulated - Observed)
  
  # calculate all the metrics
  x |>
    group_by(waterbody) |>
    summarise(
      r = ifelse(
        !all(is.na(Observed)), 
        cor(Simulated, Observed, use = "complete"),
        NA
      ),
      md = mean(eps, na.rm = TRUE),
      rmse = sqrt(mean(eps ^ 2, na.rm = TRUE)),
      sign = sum(sign(diff(Observed)) == sign(diff(Simulated)), na.rm = TRUE) / length(Observed)
    )
  
}

# helpers to tidy names in a plot
tidy_names <- function(x) {
  x <- gsub("_", " ", x)
  x <- strsplit(x, split = " ")
  y <- sapply(x, \(.x) gsub("r", "Reach ", .x[3]))
  x <- lapply(x, \(.x) paste0(toupper(substr(.x[1:2], 1, 1)), substr(.x[1:2], 2, nchar(.x))))
  paste(sapply(x, paste, collapse = " "), y, sep = ": ")
}
metric_names <- c(
  "r" = "r",
  "md" = "log(MD)",
  "rmse" = "log(RMSE)",
  "sign" = "Sign"
)

# function to plot validation metrics for one or more rivers at a time
plot_metric <- function(x) {
  
  # prepare data
  x <- x %>%
    pivot_longer(
      cols = c(r, md, rmse, sign),
      names_to = "name",
      values_to = "value"
    ) |>
    mutate(
      waterbody = tidy_names(waterbody),
      metric = metric_names[name],
      metric = factor(metric, levels = c("r", "log(MD)", "log(RMSE)", "Sign"))
    )
  
  # plot and return
  x %>%
    ggplot(aes(y = value, x = waterbody, fill = waterbody)) + 
    geom_bar(position = "dodge", stat = "identity") +
    ylab("Value") +
    facet_wrap( ~ metric, scales = "free") +
    scale_fill_brewer(palette = "Set2", name = "Waterbody") + 
    ggthemes::theme_hc() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 60, hjust = 1),
      strip.text = element_text(size = 8),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    )
  
}

# function to create abundance hindcast plots from simulated and observed data
plot_hindcasts <- function(x, cpue, subset, sim_years, probs = c(0.1, 0.9)) {
  
  # use functions above to summarise the simulated population trajectories
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nrow(x$scenario)),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add estimated CPUE
  x <- add_cpue(
    sim = x,
    cpue_mod = cpue,
    sim_years = sim_years,
    probs = probs
  )
  
  # set up base plot
  p <- x %>%
    ggplot(aes(x = survey_year, y = mid, col = category, group = category)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), 
      width = 0.2, 
      position = position_dodge(width = 0.2)
    ) +
    scale_color_brewer(
      name = "",
      palette = "Set2"
    ) +
    xlab("Water year") +
    ylab("Scaled population growth rate") +
    ggthemes::theme_hc() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    ) + 
    facet_wrap( ~ waterbody, scales = "free")
  
  # and return
  p
  
}

# function to plot abundance trajectories from two scenarios
plot_trajectories <- function(
    x, 
    y,
    subset,
    sim_years,
    probs = c(0.1, 0.9),
    scenarios = c("Observed", "Counterfactual")
) {
  
  # use functions above to summarise the simulated population trajectories
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nrow(x$scenario)),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = FALSE
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  y <- mapply(
    summarise_sim, 
    x = y$sims,
    y = lapply(
      seq_len(nrow(y$scenario)),
      \(i) y$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = FALSE
    ),
    SIMPLIFY = FALSE
  )
  y <- bind_rows(y)
  
  # add in survey year info to simulated values
  nwaterbody <- x |> pull(waterbody) |> unique() |> length()
  x <- x |> mutate(survey_year = rep(c(sim_years[1] - 1, sim_years), nwaterbody))
  nwaterbody <- y |> pull(waterbody) |> unique() |> length()
  y <- y |> mutate(survey_year = rep(c(sim_years[1] - 1, sim_years), nwaterbody))
  
  # combine both scenarios into a single object
  x <- bind_rows(
    x |> mutate(Scenario = scenarios[1]),
    y |> mutate(Scenario = scenarios[2])
  )
  
  # set up base plot
  p <- x %>%
    ggplot(aes(x = survey_year, y = mid, col = Scenario, group = Scenario)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), 
      width = 0.2, 
      position = position_dodge(width = 0.2)
    ) +
    scale_color_brewer(
      name = "",
      palette = "Set2"
    ) +
    xlab("Water year") +
    ylab("Abundance") +
    ggthemes::theme_hc() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    ) + 
    facet_wrap( ~ waterbody, scales = "free")
  
  # and return
  p
  
}

# function to plot one-step-ahead forecasts from start to final observed year
plot_forecasts <- function(x, subset, probs, system) {
  
  # use functions above to summarise the simulated population trajectories
  nscn <- nrow(x$scenario)
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nscn),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = FALSE
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add year information
  x <- x |>
    mutate(survey_year = rep(c(2023:2025), times = nscn))
  
  # clean up variable values
  x <- x |>
    mutate(
      future = factor(
        future,
        levels = c("dry", "ave", "wet"), 
        labels = c("Dry (2023/2024)", "Ave. (2023/2024)", "Wet (2023/2024)")
      ),
      future_next = factor(
        future_next,
        levels = c("dry", "ave", "wet"),
        labels = c("Dry (2024/2025)", "Ave. (2024/2025)", "Wet (2024/2025)")
      ),
      scenario = factor(
        scenario,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      ),
      scenario_next = factor(
        scenario_next,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      )
    )
  
  #  plot it
  p <- x |>
    filter(
      waterbody == system,
      survey_year == max(survey_year)
    ) |>
    mutate(survey_year = factor(survey_year)) |>
    ggplot(aes(y = mid, x = scenario, fill = scenario_next)) +
    geom_bar(position = position_dodge(0.9), stat = "identity") +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = position_dodge(0.9),
      col = "black",
      width = 0.2
    ) +
    xlab("Flow priority (2023/2024)") +
    ylab("Abundance") +
    scale_fill_brewer(name = "Flow priority (2024/2025)", palette = "Set2") +
    facet_grid(future_next ~ future, scales = "free") +
    ggthemes::theme_hc() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    )
  
  # and return
  p
  
}
