# function to load flow data or download it from WMIS if required
fetch_flow <- function(start, end, recompile = FALSE) {
  
  # load flow data already there, download otherwise
  flow_exists <- grepl("flow-daily-compiled", dir("data"))
  if (any(flow_exists) & !recompile) {
    
    # list all and sort to newest
    flow_files <- dir("data")[flow_exists]
    flow_files <- sort(flow_files, decreasing = TRUE)[1]
    flow <- qread(paste0("data/", flow_files))
    
  } else {
    
    # check if the WMIS data have been downloaded already, load saved
    #   version if so
    flow_list_exists <- grepl("flow-list", dir("data"))
    if (any(flow_list_exists) & !recompile) {
      
      # list all and load newest      
      flow_list_file <- dir("data")[flow_list_exists]
      flow_list_file <- sort(flow_list_file, decreasing = TRUE)[1]
      flow <- qread(paste0("data/", flow_list_file))
      
    } else {
      
      # list target reaches and associated flow/temp gauges
      targets <- c(
        "broken_creek_r4" = 404210,
        "broken_river_r3" = 404224,
        "campaspe_river_r4" = 406202,
        "glenelg_river_r1" = 238210,
        "glenelg_river_r2" = 238211,
        "glenelg_river_r3" = 238206,
        "goulburn_river_r4" = 405200,
        "loddon_river_r2" = 407248,  # Tullaroop Creek, add McCallums Creek to this (following gauge)
        "loddon_river_r2b" = 407213,
        "loddon_river_r2_temp" = 407322, ## ONLY UNTIL 2015, full data at 407203
        "loddon_river_r4" = 407224,
        "loddon_river_r4_temp" = 407323,
        "loddon_river_r5" = 407202,
        "mackenzie_river_r1" = 415202,   ## NO TEMP DATA ANYWHERE NEARBY
        "macalister_river_r1" = 225247,  ## USE Licola (225209) if above Glenmaggie
        "macalister_river_r1_temp" = 225256,  ## ONLY UNTIL 2014
        "moorabool_river_r3" = 232204,
        "ovens_river_r5" = 403241,
        "thomson_river_r3" = 225231,
        "thomson_river_r3_temp" = 225212  ## WANDOCKA, DOWNSTREAM OF TARGET REACHES
      )
      
      # download data from WMIS
      flow <- vector("list", length = length(targets))
      names(flow) <- names(targets)
      for (i in seq_along(flow)) {
        
        # set variable codes by site
        varfrom <- c("100.00", "450.00")
        varto <- c("141.00", "450.00")
        
        # update for specific sites
        if (targets[i] == 404210)
          varfrom[1] <- "141.00"
        
        # grab data for each variable
        flow[[i]] <- fetch_hydro(
          sites = targets[i],
          start = dmy(paste0("01-01-", start)),
          end = dmy(paste0("31-12-", end)),
          options = list(varfrom = varfrom, varto = varto),
          include_missing = TRUE
        )
        
      }
      
      # widen to put each variable in its own column
      flow <- lapply(
        flow,
        function(x) x %>% pivot_wider(
          id_cols = c("date_formatted"),
          names_from = variable_name,
          values_from = c(value)
        )
      )
      
      # save to file
      qsave(flow, file = paste0("data/flow-list-", Sys.Date(), ".qs"))
      
    }
    
    # merge temperature data from backup gauges into main file
    flow$loddon_river_r2$stream_discharge_mld <-
      flow$loddon_river_r2$stream_discharge_mld +
      flow$loddon_river_r2b$stream_discharge_mld
    flow$loddon_river_r2$water_temperature_c <-
      flow$loddon_river_r2_temp$water_temperature_c
    flow$loddon_river_r2$`NA` <- NULL
    flow$macalister_river_r1$water_temperature_c <-
      flow$macalister_river_r1_temp$water_temperature_c
    flow$macalister_river_r1$`NA` <- NULL
    flow$thomson_river_r3$water_temperature_c <-
      flow$thomson_river_r3_temp$water_temperature_c
    flow$thomson_river_r3$`NA` <- NULL
    
    # add temperature information for Mackenzie River
    flow$mackenzie_river_r1$water_temperature_c <-
      flow$loddon_river_r2$water_temperature_c
    flow$mackenzie_river_r1$`NA` <- NULL
    
    # and rename mackenzie as R3 not R1
    flow$mackenzie_river_r3 <- flow$mackenzie_river_r1
    flow$mackenzie_river_r1 <- NULL
    
    # fill gaps in gauges for a few remaining gauges
    flow$loddon_river_r4$water_temperature_c <-
      flow$loddon_river_r5$water_temperature_c
    flow$loddon_river_r4$`NA` <- NULL
    flow$loddon_river_r5 <- NULL
    
    # and remove the _temp information that is now not needed
    flow <- flow[!grepl("_temp", names(flow))]
    flow$loddon_river_r2b <- NULL
    
    # remove any NA dates
    flow <- lapply(flow, \(x) x[!is.na(x$date_formatted), ])
    
    # and fill gaps and cut off extreme values
    for (i in seq_along(flow)) {
      flow[[i]]$stream_discharge_mld <- impute_rolling(flow[[i]]$stream_discharge_mld, recursive = TRUE)
      flow[[i]]$water_temperature_c <- ifelse(flow[[i]]$water_temperature_c == 0, NA, flow[[i]]$water_temperature_c)
      flow[[i]]$water_temperature_c <- impute_year(
        flow[[i]]$water_temperature_c,
        date = flow[[i]]$date_formatted,
        threshold = 100
      )
      flow[[i]]$water_temperature_c <- impute_rolling(flow[[i]]$water_temperature_c, recursive = TRUE)
      flow[[i]]$water_temperature_c <- rm_extremes(flow[[i]]$water_temperature_c, range = c(3, 32))
    }
    
    # save to file
    qsave(flow, file = paste0("data/flow-daily-compiled-", Sys.Date(), ".qs"))
    
  }
  
  # and return
  flow
  
}

# function to remove extreme values
rm_extremes <- function(x, range) {
  x <- ifelse(x < range[1], range[1], x)
  x <- ifelse(x > range[2], range[2], x)
  x
}

# function calculate flow metrics from daily discharge and temperature data
calculate_flow_metrics <- function(
    data, carrying_capacity, hypoxia_risk = NULL
) {
  
  # which year range do we want?  
  year <- min(year(data$date_formatted)):max(year(data$date_formatted))
  
  # and grab some reference values for rescaling
  reference_discharge <- median(data$stream_discharge_mld, na.rm = TRUE)
  reference_temperature <- median(data$water_temperature_c, na.rm = TRUE)
  reference_max <- max(data$stream_discharge_mld, na.rm = TRUE)
  
  # pull out relevant variables
  data <- data %>% 
    select(
      stream_discharge_mld,
      date_formatted,
      water_temperature_c
    ) %>%
    rename(
      flow = stream_discharge_mld,
      date = date_formatted,
      temp = water_temperature_c
    )
  
  # calculate metrics
  out <- data %>% reframe(
    water_year = calculate(
      flow,
      date,
      resolution = survey(season = 9:11, subset = year), 
      na.rm = TRUE
    )$date,
    proportional_spring_flow = calculate(
      flow,
      date,
      resolution = survey(season = 9:11, subset = year), 
      na.rm = TRUE
    )$metric,
    relative_max_antecedent = calculate(
      flow, 
      date,
      resolution = survey(season = 7:18, lag = 1, subset = year + 1),
      fun = max,
      na.rm = TRUE
    )$metric,
    proportional_antecedent_flow = calculate(
      flow,
      date,
      resolution = survey(season = c(7:18), subset = year, lag = 1),
      na.rm = TRUE
    )$metric,
    proportional_summer_flow = calculate(
      flow,
      date,
      na.rm = TRUE,
      resolution = survey(season = c(12:15), subset = year)
    )$metric,
    proportional_winter_flow = calculate(
      flow, 
      date, 
      na.rm = TRUE, 
      resolution = survey(season = 6:8, subset = year)
    )$metric,
    spawning_temperature = calculate(
      temp, date, resolution = survey(season = 10:12, subset = year)
    )$metric,
    spawning_flow_variability = calculate(
      flow,
      date,
      resolution = survey(season = 10:12, subset = year),
      fun = rolling_range,
      lag = 3,
      na.rm = TRUE
    )$metric,
    minimum_daily_flow = calculate(
      flow, 
      date, 
      na.rm = TRUE, 
      resolution = survey(season = 7:18, subset = year),
      fun = min
    )$metric,
    nday_gt16 = calculate(
      temp,
      date,
      resolution = survey(season = 9:12, subset = year),
      fun = days_above,
      threshold = 16
    )$metric,
    nday_lt5 = calculate(
      temp,
      date,
      resolution = survey(season = 7:18, subset = year),
      fun = days_below,
      threshold = 5
    )$metric,
    nday_gt20 = calculate(
      temp,
      date,
      resolution = survey(season = 10:15, subset = year),
      fun = days_above,
      threshold = 20
    )$metric,
    nday_lt10 = calculate(
      temp,
      date,
      resolution = survey(season = 7:18, subset = year),
      fun = days_below,
      threshold = 10
    )$metric
  )
  
  # rescale based on reference values
  idx <- grep("proportional", colnames(out))
  out[, idx] <- out[, idx] / reference_discharge
  out[, "spawning_temperature"] <- out[, "spawning_temperature"] / reference_temperature
  out[, "relative_max_antecedent"] <- out[, "relative_max_antecedent"] / reference_max
  
  # rename max antecedent
  colnames(out) <- gsub("relative", "proportional", colnames(out))
  
  # calculate dynamic carrying capacity
  annual_median_discharge <- calculate(
    data$flow, 
    data$date,
    standardise = by_median(subset = year, na.rm = TRUE),
    na.rm = TRUE,
    resolution = survey(season = 7:18, subset = year)
  )$metric
  lagged_median <- c(
    annual_median_discharge[1:3],
    rolling_median(annual_median_discharge, lag = 3, na.rm = TRUE)
  )
  lagged_median <- lagged_median[-c(1, (length(lagged_median) - 1L):length(lagged_median))]
  kdyn <- sapply(
    lagged_median,
    function(x) x * carrying_capacity
  )
  kdyn[kdyn > 3 * carrying_capacity] <- 3 * carrying_capacity
  out <- cbind(out, "kdyn" = kdyn)
  
  # and calculate blackwater risk if triggers known
  if (!is.null(hypoxia_risk)) {
    blackwater_risk_flow <- calculate(
      data$flow,
      data$date,
      resolution = survey(season = 12:15, subset = year), 
      fun = blackwater_fn,
      trigger = hypoxia_risk[1],
      na.rm = TRUE
    )$metric
    blackwater_risk_temp <- calculate(
      data$temp,
      data$date,
      resolution = survey(season = 12:15, subset = year), 
      fun = blackwater_fn,
      trigger = hypoxia_risk[2],
      type = days_above,
      na.rm = TRUE
    )$metric
    blackwater_risk <- ifelse(
      blackwater_risk_flow & blackwater_risk_temp,
      1,
      0
    )
    out <- cbind(out, "blackwater_risk" = blackwater_risk)
  }
  
  # return
  as.data.frame(out)
  
}

# function to calculate a rolling median over specified lag
rolling_median <- function(x, lag, ...) {
  n <- length(x)
  idx <- sapply(rev(seq_len(lag)) - 1, function(x) rep(x, n))
  idx <- sweep(idx, 1, seq_len(n), "+")
  idx <- ifelse(idx > n, NA, idx)
  df <- matrix(x[idx], nrow = n)
  apply(df, 1, median, ...)
}

# function to calculate blackwater risk from a set of flow and
#    temperature triggers
blackwater_fn <- function(x, trigger, type = days_below, ...) {
  out <- 0
  if (type(x, trigger, ...) > 0)
    out <- 1
  out
}

# helper functions to define a flow scenario based on freshes or baseflows
get_season <- function(type, season) {
  season <- season
  if (type %in% c("none", "baseflow"))
    season <- c(1:12)
  season
}
get_duration <- function(type, duration) {
  duration <- duration
  if (type %in% c("none", "baseflow"))
    duration <- rep(Inf, 12)
  duration
}
get_magnitude <- function(type, magnitude) {
  magnitude <- magnitude
  if (type %in% c("none", "baseflow") & length(magnitude) == 1)
    magnitude <- rep(magnitude, 12)
  magnitude
}

# function to add scenarios to the of an observed flow sequence based
#   on years to resample and events to include
add_scenario <- function(x, scn, template_years, cutoff = NULL) {
  
  # set a default cutoff
  if (is.null(cutoff))
    cutoff <- max(x$date_formatted[!is.na(x$stream_discharge_mld)])
  
  # drop all values after the cutoff if specified
  x <- x %>% filter(date_formatted <= cutoff)
  
  # apply each flow event for a given starting year
  apply_flow_event(
    x$stream_discharge_mld,
    x$date_formatted,
    scn,
    template_years,
    x$water_temperature_c,
    cutoff = cutoff
  )
  
}

# internal function to apply a series of flow events based on specified
#   season, duration, and magnitude
apply_flow_event <- function(x, y, scn, template_years, w, cutoff) {
  
  # initialise outer list (for different years)
  out <- vector("list", length = length(template_years))
  
  # add some identifiers to this list
  names(out) <- names(template_years)
  
  # work through each added year
  for (i in seq_along(template_years)) {
    
    # initialise inner list (for different e-flows)
    out[[i]] <- vector("list", length = length(scn))
    
    # and add some identifiers to this
    names(out[[i]]) <- names(scn)
    
    # pull out data for the year to be replicated
    cutoff_start <- cutoff + days(1)
    year(cutoff_start) <- template_years[i]
    cutoff_end <- cutoff_start + years(1) - days(1)
    idy <- y %within% interval(cutoff_start, cutoff_end)
    xbase <- x[idy]
    ybase <- y[idy]
    wbase <- w[idy]
    
    # work through each different scenario
    for (j in seq_along(scn)) {
      
      # set an initial data set for this scenario
      xsub <- xbase
      ysub <- ybase
      wsub <- wbase
      
      # reformat scenario definitions so they can be applied generically        
      season_set <- get_season(names(scn)[j], season = scn[[j]]$season)
      duration_set <- get_duration(names(scn)[j], duration = scn[[j]]$duration)
      magnitude_set <- get_magnitude(names(scn)[j], magnitude = scn[[j]]$magnitude)
      
      # work through each flow "event" and modify flows accordingly
      for (k in seq_along(season_set)) {
        
        # filter to target month
        idx <- month(ysub) == season_set[k]
        
        # and subsample the days if event has a duration shorter than the month
        if (duration_set[k] < sum(idx)) {
          start <- sample(which(idx), size = 1)
          if ((start + duration_set[k]) > max(which(idx)))
            start <- max(which(idx)) - duration_set[k] + 1L
          idx <- idx & (seq_along(idx) %in% start:(start + duration_set[k] - 1L))
        }
        
        xsub[idx] <- ifelse(
          xsub[idx] < magnitude_set[k],
          magnitude_set[k],
          xsub[idx]
        )
        
      }
      
      # combine xsub into a clean data.frame, adding 6 months of NAs at the end
      #   to keep things complete (sequence ends in June)
      out[[i]][[j]] <- data.frame(
        date_formatted = ysub,
        stream_discharge_mld = xsub,
        water_temperature_c = wsub
      )
      
      # update the date to run from the day after the cutoff
      #  to one year later
      cutoff_start <- cutoff + days(1)
      cutoff_end <- cutoff_start + years(1) - days(1)
      new_date <- seq(cutoff_start, cutoff_end, by = 1)
      
      # but drop Feb 29 from the resampled values or the imputed dates
      #  if the number of days doesn't match
      if (length(new_date) > nrow(out[[i]][[j]]))
        new_date <- new_date[!(month(new_date) == 2 & day(new_date) == 29)]
      if (length(new_date) < nrow(out[[i]][[j]])) {
        out[[i]][[j]] <- out[[i]][[j]] %>%
          filter(!(month(date_formatted) == 2 & day(date_formatted) == 29))
      }
      
      # and add this new date to the output
      out[[i]][[j]]$date_formatted <- new_date
      
      # pull out any extra days to the end of the calendar year if required
      #   and fill with NAs (to keep aae.hydro::calculate happy)
      if (!(month(cutoff) == 12 & day(cutoff) == 31)) {
        fill_start <- cutoff + days(1)
        year(fill_start) <- max(template_years[i])
        fill_end <- fill_start
        month(fill_end) <- 12
        day(fill_end) <- 31
        idz <- y %within% interval(fill_start, fill_end)
        yextra <- y[idz]
        year(yextra) <- max(year(out[[i]][[j]]$date_formatted))
        na_extra <- rep(NA, sum(idz))
        out[[i]][[j]] <- rbind(
          out[[i]][[j]],
          data.frame(
            date_formatted = yextra,
            stream_discharge_mld = na_extra,
            water_temperature_c = na_extra
          )
        )
      }
      
    }
    
  }
  
  # return
  out
  
}

# define flow futures (scenarios to append to the initial data set to
#   demonstrate possible outcomes)
specify_flow_future <- function(x) {
  
  # specify flow events for each reach  
  system_details <- list(
    "broken_creek_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(250, 5), rep(40, 3), rep(250, 5))),
      "fresh" = list(
        magnitude = c(300, 450, 400), 
        duration = c(10, 8, 14),
        season = c(7, 10, 11)
      )
    ),
    "broken_river_r3" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 30),
      "fresh" = list(
        magnitude = 500, 
        duration = 5,
        season = 3
      )
    ),
    "campaspe_river_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(50, 5), rep(200, 6), 50)),
      "fresh" = list(
        magnitude = c(200, 1800, 1800, 200), 
        duration = c(3, 5, 7, 3),
        season = c(1, 7, 10, 12)
      )
    ),
    "glenelg_river_r1" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(15, 5), rep(100, 6), 15)),
      "fresh" = list(
        magnitude = c(250, 250, 250), 
        duration = c(2, 5, 4),
        season = c(7, 9, 11)
      )
    ),
    "glenelg_river_r2" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(25, 5), rep(160, 6), 25)),
      "fresh" = list(
        magnitude = c(300, 300, 300), 
        duration = c(2, 5, 4),
        season = c(7, 9, 11)
      )
    ),
    "glenelg_river_r3" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(80, 5), rep(400, 6), 80)),
      "fresh" = list(
        magnitude = c(500, 500, 500), 
        duration = c(2, 5, 4),
        season = c(7, 9, 11)
      )
    ),
    "goulburn_river_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 800),
      "fresh" = list(
        magnitude = c(5700, 7300, 10500, 8000),
        duration = c(2, 2, 8, 3),
        season = c(5, 7, 10, 11)
      )
    ),
    "loddon_river_r2" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 0),
      "fresh" = list(magnitude = 0)
    ),
    "loddon_river_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(50, 5), rep(100, 6), 50)),
      "fresh" = list(
        magnitude = c(450),
        duration = c(9),
        season = c(10)
      )
    ),
    "macalister_river_r1" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(90, 6), rep(300, 4), rep(90, 2))),
      "fresh" = list(
        magnitude = c(750, 1350),
        duration = c(8, 5),
        season = c(12, 10)
      )
    ),
    "mackenzie_river_r3" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 1),
      "fresh" = list(
        magnitude = c(55, 55, 120),
        duration = c(4, 6, 3),
        season = c(7, 11, 9)
      )
    ),
    "moorabool_river_r3" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(10, 5), rep(60, 6), 10)),
      "fresh" = list(
        magnitude = c(80, 160, 120),
        duration = c(5, 5, 5),
        season = c(7, 9, 11)
      )
    ),
    "ovens_river_r5" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 0),
      "fresh" = list(
        magnitude = c(80, 430, 260), 
        duration = c(2, 3, 3),
        season = c(2, 3, 4)
      )
    ),
    "thomson_river_r3" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 125),
      "fresh" = list(
        magnitude = c(800, 900), 
        duration = c(5, 7),
        season = c(9, 11)
      )
    )  
  )
  
  # define a annual flow scenarios for each future type (can chop and
  #   change these to create multi-year sequences)
  template_years <- c("ave" = 2019, "dry" = 2009, "wet" = 2011)
  out <- vector("list", length = length(system_details))
  for (i in seq_along(system_details)) {
    
    # pull out flow for a specific system
    xsub <- x[[names(system_details)[i]]]
    
    # and append the relevant scenarios
    out[[i]] <- add_scenario(
      xsub, 
      system_details[[i]],
      template_years = template_years, 
      cutoff = dmy("30-06-2023")
    )
    
  }
  
  # add some identifiers to the output
  names(out) <- names(system_details)
  
  # return
  out
  
}
