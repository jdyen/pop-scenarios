# function to tidy up lists by adding names as a column and then binding all rows
add_col_bind <- function(x, name = "scenario") {
  id <- names(x)
  out <- mapply(\(.x, .y) .x |> mutate(category = .y), .x = x, .y = id, SIMPLIFY = FALSE)
  bind_rows(out) |> rename(!!name := category)
}

# function to set initial conditions based on assumed constant proportional
#   survival over all age classes, but modified based on survey data
#   to account for shifts to older/younger age classes
set_initial <- function(
    species,
    cpue,
    cpue_max,
    effort_h,
    n = 10000,
    nsim = 1,
    rescale = NULL
) {
  
  # define a population matrix and work out steady state age frequency
  if (species == "gadopsis_marmoratus") {
    adults <- 2:11
    mat <- river_blackfish(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  if (species == "maccullochella_peelii") {
    adults <- 5:30
    mat <- murray_cod(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  if (species == "melanotaenia_fluviatilis") {
    adults <- 2:5
    mat <- murray_rainbowfish(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  
  # work out proportion of fish relative to max catch for a reach
  n_scaled <- 0
  if (cpue_max > 0)
    n_scaled <- n * (cpue / cpue_max)
  
  # rescale if required
  if (!is.null(rescale))
    n_scaled <- 0.6 * cpue * rescale  # current assumption is 60% of catch is adults
  
  # set a min value if cpue_adult = 0
  if (n_scaled == 0) {
    cpue_min <- min(c(cpue_max, 1)) / effort_h
    n_scaled <- n * (cpue_min / cpue_max)
  }
  
  # standardise this for adult abundances only
  initial_age_frequency <- initial_age_frequency / sum(initial_age_frequency[adults])
  
  # return
  matrix(
    rpois(
      nsim * length(initial_age_frequency),
      lambda = n_scaled * initial_age_frequency
    ),
    nrow = nsim,
    byrow = TRUE
  )
  
}
