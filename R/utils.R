# function to tidy up lists by adding names as a column and then binding all rows
add_col_bind <- function(x, name = "scenario") {
  id <- names(x)
  out <- mapply(\(.x, .y) .x |> mutate(category = .y), .x = x, .y = id, SIMPLIFY = FALSE)
  bind_rows(out) |> rename(!!name := category)
}
