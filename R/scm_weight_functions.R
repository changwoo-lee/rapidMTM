### Weight functions for spatial clustering model ###

# sqrt weight function
lweight_fcn_sqrt <- function(log_post_new, log_post_curr) {
  return(0.5 * (log_post_new - log_post_curr))
}

# ordinal weight function
lweight_fcn_ord <- function(log_post_new, log_post_curr) { return(log_post_new) }

# min weight function
lweight_fcn_min <- function(log_post_new, log_post_curr) {
  return(min(0, log_post_new - log_post_curr))
}

# max weight function
lweight_fcn_max <- function(log_post_new, log_post_curr) {
  return(max(0, log_post_new - log_post_curr))
}
