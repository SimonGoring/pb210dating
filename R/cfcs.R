#' @title Empty dating function
#' @description Add stuff, but include that all uncertinty should be reported as one sigma, or instrument resolution.
#' @param constant A NULL value when measured supported values are used to calculate excess 210Pb.  If a constant background value is used, constant must be a two element numeric vector with the value and uncertainty.
#' @param inv_top In cases where a user may believe that some of the upper core is missing, approximate the missing proportion inventory.
#' @param inv_bot In cases where 210Pb does not meet the support values, approximate the missing proportion inventory to reach equilibrium.
#' @param runs The number of simulations to run for the MCMC.
#'
cfcs <- function(total_col,
                 support_col,
                 sample_date = NULL,
                 constants = NULL,
                 n_sect = length(s_mass),
                 constant = NULL,
                 runs = 100000) {

  # test, constant is either a null, or it is a two element vector, with a value and uncertainty.
  #sample_time needs to be a decimal year

  # All uncertainties should be one sigma, or instrument resolution.

  reduction_factor <- 100



}
