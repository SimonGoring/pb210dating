#' @title Fill missing data with exponential fit.
#' @description When \eqn{^210}Pb or other isotopic data is missing at a depth, fill in values using an eponential fit.
#' @param x A \code{vector} or \code{data.frame} to interpolate.
#' @param depth The core depth.
#' @param s_mass The section mass for the sample, in order by depth from highest sample to deepest.
#' @import dplyr
#' @export

exp_interp <- function(x, depth, s_mass) {

  # Data checks
  ( c("Pb210", "Po210")  %in% colnames(x) &
      c("Pb214", "Ra226") %in% colnames(x))

  fixed_cols <- c("SampleCode","DepthLayer","DepthLayerUnc",
                  "MassSection","MassSectionUnc")

  assertthat::assert_that(!any(is.na(x$MassSection)),
                          msg = "Interpolation requires mass at all depths.")

  accum_d <- c(0, cumsum(x$MassSection))
  accum_mean <- (accum_d[1:(length(accum_d) - 1)] + accum_d[2:length(accum_d)]) / 2

  fill <- function(x){
    if(sum(!is.na(x)) == 1) {
      filled <- rep(x, length(x))
    } else {
      filled <- approx(x = accum_mean[!is.na(x)],
                       y = x[!is.na(x)] %>% log,
                       xout = accum_mean,
                       rule = 1)$y %>% exp
    }
    return(filled)
  }

  x[!colnames(x) %in% fixed_cols] <- apply( x[!colnames(x) %in% fixed_cols], 2, fill)

  # Needs to be modified to then return things with class "pbpo".

  return(x)
}
