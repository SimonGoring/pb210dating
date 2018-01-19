#' @title Load constants for dating evaluation.
#' @description Load in datasets and associated data to help construct age models.
#' @param file

constants <- function(file, raw_data,
                      total_col, support_col,
                      NSectionDating = NULL) {

  if(is.null(NSectionDating)) NSectionDating <- nrow(raw_data)

  # Transform & fill missing values
  transformed_table <- exp_interp(raw_data)

  meta <- validate_textfile(file)

 constants = list(transformed_table,
                  meta,
                  NSectionDating,
               Lambda <- log(2) / use_lara()[1],
          DiameterUnc <- meta$DiameterUnc / 2 / sqrt(3),
              Surface <- pi * (meta$Diameter / 2)^2,
           DepthDelta <- diff(raw_data$DepthLayer),
         DepthSection <- rollmean(raw_data$DepthLayer, 2),
       MassSectionUnc <- raw_data$MassSectionUnc / 2 / sqrt(3),
              Density <- raw_data$MassSection / Surface / DepthDelta,
     MassDepthSection <- raw_data$MassSection / Surface,
    MassDepthLayerAcc <- c(0,cumsum(MassDepthSection)),
  MassDepthSectionAcc <- rollmean(MassDepthLayerAcc,2),
        Concentration <- raw_data[1 : NSectionsDating, total_col] - raw_data[1 : NSectionsDating, support_col]
  )

 class(constants) <- c('list', 'constant')

}
