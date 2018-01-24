# pb210dating

The `pb210dating` package is a package in development for CRAN and the R environment.  The package is designed to take data directly from instruments and calibrate and quantify time using both `cf` and `cfcs` methods of age calibration.

## Overview

The package uses several objects to manage the workflow from raw data to estimated ages.  These include a `constants` object.  This object may be read directly from file or may be added by a user.  These data include half-life constants, section depth, volume and mass and other variables that describe the samples and constants used to estimate age.

## Package Developers

  * Simon J Goring - University of Wisconsin - Madison
  * Joan Albert Sánchez Cabeza - UNAM

## Citations

Sanchez-Cabeza, J. A., & Ruiz-Fernández, A. C. (2012). ^210^Pb sediment radiochronology: an integrated formulation and classification of dating models. *Geochimica et Cosmochimica Acta*, **82**, 183-200. [DOI](http://dx.doi.org/10.1016/j.gca.2010.12.024)

Sanchez-Cabeza, J. A., & Ruiz-Fernández, A. C., Ontiveros-Cuadras, J. F., Bernal, L. H. P., Olid, C. (2014). Monte Carlo uncertainty calculation of 210 Pb chronologies and accumulation rates of # sediments and peat bogs. *Quaternary Geochronology*, **23**, 80-93.
