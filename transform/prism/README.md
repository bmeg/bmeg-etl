# Notes on PRISM Data

## secondary-screen-dose-response-curve-parameters.csv

- Some AUCs are > 1. This is only the case when the ic50 value is NaN.
  - https://github.com/broadinstitute/repurposing/blob/master/secondary_processing_script.R#L32-L37

## primary-screen-replicate-collapsed-logfold-change.csv / secondary-screen-replicate-collapsed-logfold-change.csv

- fold change values are actually log2(Viability)
  - https://github.com/broadinstitute/repurposing/blob/master/secondary_processing_script.R#L238-L240
  
