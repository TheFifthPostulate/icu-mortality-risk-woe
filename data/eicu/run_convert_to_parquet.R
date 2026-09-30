library(arrow)
library(here)
source(here::here("convert_to_parquet.R"))

file_name <- "v2_features_signals_eicu"

convert_to_parquet(file_name)


