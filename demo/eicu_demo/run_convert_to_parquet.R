library(arrow)
library(here)
source(here::here("convert_to_parquet.R"))

file_name <- "v2_severity_eicu_demo"

convert_to_parquet(file_name)


