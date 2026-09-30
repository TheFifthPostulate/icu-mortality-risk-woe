
convert_to_parquet <- function(file) {
    data <- arrow::read_csv_arrow(here::here(paste0(file, ".csv")))
    arrow::write_parquet(data, here::here(paste0(file, ".parquet")))
}
