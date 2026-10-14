pull_site_data <- function(site_id, first_meas, last_meas,
                           out_dir     = "/nas/cee-water/cjgleason/colin/Batesian_Validation/stage_chunks",
                           param_stage = "00065",
                           param_q     = "00060",
                           pause       = 1) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(gsub("USGS-", "", site_id), ".parquet"))
  if (file.exists(out_file)) return(arrow::read_parquet(out_file))
  
  chunks     <- seq(as.Date(first_meas), as.Date(last_meas), by = "3 years")
  chunk_ends <- pmin(chunks + years(3) - days(1), as.Date(last_meas))
  
  # ── Continuous (15-min) ───────────────────────────────────────────────────
  continuous <- do.call(rbind, lapply(seq_along(chunks), function(j) {
    Sys.sleep(pause)
    tryCatch(
      read_waterdata_continuous(
        monitoring_location_id = site_id,
        parameter_code         = c(param_stage, param_q),
        time                   = c(format(chunks[j], "%Y-%m-%d"), format(chunk_ends[j], "%Y-%m-%d"))
      ),
      error = function(e) { message("FAILED continuous ", site_id, " chunk ", j, ": ", e$message); NULL }
    )
  }))
  
  # ── Daily ─────────────────────────────────────────────────────────────────
daily <- do.call(rbind, lapply(seq_along(chunks), function(j) {
    Sys.sleep(pause)
    tryCatch(
      read_waterdata_daily(
        monitoring_location_id = site_id,
        parameter_code         = param_q,
        statistic_id           = "00003",
        time                   = c(format(chunks[j], "%Y-%m-%d"), format(chunk_ends[j], "%Y-%m-%d")),
        skipGeometry           = TRUE
      ),
      error = function(e) { message("FAILED daily ", site_id, " chunk ", j, ": ", e$message); NULL }
    )
  }))
  # ── Pivot continuous to wide ──────────────────────────────────────────────
  continuous_wide <- continuous %>%
    mutate(
      date          = as.Date(time),
      unit_of_measure = as.character(unit_of_measure)  # strip units class
    ) %>%
    select(monitoring_location_id, time, date, parameter_code, value,
           unit_of_measure, approval_status, qualifier) %>%
    pivot_wider(names_from  = parameter_code,
                values_from = c(value, unit_of_measure),
                names_glue  = "{.value}_{parameter_code}") %>%
    rename(stage_ft = value_00065, q_cfs_15min = value_00060)
  
  # ── Join daily to continuous ──────────────────────────────────────────────
  daily_slim <- daily %>%
    mutate(unit_of_measure = as.character(unit_of_measure)) %>%
    select(monitoring_location_id, date = time, q_cfs_daily = value)
  
  result <- continuous_wide %>%
    left_join(daily_slim, by = c("monitoring_location_id", "date"))

    result <- result %>% mutate(across(where(is.numeric), as.numeric))
  
  arrow::write_parquet(result, out_file)
  result
}

param_stage  <- "00065"
param_q      <- "00060"
pause_sec    <- 0.5

site_ranges <- ADCP %>%
  group_by(site_id) %>%
  summarise(
    first_meas = min(as.Date(measurement_dateTime)),
    last_meas  = max(as.Date(measurement_dateTime)),
    .groups = "drop"
  )

#pull gauge Q and stage from gauges at sites with ADCP on same day as ADCP measurements. This limits API and data calls
for (i in seq_len(nrow(site_ranges))) {
  suppressMessages(pull_site_data(
    site_id    = site_ranges$site_id[i],
    first_meas = site_ranges$first_meas[i],
    last_meas  = site_ranges$last_meas[i]
  ))
}