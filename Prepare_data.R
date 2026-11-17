base <- '/nas/cee-water/cjgleason/colin/Batesian_Validation/'

get_USGS <- function(make, usgs_dir, out_file) {
  if (make == 1) {
    usgs <- open_dataset(list.files(usgs_dir, full.names = TRUE, pattern = "\\.parquet$")) %>%
      mutate(source_file = add_filename()) %>%
      select(-parameter_code, -statistic_id, -last_modified, -qualifier, -geometry) %>%
      rename(USGS_ID = monitoring_location_id) %>%
      collect() %>%
      mutate(file_ID = paste0("USGS-", sub("_discharge_uncertainty\\.parquet$", "", basename(source_file))))

    # each file is one gauge: any ID already present must match the file name
    stopifnot(all(is.na(usgs$USGS_ID) | usgs$USGS_ID == usgs$file_ID))

    usgs %>%
      mutate(USGS_ID = file_ID) %>%
      select(-file_ID, -source_file) %>%
      write_parquet(out_file)
  }
  read_parquet(out_file)
}

# map USGS IDs to SWORD v17 IDs; drop reaches with more than one gauge
make_xwalk <- function(adcp_file) {
  read.csv(adcp_file) %>%
    distinct(gage_id, v17_reach_id) %>%
    transmute(USGS_ID  = paste0("USGS-", sprintf("%08.0f", as.numeric(gage_id))),
              reach_id = sprintf("%.0f", v17_reach_id)) %>%
    group_by(reach_id) %>%
    filter(n() == 1) %>%
    ungroup()
}

get_SOS <- function(make, sos_dir, algos, out_file) {
  if (make == 1) {
      read_algo_q <- function(sos_file,
                        algos = c("hivdi", "metroman", "busboi", "sic4dvar", "momma", "consensus")) {

    library(RNetCDF)

  algo_getters <- list(
    hivdi     = function(s) s$hivdi$Q,
    sad       = function(s) s$sad$Qa,
    metroman  = function(s) s$metroman$allq,
    busboi    = function(s) s$busboi$q,
    sic4dvar  = function(s) s$sic4dvar$Q_da,
    momma     = function(s) s$momma$Q,
    consensus = function(s) s$consensus$consensus_q
  )


  cat("Reading:", basename(sos_file), "\n")
  nc  <- open.nc(sos_file)
  sos <- read.nc(nc, recursive = TRUE)
  close.nc(nc)

  epoch     <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
  reach_ids <- sos$reaches$reach_id
  t_lengths <- sapply(sos$reaches$time, length)

  pad_to <- function(q_list, target_lengths) {
    mapply(function(q, n) {
      if (length(q) == n) q else rep(NA_real_, n)
    }, q_list, target_lengths, SIMPLIFY = FALSE)
  }

swot_time <- epoch + unlist(sos$reaches$time)

  df <- data.frame(
    time     = swot_time,
    date     = as.Date(swot_time, tz = "UTC"),
        reach_id = rep(as.character(reach_ids), t_lengths)
  )

  for (a in algos) {
    df[[a]] <- unlist(pad_to(algo_getters[[a]](sos), t_lengths))
  }

  # fill values -> NA on the Q columns only
  df[algos] <- lapply(df[algos], function(x) replace(x, x < -1e11, NA))

  # row mask: keep a row if at least one requested algo has a value
  keep_rows <- rowSums(!is.na(df[algos])) > 0 


  df <- df[keep_rows, , drop = FALSE]


  rownames(df) <- NULL
  df
}



                      
    SoS_files <- list.files(sos_dir, full.names = TRUE, pattern = 'results')
    clust <- makeCluster(length(SoS_files))
    sos_hydrographs <- do.call(rbind, parLapply(clust, SoS_files, read_algo_q, algos = algos))
    stopCluster(clust)
    write_parquet(sos_hydrographs, out_file)
  } else {
    sos_hydrographs <- read_parquet(out_file)
  }
  sos_hydrographs
}



# keep only rows where every algorithm column in the table has a value (strict common dates)
keep_common_dates <- function(joined_df, algos) {
  algos <- intersect(algos, names(joined_df))
  joined_df %>%
    filter(if_all(all_of(algos), ~ !is.na(.x)))
}

                      
# SoS error bands from the V4 UAT validation: one nRMSE per algorithm (median across UAT reaches),
# applied as a constant absolute half-width per reach: nRMSE x mean(algorithm Q at that reach).
# Constant width means relatively wider bands at low flow. Algorithms with no UAT stat use metroman's.
add_sos_CI <- function(sos_hydrographs, uat_file, sos_algos) {
  library(dplyr)
  library(arrow)
 
  nrmse <- read_parquet(uat_file) %>%
    filter(stat == "nRMSE") %>%
    group_by(algo) %>%
    summarise(nRMSE = median(value), .groups = "drop")
 
  nrmse_algo <- setNames(nrmse$nRMSE[match(sos_algos, nrmse$algo)], sos_algos)
  nrmse_algo[is.na(nrmse_algo)] <- nrmse$nRMSE[nrmse$algo == "metroman"]
 
  for (a in sos_algos) {
    reach_mean <- ave(sos_hydrographs[[a]], sos_hydrographs$reach_id,
                      FUN = function(x) mean(x, na.rm = TRUE))
    half <- nrmse_algo[[a]] * reach_mean
    sos_hydrographs[[paste0(a, "_lower_CI")]] <- pmax(0, sos_hydrographs[[a]] - half)   # flow cannot be negative
    sos_hydrographs[[paste0(a, "_upper_CI")]] <- sos_hydrographs[[a]] + half
  }
  sos_hydrographs
}
 
