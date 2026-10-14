base <- '/nas/cee-water/cjgleason/colin/Batesian_Validation/'

get_USGS <- function(make, usgs_dir, out_file) {
  if (make == 1) {
    read_USGS <- function(file) {
      read_parquet(file) %>%
        select(-parameter_code, -statistic_id, -last_modified, -qualifier, -geometry) %>%
        mutate(USGS_ID = as.character(monitoring_location_id)) %>%
        select(-monitoring_location_id)
    }
    USGS_df <- do.call(rbind, lapply(list.files(usgs_dir, full.names = TRUE), read_USGS))
    write_parquet(USGS_df, out_file)
  } else {
    USGS_df <- read_parquet(out_file)
  }
  USGS_df
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
    source(paste0(base, 'make_SOS_parquet.R'))
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

get_models <- function(make, xwalk, model_dir, out_file, n_cores) {
  if (make == 1) {
    source(paste0(base, 'make_ml_hydrographs_parquet.R'))
    model_files     <- list.files(model_dir, pattern = 'daily_ensemble\\.rds$', full.names = TRUE)
    model_reach_ids <- gsub('daily_ensemble\\.rds$', '', basename(model_files))
    model_files     <- model_files[model_reach_ids %in% xwalk$reach_id]
    model_hydrographs <- read_model_daily(model_files, n_cores = n_cores)
    write_parquet(model_hydrographs, out_file)
  } else {
    model_hydrographs <- read_parquet(out_file)
  }
  model_hydrographs
}

