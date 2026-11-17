# Model hydrographs with empirical, flow-quantile-dependent error bands.
# One call: get_models(). When make == 1 it rebuilds the model parquet and the SVS-based errors;
# otherwise it reads both back from disk.

get_models <- function(make, xwalk, model_dir, out_file, error_file, svs_dir, SWORD_version,
                       probs, n_cores, bin_width = 0.05) {

  library(dplyr)
  library(tidyr)
  library(RNetCDF)
  library(arrow)
  library(parallel)


  # 1a. Download the newest SVS file (if not already present) and return its path -----------
  get_latest_svs <- function(dest_dir, doi = "doi:10.18419/DARUS-5843") {
    meta  <- jsonlite::fromJSON(paste0(
      "https://darus.uni-stuttgart.de/api/datasets/:persistentId/?persistentId=", doi))
    files <- meta$data$latestVersion$files$dataFile
    svs   <- files[grepl("^SVS_v[0-9_]+\\.nc$", files$filename), ]
    vers  <- numeric_version(gsub("_", ".", gsub("^SVS_v|\\.nc$", "", svs$filename)))
    newest <- svs[which(vers == max(vers))[1], ]

    dest <- file.path(dest_dir, newest$filename)
    if (!file.exists(dest)) {
      cat("Downloading", newest$filename, "\n")
      download.file(paste0("https://darus.uni-stuttgart.de/api/access/datafile/", newest$id),
                    dest, mode = "wb")
    } else cat("Already have", newest$filename, "\n")
    dest
  }


  # 1b. Read the SVS into one row per station x reach x day ----------------------------------
  read_svs <- function(svs_file, SWORD_version, start = as.Date("2023-01-01")) {
    nc <- open.nc(svs_file)
    on.exit(close.nc(nc))
    q   <- var.get.nc(nc, "Q")
    ymd <- var.get.nc(nc, "date_ymd")
    rid <- var.get.nc(nc, paste0("reach_id_v", SWORD_version))
    ids_raw <- var.get.nc(nc, "station_id", rawchar = TRUE)   # [station, strlen]
    ids <- trimws(apply(ids_raw, 1, function(r) rawToChar(r[r != as.raw(0)])))

    if (nrow(ymd) == 3) ymd <- t(ymd)
    dates <- as.Date(sprintf("%d-%02d-%02d", ymd[, 1], ymd[, 2], ymd[, 3]))
    if (nrow(q) != length(dates)) q <- t(q)          # q is [time, station]
    rid <- matrix(rid, ncol = ncol(q))                # one column per station

    reach_map <- data.frame(station_id = rep(ids, each = nrow(rid)),
                            reach_id   = as.vector(rid)) %>%
      filter(!is.na(reach_id))

    keep <- dates >= start
    q_long <- data.frame(station_id = rep(ids, each = sum(keep)),
                         date       = rep(dates[keep], ncol(q)),
                         gauge_q    = as.vector(q[keep, ])) %>%
      filter(!is.na(gauge_q))

    # station ids are not unique in the SVS: drop any station that duplicates a reach x day
    inner_join(reach_map, q_long, by = "station_id", relationship = "many-to-many") %>%
      group_by(station_id) %>%
      filter(!any(duplicated(paste(reach_id, date)))) %>%
      ungroup()
  }


  # Read the daily model ensembles ------------------------------------------------------------
  read_model_daily <- function(files,
                               cols = c(
                                 GradesHydroDL    = "ming_Q",
                                 dHBV2            = "chaopeng_Q",
                                 ML_ensemble      = "ML_ensemble",
                                 ERA5_rapid       = "ziyun_Q",
                                 CaMaFlood        = "cama_Q",
                                 GLOFAS_HRR       = "dongmei_Q",
                                 hymap            = "hymap_Q",
                                 physics_ensemble = "Trad_ensemble"
                               ),
                               n_cores) {

    ml_members   <- c("GradesHydroDL", "dHBV2")
    phys_members <- c("ERA5_rapid", "CaMaFlood", "GLOFAS_HRR", "hymap")
    members      <- c(ml_members, phys_members)
    value_cols   <- names(cols)

    first_non_na <- function(x) x[!is.na(x)][1]

    read_one <- function(f) {
      dt <- readRDS(f) |>
        dplyr::ungroup() |>
        dplyr::select(date, reach_id, dplyr::all_of(cols)) |>
        data.table::as.data.table()

      # two rows for one reach/date must never disagree on a model value
      conflicts <- dt[, lapply(.SD, data.table::uniqueN, na.rm = TRUE),
                      by = .(reach_id, date), .SDcols = members]
      if (any(conflicts[, ..members] > 1)) {
        stop("Conflicting values for the same reach/date in ", basename(f))
      }

      dt <- dt[, c(lapply(.SD, first_non_na), list(n_rows = .N)),
               by = .(reach_id, date), .SDcols = value_cols]

      # file ensembles on collapsed dates were built from partial rows: recompute
      dt[n_rows > 1, ML_ensemble      := rowMeans(.SD, na.rm = TRUE), .SDcols = ml_members]
      dt[n_rows > 1, physics_ensemble := rowMeans(.SD, na.rm = TRUE), .SDcols = phys_members]
      dt[is.nan(ML_ensemble),      ML_ensemble      := NA_real_]
      dt[is.nan(physics_ensemble), physics_ensemble := NA_real_]

      dt[, n_rows := NULL]
      dt
    }

    clust <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(clust))
    invisible(parallel::clusterEvalQ(clust, {
      library(dplyr)
      library(data.table)
      setDTthreads(1)
    }))

    out <- parallel::parLapply(clust, files, read_one)

    dplyr::as_tibble(data.table::rbindlist(out))
  }


  # 2-5. Pair models with gauges, index by model quantile, compute and save the errors -------
  build_quantile_errors <- function(svs_file, SWORD_version, ensemble_path, out_file, n_cores) {

    # gauges: one gauge per reach (reaches with more than one gauge are dropped)
    gauge_df <- read_svs(svs_file, SWORD_version) %>%
      transmute(reach_id = as.character(reach_id), station_id, date, gauge_Q = gauge_q) %>%
      group_by(reach_id) %>%
      filter(n_distinct(station_id) == 1) %>%
      ungroup() %>%
      select(-station_id)

    # models: only the reaches that have a gauge
    model_files <- list.files(ensemble_path, pattern = 'daily_ensemble\\.rds$', full.names = TRUE)
    model_ids   <- gsub('daily_ensemble\\.rds$', '', basename(model_files))
    model_files <- model_files[model_ids %in% unique(gauge_df$reach_id)]

    # quantile of each day on the model's own flow-duration curve (full record)
    models <- read_model_daily(model_files, n_cores = n_cores) %>%
      mutate(reach_id = as.character(reach_id)) %>%
      pivot_longer(-c(reach_id, date), names_to = 'model', values_to = 'model_Q') %>%
      filter(!is.na(model_Q)) %>%
      group_by(reach_id, model) %>%
      mutate(model_pct = ecdf(model_Q)(model_Q)) %>%
      ungroup()

    # same-day pairs in the science orbit, at least 20 per reach x model; directional error
    errors <- models %>%
      inner_join(gauge_df, by = c('reach_id', 'date'), relationship = 'many-to-one') %>%
      filter(date > as.Date('2023-06-21'), date < as.Date('2026-12-31'), gauge_Q > 0) %>%
      group_by(reach_id, model) %>%
      filter(n() >= 20) %>%
      ungroup() %>%
      mutate(rel_err = (model_Q - gauge_Q) / gauge_Q)

    # save every paired error so the summary can change without re-reading
    write_parquet(errors, out_file)
    errors
  }


  # 6. Pool errors across gauges by model quantile into a CI ---------------------------------
  summarise_quantile_errors <- function(errors, bin_width, probs) {
    errors %>%
      mutate(pct_bin = cut(model_pct, breaks = seq(0, 1, by = bin_width), include.lowest = TRUE)) %>%
      group_by(model, pct_bin) %>%
      summarise(n_pairs   = n(),
                n_reaches = n_distinct(reach_id),
                lower     = quantile(rel_err, probs[1]),
                median    = median(rel_err),
                upper     = quantile(rel_err, probs[2]),
                .groups   = 'drop')
  }


  # Apply the error band to every model day -------------------------------------------------
  add_model_CI <- function(model_hydrographs, error_ci, bin_width) {
    model_hydrographs %>%
      pivot_longer(-c(reach_id, date), names_to = "model", values_to = "model_Q") %>%
      filter(!is.na(model_Q)) %>%
      group_by(reach_id, model) %>%
      mutate(model_pct = ecdf(model_Q)(model_Q)) %>%
      ungroup() %>%
      mutate(pct_bin = cut(model_pct, breaks = seq(0, 1, by = bin_width), include.lowest = TRUE)) %>%
      left_join(select(error_ci, model, pct_bin, lower, upper),
                by = c("model", "pct_bin"), relationship = "many-to-one") %>%
      mutate(lower_CI = model_Q / (1 + upper),
             upper_CI = model_Q / (1 + lower)) %>%
      select(reach_id, date, model, model_Q, lower_CI, upper_CI)
  }


  # Run ---------------------------------------------------------------------------------------
  # make: 2 = rebuild model parquet + error parquet, 1 = rebuild error parquet only, 0 = read both
  if (make >= 1) {
    svs_file <- get_latest_svs(svs_dir)
    errors   <- build_quantile_errors(svs_file, SWORD_version, model_dir, error_file, n_cores)
  } else {
    errors <- read_parquet(error_file)
  }

  if (make == 2) {
    model_files     <- list.files(model_dir, pattern = 'daily_ensemble\\.rds$', full.names = TRUE)
    model_reach_ids <- gsub('daily_ensemble\\.rds$', '', basename(model_files))
    model_files     <- model_files[model_reach_ids %in% xwalk$reach_id]
    model_hydrographs <- read_model_daily(model_files, n_cores = n_cores)
    write_parquet(model_hydrographs, out_file)
  } else {
    model_hydrographs <- read_parquet(out_file)
  }

  error_ci <- summarise_quantile_errors(errors, bin_width, probs)

  add_model_CI(model_hydrographs, error_ci, bin_width)
}