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

