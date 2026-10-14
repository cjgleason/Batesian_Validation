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

