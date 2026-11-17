validate_algo <- function(algo_name, joined_df, truth_col, lower_col, upper_col, xwalk, val_type, truth_type) {
  val_type <- match.arg(val_type, c("trad", "Batesian"))

  # an algorithm with no data in this table is an all-NA column: it fails the pair minimum below
  if (!algo_name %in% names(joined_df)) joined_df[[algo_name]] <- NA_real_

  # the algorithm's own interval, if it has one (models: <algo>_lower_CI / <algo>_upper_CI);
  # algorithms without one (SoS) get NA, so their interval metrics come out NA
  algo_lower <- paste0(algo_name, "_lower_CI")
  algo_upper <- paste0(algo_name, "_upper_CI")
  if (!algo_lower %in% names(joined_df)) joined_df[[algo_lower]] <- NA_real_
  if (!algo_upper %in% names(joined_df)) joined_df[[algo_upper]] <- NA_real_

  algo_df <- joined_df %>%
    select(reach_id,
           time,
           Q_truth       = all_of(truth_col),
           Q_lower       = all_of(lower_col),
           Q_upper       = all_of(upper_col),
           Q_algo        = all_of(algo_name),
           Q_algo_lower  = all_of(algo_lower),
           Q_algo_upper  = all_of(algo_upper))

  per_reach <- lapply(split(algo_df, algo_df$reach_id), function(reach_df) {
    switch(val_type,
      trad     = trad_validate(truth = reach_df$Q_truth, prediction = reach_df$Q_algo),
      Batesian = batesian_validate(time       = reach_df$time,
                                   truth      = reach_df$Q_truth,
                                   prediction = reach_df$Q_algo,
                                   lower      = reach_df$Q_lower,
                                   upper      = reach_df$Q_upper,
                                   pred_lower = reach_df$Q_algo_lower,
                                   pred_upper = reach_df$Q_algo_upper))
  })

  # no reach had enough pairs: nothing to report for this algorithm (do.call(rbind) skips NULL)
  if (all(vapply(per_reach, is.null, logical(1)))) return(NULL)

  bind_rows(per_reach, .id = "reach_id") %>%
    left_join(xwalk, by = "reach_id", relationship = "many-to-one") %>%
    mutate(algo = algo_name, type = truth_type, .after = reach_id) %>%
    relocate(USGS_ID, .after = reach_id)
}



trad_validate <- function(truth, prediction) {

  ok         <- !is.na(truth) & !is.na(prediction)
  truth      <- truth[ok]
  prediction <- prediction[ok]
  if (length(truth) < 20) return(NULL)

  sd_pop <- function(x) sqrt(mean((x - mean(x))^2))   # np.std (divide by n)

  r     <- cor(prediction, truth)
  RMSE  <- sqrt(mean((prediction - truth)^2))
  nBIAS <- (sum(prediction - truth) / length(truth)) / mean(truth)

  tibble(
    r         = r,
    Rsq       = r^2,
    NSE       = if (sum((truth - mean(truth))^2) == 0) NA_real_ else
                  1 - sum((truth - prediction)^2) / sum((truth - mean(truth))^2),
    KGE       = if (is.na(r) || sd_pop(truth) == 0) NA_real_ else
                  1 - sqrt((r - 1)^2 + (sd_pop(prediction) / sd_pop(truth) - 1)^2 +
                           (mean(prediction) / mean(truth) - 1)^2),
    RMSE      = RMSE,
    nRMSE     = RMSE / mean(truth),
    nBIAS     = nBIAS,
    abs_nBIAS = abs(nBIAS)
  )
}



batesian_validate <- function(time, truth, prediction, lower, upper, pred_lower, pred_upper) {
  ok         <- !is.na(truth) & !is.na(prediction) & !is.na(lower) & !is.na(upper)
  time       <- time[ok]
  truth      <- truth[ok]
  prediction <- prediction[ok]
  lower      <- lower[ok]
  upper      <- upper[ok]
  pred_lower <- pred_lower[ok]
  pred_upper <- pred_upper[ok]
  if (length(prediction) < 20) return(NULL)

  is_OP <- prediction > upper
  is_UP <- prediction < lower
  is_IP <- prediction >= lower & prediction <= upper
  miss  <- is_OP | is_UP

  # Accuracy  = IP / (IP + OP + UP): hit rate of the point prediction against the gauge band.
  # Miss_bias = (OP - UP) / N: 0 = balanced or perfect, +1 = always over, -1 = always under
  # Miss_magnitude = mean log distance outside the band edge crossed, over misses only:
  #   ln(model / upper) for OP, ln(lower / model) for UP. Symmetric: a factor-of-k miss scores
  #   ln k in either direction. 0 if there are no misses.
  tibble(
    Accuracy       = sum(is_IP) / (sum(is_OP) + sum(is_UP) + sum(is_IP)),
    Miss_bias      = (sum(is_OP) - sum(is_UP)) / length(prediction),
    Miss_magnitude = if (!any(miss)) 0 else
                       mean(c(log(prediction[is_OP] / upper[is_OP]), log(lower[is_UP] / prediction[is_UP]))),
    Bcorr_rho0     = batesian_corr(time, truth, prediction, lower, upper, rho = 0),
    Bcorr_rho08    = batesian_corr(time, truth, prediction, lower, upper, rho = 0.8)
  ) %>%
    bind_cols(interval_skill(lower, upper, pred_lower, pred_upper))
}



# Interval skill: the model interval vs the gauge interval, treated as a classification of the
# whole record (like summing pixels over an image), in linear discharge units, pooled per reach:
#   Precision = sum(overlap) / sum(model width)    how much of what the model claims is supported
#   Recall    = sum(overlap) / sum(gauge width)    how much of the plausible truth the model covers
#   F1        = 2 sum(overlap) / (sum(model width) + sum(gauge width))
#   IoU       = sum(overlap) / sum(union)
# Uses only days where the model has an interval; NA when it has none (e.g. SoS algorithms).
interval_skill <- function(lower, upper, pred_lower, pred_upper) {
  has <- !is.na(pred_lower) & !is.na(pred_upper)
  if (!any(has)) {
    return(tibble(Precision = NA_real_, Recall = NA_real_, F1 = NA_real_, IoU = NA_real_))
  }

  gauge_w <- upper[has] - lower[has]
  model_w <- pred_upper[has] - pred_lower[has]
  overlap <- pmax(0, pmin(upper[has], pred_upper[has]) - pmax(lower[has], pred_lower[has]))
  union   <- gauge_w + model_w - overlap

  tibble(Precision = sum(overlap) / sum(model_w),
         Recall    = sum(overlap) / sum(gauge_w),
         F1        = 2 * sum(overlap) / (sum(model_w) + sum(gauge_w)),
         IoU       = sum(overlap) / sum(union))
}



# Batesian correlation: fraction of overpass-to-overpass log changes where the model's
# change falls inside the 95% CI of the gauge's change. The gauge band is symmetric in
# log space, so sigma = (log(upper) - log(lower)) / (2 * 1.96) at each time, and
# SE(change) = sqrt(s1^2 + s2^2 - 2 * rho * s1 * s2). Bias cancels in the differencing.
batesian_corr <- function(time, truth, prediction, lower, upper, rho) {
  ord <- order(time)

  sig     <- (log(upper[ord]) - log(lower[ord])) / (2 * 1.96)
  d_gauge <- diff(log(truth[ord]))
  d_model <- diff(log(prediction[ord]))
  s1      <- head(sig, -1)
  s2      <- tail(sig, -1)

  se <- sqrt(s1^2 + s2^2 - 2 * rho * s1 * s2)
  mean(abs(d_model - d_gauge) <= 1.96 * se)
}



rank_algos <- function(val_df, tie_tol = 0.025) {
  higher_better <- c("r", "Rsq", "NSE", "KGE", "Accuracy", "Bcorr_rho0", "Bcorr_rho08",
                     "Precision", "Recall", "F1", "IoU")
  lower_better  <- c("RMSE", "nRMSE", "abs_nBIAS", "Miss_magnitude")
  target_better <- c(nBIAS = 0, Miss_bias = 0)
  # tie_tol: these tie on an absolute difference of tie_tol; all others on a relative difference of tie_tol
  absolute_tie  <- c("nBIAS", "abs_nBIAS", "Miss_bias")

  ranked <- val_df %>%
    pivot_longer(-c(reach_id, USGS_ID, algo, type), names_to = "metric", values_to = "value")

  stopifnot(all(ranked$metric %in% c(higher_better, lower_better, names(target_better))))

  # score: lower is always better; rank 1 = best algo for that reach/type/metric
  ranked %>%
    mutate(score = case_when(metric %in% higher_better        ~ -value,
                             metric %in% lower_better         ~  value,
                             metric %in% names(target_better) ~ abs(value - target_better[metric]))) %>%
    group_by(reach_id, type, metric) %>%
    mutate(rank = tie_rank(score, absolute = first(metric) %in% absolute_tie, tol = tie_tol)) %>%
    ungroup()
}



# Competition ranking with ties: sort scores (lower = better); each tie group is led by
# its best member, and anything within tol of that leader shares the leader's rank.
# The next group's rank is its position, so ranks 1,1 are followed by 3.
tie_rank <- function(score, absolute, tol = 0.025) {
  ord  <- order(score, na.last = NA)
  rank <- rep(NA_integer_, length(score))

  for (i in seq_along(ord)) {
    s <- score[ord[i]]
    if (i == 1 || s - leader > (if (absolute) tol else tol * abs(leader))) {
      leader  <- s
      current <- i
    }
    rank[ord[i]] <- current
  }
  rank
}



# ---- Monte Carlo validation ------------------------------------------------------------------
# Each draw: ONE gauge realization, shared by every algorithm, and one realization of each algorithm,
# all drawn per observation from their bands. Stats are computed against the gauge draw (no bands), so
# Miss_bias = (over - under) / N and Miss_magnitude = mean |ln(algo / gauge)|. The algorithms are then
# ranked on each stat within the draw (exact ties only) and only the ranks are kept: close scores that
# would be ties in a single comparison sort themselves out across the draws.

# Split-normal draw in linear space: below the centre the sd is (centre - lower) / z, above it is
# (upper - centre) / z. z converts the band edge to 1 sd (1.96 for a 95% band, 1 for a +/-1 sd band).
# Any draw at or below zero is set to 0.01.
draw_linear <- function(centre, lower, upper, z) {
  e <- rnorm(length(centre))
  d <- centre + ifelse(e < 0, (centre - lower) / z, (upper - centre) / z) * e
  ifelse(d <= 0, 0.01, d)
}


# Per-draw ranks for every reach x algorithm x stat. Columns: reach_id, USGS_ID, algo, type, draw,
# metric, rank. Algorithms not in joined_df are skipped; an algorithm needs at least 20 observations
# with truth, value and both CI edges at a reach. gauge_z / algo_z: band-edge z for gauge / algorithm bands.
# Parallel over reaches; each reach uses seed + its index, so results don't depend on n_cores.
mc_validate <- function(joined_df, algos, central_col, lower_col, upper_col, xwalk, truth_type,
                        n_draw, gauge_z, algo_z, n_cores, seed) {
  algos   <- algos[algos %in% names(joined_df)]
  reaches <- split(joined_df, joined_df$reach_id)

  mc_reach <- function(i) {
    set.seed(seed + i)
    reach_df <- reaches[[i]]

    gauge_ok <- !is.na(reach_df[[central_col]]) & !is.na(reach_df[[lower_col]]) & !is.na(reach_df[[upper_col]])
    algo_ok  <- lapply(setNames(algos, algos), function(a)
      gauge_ok & !is.na(reach_df[[a]]) & !is.na(reach_df[[paste0(a, "_lower_CI")]]) & !is.na(reach_df[[paste0(a, "_upper_CI")]]))
    algo_ok  <- algo_ok[vapply(algo_ok, sum, numeric(1)) >= 20]
    if (length(algo_ok) == 0) return(NULL)

    bind_rows(lapply(seq_len(n_draw), function(k) {
      # one gauge draw for this draw, shared by all algorithms
      g <- draw_linear(reach_df[[central_col]], reach_df[[lower_col]], reach_df[[upper_col]], gauge_z)

      bind_rows(lapply(names(algo_ok), function(a) {
        ok <- algo_ok[[a]]
        m  <- draw_linear(reach_df[[a]][ok], reach_df[[paste0(a, "_lower_CI")]][ok], reach_df[[paste0(a, "_upper_CI")]][ok], algo_z)
        na <- rep(NA_real_, sum(ok))
        # USGS_ID placeholder: rank_algos expects the column; it is filled from xwalk at the end
        bind_cols(tibble(reach_id = names(reaches)[i], USGS_ID = NA_character_, algo = a, type = truth_type),
                  trad_validate(g[ok], m),
                  batesian_validate(reach_df$time[ok], g[ok], m, g[ok], g[ok], na, na) %>% select(Miss_bias, Miss_magnitude))
      })) %>%
        rank_algos(tie_tol = 0) %>%
        transmute(reach_id, algo, type, draw = k, metric, rank)
    }))
  }

  bind_rows(parallel::mclapply(seq_along(reaches), mc_reach, mc.cores = n_cores)) %>%
    left_join(xwalk, by = "reach_id", relationship = "many-to-one") %>%
    relocate(USGS_ID, .after = reach_id)
}