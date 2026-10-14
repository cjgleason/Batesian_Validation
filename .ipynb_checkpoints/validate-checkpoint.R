validate_algo <- function(algo_name, joined_df, truth_col, lower_col, upper_col, xwalk, val_type, truth_type) {
  val_type <- match.arg(val_type, c("trad", "Batesian"))

  algo_df <- joined_df %>%
    select(reach_id,
           time,
           Q_truth = all_of(truth_col),
           Q_lower = all_of(lower_col),
           Q_upper = all_of(upper_col),
           Q_algo  = all_of(algo_name))

  per_reach <- lapply(split(algo_df, algo_df$reach_id), function(reach_df) {
    switch(val_type,
      trad     = trad_validate(truth = reach_df$Q_truth, prediction = reach_df$Q_algo),
      Batesian = batesian_validate(time       = reach_df$time,
                                   truth      = reach_df$Q_truth,
                                   prediction = reach_df$Q_algo,
                                   lower      = reach_df$Q_lower,
                                   upper      = reach_df$Q_upper))
  })

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



batesian_validate <- function(time, truth, prediction, lower, upper) {
  ok         <- !is.na(truth) & !is.na(prediction) & !is.na(lower) & !is.na(upper)
  time       <- time[ok]
  truth      <- truth[ok]
  prediction <- prediction[ok]
  lower      <- lower[ok]
  upper      <- upper[ok]
  if (length(prediction) < 20) return(NULL)

  is_OP <- prediction > upper
  is_UP <- prediction < lower
  is_IP <- prediction >= lower & prediction <= upper
  miss  <- is_OP | is_UP

  # Accuracy  = IP / (IP + OP + UP), equivalent to IoU / CSI / threat score.
  #   Precision, Recall and F1 were dropped: a point prediction has no extent, so every miss is
  #   both a false positive and a false negative, and true precision = recall = Accuracy.
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
  )
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



rank_algos <- function(val_df) {
  higher_better <- c("r", "Rsq", "NSE", "KGE", "Accuracy", "Bcorr_rho0", "Bcorr_rho08")
  lower_better  <- c("RMSE", "nRMSE", "abs_nBIAS", "Miss_magnitude")
  target_better <- c(nBIAS = 0, Miss_bias = 0)
  # these tie on an absolute difference of 0.025; all others on a 2.5% relative difference
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
    mutate(rank = tie_rank(score, absolute = first(metric) %in% absolute_tie)) %>%
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