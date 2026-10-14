# fixed anonymization: never changes with a seed
algo_key <- tibble(algo  = c("metroman", "sad", "consensus", "momma", "busboi", "sic4dvar"),
                   model = c("A",        "B",   "C",         "D",     "E",      "F"))

n_algo    <- nrow(algo_key)
rank_cols <- setNames(rev(viridisLite::viridis(n_algo)), seq_len(n_algo))   # rank 1 = brightest

plot_rank_group <- function(rank_df, trad_metrics, bates_metrics, title, formulae) {
  metrics <- c(trad_metrics, bates_metrics)

  plot_df <- rank_df %>%
    filter(metric %in% metrics, !is.na(rank)) %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one") %>%
    mutate(metric = factor(metric, levels = metrics),
           rank   = factor(rank, levels = seq_len(n_algo)))

  n_gauges <- n_distinct(plot_df$USGS_ID)

  ggplot(plot_df, aes(x = metric, fill = rank)) +
    geom_bar(position = "fill") +
    geom_vline(xintercept = length(trad_metrics) + 0.5, linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ model, nrow = 1, labeller = labeller(model = function(x) paste("Model", x))) +
    scale_fill_manual(values = rank_cols, drop = FALSE, name = "Rank") +
    labs(title = paste0(title, " (n = ", n_gauges, " gauges)"),
         x = NULL, y = "Fraction of reaches", caption = formulae) +
    theme_bw(base_size = 22) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.caption = element_text(hjust = 0, family = "mono", size = 16, lineheight = 1.2),
          plot.caption.position = "plot")
}

band_def <- "IP: model inside gauge 95% band [LowerCI, UpperCI]   OP: above it   UP: below it   N: all pairs\n"

plot_rank_group(sos_rank_approvedQ_15min,
                trad_metrics  = c("KGE", "NSE", "Rsq"),
                bates_metrics = "Accuracy",
                title = "Overall skill",
                formulae = paste0(band_def,
                  "Accuracy = IP / (IP + OP + UP)"))

plot_rank_group(sos_rank_approvedQ_15min,
                trad_metrics  = c("nBIAS", "abs_nBIAS"),
                bates_metrics = c("Miss_bias", "Miss_magnitude"),
                title = "Bias",
                formulae = paste0(band_def,
                  "Miss_bias      = (OP - UP) / N      0 = balanced, +1 = always over, -1 = always under\n",
                  "Miss_magnitude = mean over misses of ln(model / UpperCI) [OP] or ln(LowerCI / model) [UP]"))

plot_rank_group(sos_rank_approvedQ_15min,
                trad_metrics  = "r",
                bates_metrics = c("Bcorr_rho0", "Bcorr_rho08"),
                title = "Correlation",
                formulae = paste0(
                  "Bcorr = fraction of overpass steps where |dln(model) - dln(gauge)| <= 1.96*SE\n",
                  "SE    = sqrt(s1^2 + s2^2 - 2*rho*s1*s2),   s = (ln UpperCI - ln LowerCI) / 3.92\n",
                  "rho = 0 (Bcorr_rho0) and rho = 0.8 (Bcorr_rho08)"))