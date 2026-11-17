# Rank figures. Needs algo_key (defined once, in the notebook) and the rank dataframes.
# Defines functions only; call make_rank_figs() from the notebook.

plot_rank_group <- function(rank_df, trad_metrics, bates_metrics, title, formulae) {
  metrics <- c(trad_metrics, bates_metrics)

  # number of algorithms ranked together in this figure, not the length of algo_key
  n_algo    <- n_distinct(rank_df$algo)
  rank_cols <- setNames(rev(viridisLite::viridis(n_algo)), seq_len(n_algo))   # rank 1 = brightest

  plot_df <- rank_df %>%
    filter(metric %in% metrics, !is.na(rank)) %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one") %>%
    mutate(metric = factor(metric, levels = metrics),
           rank   = factor(rank, levels = seq_len(n_algo)))

  n_gauges   <- n_distinct(plot_df$USGS_ID)
  truth_type <- paste(unique(plot_df$type), collapse = ", ")

  ggplot(plot_df, aes(x = metric, fill = rank)) +
    geom_bar(position = "fill") +
    geom_vline(xintercept = length(trad_metrics) + 0.5, linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ model, nrow = 1, labeller = labeller(model = function(x) paste("Model", x))) +
    scale_fill_manual(values = rank_cols, drop = FALSE, name = "Rank") +
    labs(title = paste0(title, ", ", truth_type, " (n = ", n_gauges, " gauges)"),
         x = NULL, y = "Fraction of reaches", caption = formulae) +
    theme_bw(base_size = 22) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.caption = element_text(hjust = 0, family = "mono", size = 16, lineheight = 1.2),
          plot.caption.position = "plot")
}


# Rank heatmap: one panel per model, metrics as rows (traditional above the line, Batesian below),
# ranks as columns, cell = share of that model's ranked gauges on that metric (%).
plot_rank_heatmap <- function(rank_df, trad_metrics, bates_metrics, ncol = 7) {
  metrics <- c(trad_metrics, bates_metrics)
  n_algo  <- n_distinct(rank_df$algo)

  plot_df <- rank_df %>%
    filter(metric %in% metrics, !is.na(rank)) %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one") %>%
    count(model, metric, rank) %>%
    complete(model, metric, rank = seq_len(n_algo), fill = list(n = 0)) %>%
    group_by(model, metric) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(metric = factor(metric, levels = rev(metrics)))

  n_gauges   <- n_distinct(rank_df$USGS_ID)
  truth_type <- paste(unique(rank_df$type), collapse = ", ")

  ggplot(plot_df, aes(x = factor(rank), y = metric, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    # numbers only when cells are wide enough to read them
    { if (n_algo <= 8) geom_text(aes(label = ifelse(pct > 0, round(pct), "")), size = 4) } +
    geom_hline(yintercept = length(bates_metrics) + 0.5, colour = "#E07B00", linewidth = 1) +
    facet_wrap(~ model, ncol = ncol, scales = "free_x",   # free_x: rank labels under every panel
               labeller = labeller(model = function(x) paste("Model", x))) +
    scale_x_discrete(breaks = function(x) x[as.integer(x) %% 2 == 1 | n_algo <= 8]) +
    scale_fill_gradient(low = "white", high = "#08306b", name = "Share of\ngauges (%)") +
    labs(title = paste0("Rank distribution by metric, ", truth_type, " (n = ", n_gauges, " gauges)"),
         x = "Rank", y = NULL) +
    theme_minimal(base_size = 18) +
    theme(panel.grid = element_blank(), panel.spacing = unit(1.2, "lines"),
          strip.text = element_text(face = "bold"))
}


# Rank shift: one panel per model, rank on metric xrank (x) against rank on metric yrank (y),
# cell = share of all gauges (%). Diagonal outlined: on it the two metrics agree. Header gives
# the share of gauges where the yrank rank is better (improved), equal, or worse (declined).
plot_rank_shift <- function(rank_df, xrank, yrank) {
  n_algo <- n_distinct(rank_df$algo)

  paired <- rank_df %>%
    filter(metric %in% c(xrank, yrank)) %>%
    select(reach_id, USGS_ID, algo, metric, rank) %>%
    pivot_wider(names_from = metric, values_from = rank) %>%
    rename(x_rank = all_of(xrank), y_rank = all_of(yrank)) %>%
    filter(!is.na(x_rank), !is.na(y_rank)) %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one")

  n_gauges <- n_distinct(paired$USGS_ID)

  shift <- paired %>%
    group_by(model) %>%
    summarise(improved = 100 * mean(y_rank < x_rank),
              same     = 100 * mean(y_rank == x_rank),
              declined = 100 * mean(y_rank > x_rank),
              .groups  = "drop") %>%
    mutate(strip = sprintf("Model %s\nImproved %.0f%% | Same %.0f%% | Declined %.0f%%",
                           model, improved, same, declined))

  plot_df <- paired %>%
    count(model, x_rank, y_rank) %>%
    complete(model, x_rank = seq_len(n_algo), y_rank = seq_len(n_algo), fill = list(n = 0)) %>%
    group_by(model) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    left_join(select(shift, model, strip), by = "model", relationship = "many-to-one")

  diag_df <- distinct(plot_df, strip) %>%
    cross_join(tibble(x_rank = seq_len(n_algo), y_rank = seq_len(n_algo)))

  ggplot(plot_df, aes(x = x_rank, y = y_rank, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_tile(data = diag_df, aes(x = x_rank, y = y_rank), inherit.aes = FALSE,
              fill = NA, colour = "#E07B00", linewidth = 0.9) +
    geom_text(aes(label = ifelse(pct > 0, round(pct), "")), size = 4) +
    facet_wrap(~ strip, ncol = 3) +
    scale_x_continuous(breaks = seq_len(n_algo), expand = c(0, 0)) +
    scale_y_reverse(breaks = seq_len(n_algo), expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "#08306b",
                        name = paste0("Share of all ", n_gauges, " gauges (%)"),
                        guide = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    labs(title = paste0(xrank, " rank vs ", yrank, " rank, ",
                        paste(unique(rank_df$type), collapse = ", ")),
         x = paste(xrank, "rank"), y = paste(yrank, "rank")) +
    theme_minimal(base_size = 16) +
    theme(panel.grid = element_blank(), panel.spacing = unit(0.8, "lines"),
          strip.text = element_text(size = 12),
          legend.position = "bottom",
          legend.key.width = unit(2.5, "lines"),
          plot.margin = margin(4, 4, 4, 4))
}

# Rank shift as a table: per model, % of gauges where the yrank rank is better than (improved),
# equal to (same) or worse than (declined) the xrank rank. Drawn as a booktabs-style table grob.
rank_shift_table <- function(rank_df, xrank, yrank) {
  shift <- rank_df %>%
    filter(metric %in% c(xrank, yrank)) %>%
    select(reach_id, USGS_ID, algo, metric, rank) %>%
    pivot_wider(names_from = metric, values_from = rank) %>%
    rename(x_rank = all_of(xrank), y_rank = all_of(yrank)) %>%
    filter(!is.na(x_rank), !is.na(y_rank)) %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one") %>%
    group_by(model) %>%
    summarise(n_gauges = n_distinct(USGS_ID),
              improved = 100 * mean(y_rank < x_rank),
              same     = 100 * mean(y_rank == x_rank),
              declined = 100 * mean(y_rank > x_rank),
              .groups  = "drop") %>%
    arrange(model)

  tab <- data.frame(Model          = shift$model,
                    Gauges         = shift$n_gauges,
                    `Improved (%)` = sprintf("%.0f", shift$improved),
                    `Same (%)`     = sprintf("%.0f", shift$same),
                    `Declined (%)` = sprintf("%.0f", shift$declined),
                    check.names = FALSE)

  tg <- gridExtra::tableGrob(tab, rows = NULL, theme = gridExtra::ttheme_minimal(
    base_size = 14,
    core    = list(fg_params = list(hjust = 0.5, x = 0.5), padding = unit(c(14, 6), "pt")),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5), padding = unit(c(14, 8), "pt"))))

  # booktabs rules: thick top and bottom, thin under the header
  rule <- function(y, lwd) grid::segmentsGrob(x0 = 0, x1 = 1, y0 = y, y1 = y, gp = grid::gpar(lwd = lwd))
  n_row <- nrow(tg)
  tg <- gtable::gtable_add_grob(tg, rule(1, 2),   t = 1,     l = 1, r = ncol(tg))
  tg <- gtable::gtable_add_grob(tg, rule(0, 1),   t = 1,     l = 1, r = ncol(tg))
  tg <- gtable::gtable_add_grob(tg, rule(0, 2),   t = n_row, l = 1, r = ncol(tg))

  title <- grid::textGrob(paste0("Rank change from ", xrank, " to ", yrank, ", ",
                                 paste(unique(rank_df$type), collapse = ", ")),
                          gp = grid::gpar(fontsize = 16, fontface = "bold"))

  gridExtra::arrangeGrob(title, tg, ncol = 1, heights = grid::unit.c(unit(2.5, "lines"), sum(tg$heights)))
}

band_def <- "IP: model inside gauge band [LowerCI, UpperCI]   OP: above it   UP: below it   N: all pairs\n"

formulae_skill <- paste0(band_def,
  "Accuracy  = IP / (IP + OP + UP)                        (point prediction vs gauge band)\n",
  "Precision = sum(overlap) / sum(model band width)      (model band vs gauge band, summed over the record)")

formulae_bias <- paste0(band_def,
  "Miss_bias      = (OP - UP) / N      0 = balanced, +1 = always over, -1 = always under\n",
  "Miss_magnitude = mean over misses of ln(model / UpperCI) [OP] or ln(LowerCI / model) [UP]")

formulae_corr <- paste0(
  "Bcorr = fraction of steps where |dln(model) - dln(gauge)| <= 1.96*SE\n",
  "SE    = sqrt(s1^2 + s2^2 - 2*rho*s1*s2),   s = (ln UpperCI - ln LowerCI) / 3.92\n",
  "rho = 0 (Bcorr_rho0) and rho = 0.8 (Bcorr_rho08)")

# Draw every rank figure for each rank dataframe (e.g. list(sos_rank_A15, sos_rank_A24)).
# xrank / yrank: the two metrics compared in the rank-shift figure (x axis, y axis).
make_rank_figs <- function(rank_dfs, xrank, yrank) {
  for (rank_df in rank_dfs) {

    # print(plot_rank_group(rank_df,
    #                       trad_metrics  = c("KGE", "NSE", "Rsq"),
    #                       bates_metrics = c("Accuracy", "Precision"),
    #                       title = "Overall skill",
    #                       formulae = formulae_skill))
    #
    # print(plot_rank_group(rank_df,
    #                       trad_metrics  = c("nBIAS", "abs_nBIAS"),
    #                       bates_metrics = c("Miss_bias", "Miss_magnitude"),
    #                       title = "Bias",
    #                       formulae = formulae_bias))
    #
    # print(plot_rank_group(rank_df,
    #                       trad_metrics  = "r",
    #                       bates_metrics = c("Bcorr_rho0", "Bcorr_rho08"),
    #                       title = "Correlation",
    #                       formulae = formulae_corr))

    print(plot_rank_heatmap(rank_df,
                            trad_metrics  = c("NSE", "KGE", "abs_nBIAS"),
                            bates_metrics = c("Accuracy", "Miss_bias", "Miss_magnitude", "Precision")))

    # more than 8 ranks (e.g. A24, all 14 algorithms): the grid is unreadable, so print the table
    if (n_distinct(rank_df$algo) > 8) {
      grid::grid.newpage()
      grid::grid.draw(rank_shift_table(rank_df, xrank = xrank, yrank = yrank))
    } else {
      print(plot_rank_shift(rank_df, xrank = xrank, yrank = yrank))
    }
  }
}


# Percent change from A15 to A24, per metric and model: 100 * (A24 - A15) / |A15|, computed per
# gauge on gauges present in both, then summarised by the median across gauges.
plot_truth_change <- function(val_from, val_to, metrics) {
  long <- function(v) v %>%
    select(reach_id, USGS_ID, algo, any_of(metrics)) %>%
    pivot_longer(any_of(metrics), names_to = "metric", values_to = "value")

  from_type <- unique(val_from$type)
  to_type   <- unique(val_to$type)

  change <- inner_join(long(val_from), long(val_to),
                       by = c("reach_id", "USGS_ID", "algo", "metric"),
                       suffix = c("_from", "_to"), relationship = "one-to-one") %>%
    filter(!is.na(value_from), !is.na(value_to), value_from != 0) %>%
    mutate(pct  = 100 * (value_to - value_from) / abs(value_from),
           diff = value_to - value_from) %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one") %>%
    group_by(model, metric) %>%
    summarise(median_pct = median(pct), min_diff = min(diff), med_diff = median(diff), max_diff = max(diff),
              n_gauges = n(), .groups = "drop") %>%
    mutate(metric = factor(metric, levels = rev(metrics)))

  lim <- max(abs(change$median_pct))

  ggplot(change, aes(x = model, y = metric, fill = median_pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("max %+.2f\nmed %+.2f\nmin %+.2f", max_diff, med_diff, min_diff)),
              size = 3.6, lineheight = 0.9) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0,
                         limits = c(-lim, lim), name = "Median\nchange (%)") +
    scale_x_discrete(labels = function(x) paste("Model", x)) +
    labs(title = paste0("Change in each metric from ", from_type, " to ", to_type),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(panel.grid = element_blank())
}


# One model at one gauge against the gauge record, on a compressed time axis: only days with
# model observations are drawn, side by side, each as one slot (midnight to midnight). Within a slot:
# the 15-minute gauge record (line) and its band (light ribbon), the daily mean (black) and its daily
# band (darker), the model at overpass time (small dot, thin bar) and its daily mean (open circle at
# noon, thick faint bar). Table: stats at 15 minute and daily, and their change (daily - 15 minute).
plot_model_A15_A24 <- function(model_letter, gauge, USGS_df, usgs_A24,
                               sos_USGS_A15, sos_USGS_A24, sos_val_A15, sos_val_A24) {
  model_algo <- algo_key$algo[algo_key$model == model_letter]
  model_col  <- "black"
  gap        <- 0.25                      # space between day slots, in days

  model_pts <- function(df, sfx) df %>%
    filter(USGS_ID == gauge) %>%
    transmute(time     = if (sfx == "A24") as.POSIXct(time, tz = "UTC") + 12 * 3600 else as.POSIXct(time, tz = "UTC"),
              Q_truth  = .data[[paste0("Q_", sfx)]],
              Q_lower  = .data[[paste0("Q_lower_", sfx)]],
              Q_upper  = .data[[paste0("Q_upper_", sfx)]],
              Q_model  = .data[[model_algo]],
              Q_mlower = .data[[paste0(model_algo, "_lower_CI")]],
              Q_mupper = .data[[paste0(model_algo, "_upper_CI")]],
              type     = sfx) %>%
    filter(!is.na(Q_model), !is.na(Q_truth))

  m15 <- model_pts(sos_USGS_A15, "A15")
  m24 <- model_pts(sos_USGS_A24, "A24")

  # compressed axis: each observation day becomes a slot; x = slot start + fraction of the day
  obs_days <- sort(unique(as.Date(c(m15$time, m24$time), tz = "UTC")))
  slot_x   <- function(time) {
    day <- as.Date(time, tz = "UTC")
    (match(day, obs_days) - 1) * (1 + gap) + as.numeric(difftime(time, as.POSIXct(day, tz = "UTC"), units = "days"))
  }

  m15 <- m15 %>% mutate(x = slot_x(time))
  m24 <- m24 %>% mutate(x = slot_x(time))

  gauge15 <- USGS_df %>%
    filter(USGS_ID == gauge, !is.na(Q_approved), as.Date(time, tz = "UTC") %in% obs_days) %>%
    select(time, Q_approved, Q_lower, Q_upper) %>%
    mutate(x = slot_x(time), slot = match(as.Date(time, tz = "UTC"), obs_days))

  gauge24 <- usgs_A24 %>%
    filter(USGS_ID == gauge, date %in% obs_days) %>%
    mutate(x0 = (match(date, obs_days) - 1) * (1 + gap), x1 = x0 + 1)

  y_top    <- max(c(gauge15$Q_upper, gauge24$Q_upper_A24, m15$Q_mupper, m24$Q_mupper,
                    m15$Q_model, m24$Q_model), na.rm = TRUE)
  y_bottom <- min(c(gauge15$Q_lower, gauge24$Q_lower_A24, m15$Q_mlower, m24$Q_mlower), na.rm = TRUE)

  # date labels at slot centres, thinned so they don't collide
  label_every <- max(1, ceiling(length(obs_days) / 12))
  lab_idx     <- seq(1, length(obs_days), by = label_every)
  breaks      <- (lab_idx - 1) * (1 + gap) + 0.5
  labels      <- format(obs_days[lab_idx], "%Y-%m-%d")
  separators  <- (seq_along(obs_days)[-1] - 1) * (1 + gap) - gap / 2

  counts <- bind_rows(m15, m24) %>%
    group_by(type) %>%
    summarise(IP = sum(Q_model >= Q_lower & Q_model <= Q_upper),
              OP = sum(Q_model > Q_upper),
              UP = sum(Q_model < Q_lower),
              .groups = "drop")

  # stats table: rows = stats, columns = 15 minute, daily, change (daily - 15 minute, absolute)
  stat_names <- c("KGE", "NSE", "r", "nBIAS", "Accuracy", "Precision",
                  "Miss_bias", "Miss_magnitude", "Bcorr_rho0", "IP", "OP", "UP")
  stats <- bind_rows(sos_val_A15, sos_val_A24) %>%
    filter(algo == model_algo, USGS_ID == gauge) %>%
    left_join(counts, by = "type", relationship = "one-to-one") %>%
    select(type, all_of(stat_names)) %>%
    pivot_longer(-type, names_to = "stat", values_to = "value") %>%
    pivot_wider(names_from = type, values_from = value) %>%
    mutate(stat = factor(stat, levels = stat_names)) %>%
    arrange(stat)

  is_count <- stats$stat %in% c("IP", "OP", "UP")
  tab <- data.frame(Stat        = as.character(stats$stat),
                    `15 minute` = ifelse(is_count, sprintf("%.0f", stats$A15), sprintf("%.2f", stats$A15)),
                    Daily       = ifelse(is_count, sprintf("%.0f", stats$A24), sprintf("%.2f", stats$A24)),
                    Change      = ifelse(is_count, sprintf("%+.0f", stats$A24 - stats$A15),
                                                   sprintf("%+.2f", stats$A24 - stats$A15)),
                    check.names = FALSE)
  tab_grob <- gridExtra::tableGrob(tab, rows = NULL,
                                   theme = gridExtra::ttheme_minimal(base_size = 13,
                                     core = list(fg_params = list(hjust = 1, x = 0.95))))

  hydro <- ggplot() +
    geom_vline(xintercept = separators, colour = "grey85", linewidth = 0.3) +
    geom_ribbon(data = gauge15, aes(x, ymin = Q_lower, ymax = Q_upper, group = slot,
                                    fill = "15 minute gauge band")) +
    geom_rect(data = gauge24, aes(xmin = x0, xmax = x1, ymin = Q_lower_A24, ymax = Q_upper_A24,
                                  fill = "Daily gauge band")) +
    geom_line(data = gauge15, aes(x, Q_approved, group = slot, colour = "15 minute gauge"), linewidth = 0.4) +
    geom_segment(data = gauge24, aes(x = x0, xend = x1, y = Q_A24, yend = Q_A24, colour = "Daily gauge mean"),
                 linewidth = 0.8) +
    geom_linerange(data = m24, aes(x, ymin = Q_mlower, ymax = Q_mupper), colour = model_col,
                   alpha = 0.3, linewidth = 2.4) +
    geom_linerange(data = m15, aes(x, ymin = Q_mlower, ymax = Q_mupper), colour = model_col,
                   linewidth = 0.6) +
    geom_point(data = m24, aes(x, Q_model, shape = "Model, daily mean"), colour = model_col,
               size = 4, stroke = 1.2) +
    geom_point(data = m15, aes(x, Q_model, shape = "Model, at overpass"), colour = model_col, size = 1.8) +
    scale_fill_manual(values = c("15 minute gauge band" = alpha("#6baed6", 0.45), "Daily gauge band" = alpha("#de2d26", 0.3)),
                      name = NULL) +
    scale_colour_manual(values = c("15 minute gauge" = "#2171b5", "Daily gauge mean" = "#a50f15"), name = NULL) +
    scale_shape_manual(values = c("Model, at overpass" = 16, "Model, daily mean" = 1), name = NULL) +
    scale_x_continuous(breaks = breaks, labels = labels, expand = expansion(add = 0.3)) +
    guides(fill   = guide_legend(order = 1),
           colour = guide_legend(order = 2, override.aes = list(linewidth = 1.2)),
           shape  = guide_legend(order = 3, override.aes = list(colour = model_col, size = 3.5))) +
    coord_cartesian(ylim = c(y_bottom, y_top)) +
    labs(title = paste0("Model ", model_letter, " at ", gauge, ": 15 minute and daily"),
         x = NULL, y = expression(Q~(m^3/s))) +
    theme_bw(base_size = 16) +
    theme(legend.box = "vertical", legend.text = element_text(size = 13),
          legend.spacing.y = unit(0.2, "lines"),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  # right-hand column: the plot's legend above the table
  hydro_grob <- ggplotGrob(hydro + theme(legend.position = "right"))
  legend     <- hydro_grob$grobs[[grep("guide-box", hydro_grob$layout$name)[1]]]

  (hydro + theme(legend.position = "none")) +
    (patchwork::wrap_elements(full = legend) / patchwork::wrap_elements(full = tab_grob) +
       patchwork::plot_layout(heights = c(1, 1.5))) +
    patchwork::plot_layout(widths = c(2.6, 1))
}


# Monte Carlo: rows = stats, columns = models, cell = % of draws in which the model ranks first, averaged
# over gauges (every gauge has the same number of draws). Exact ties share first, so a row can exceed 100%.
# mc_ranks: mc_validate() output for one type.
plot_mc_best <- function(mc_ranks, metrics = c("KGE", "NSE", "r", "nBIAS", "Miss_bias", "Miss_magnitude")) {
  n_gauges <- n_distinct(mc_ranks$USGS_ID)
  lv       <- algo_key$model[algo_key$algo %in% unique(mc_ranks$algo)]

  plot_df <- mc_ranks %>%
    filter(metric %in% metrics) %>%
    group_by(USGS_ID, metric, algo) %>%
    summarise(first = mean(rank == 1), .groups = "drop") %>%
    group_by(metric, algo) %>%
    summarise(pct = 100 * mean(first), .groups = "drop") %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one") %>%
    mutate(model  = factor(model, levels = lv),
           metric = factor(metric, levels = rev(metrics)))

  ggplot(plot_df, aes(model, metric, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", pct), colour = ifelse(pct > 40, "white", "black")), size = 3.8) +
    scale_fill_distiller(palette = "Greens", direction = 1, name = "Draws ranked\nfirst (%)") +
    scale_colour_identity() +
    labs(title = paste0("Share of draws each model ranks first, ", paste(unique(mc_ranks$type), collapse = ", "),
                        " (n = ", n_gauges, " gauges, ", max(mc_ranks$draw), " draws each)"),
         x = "Model", y = NULL) +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank())
}


# compressed observation axis: one slot per observation day, labels thinned
obs_axis <- function(days) {
  every <- max(1, ceiling(length(days) / 12)); i <- seq(1, length(days), by = every)
  list(breaks = i, labels = format(days[i], "%Y-%m-%d"))
}

# MC worked example, same scenario as plot_nse_vs_precision: the same gauge, observations and two models
# (rank 1 on NSE and rank 1 on Precision in rank_df). At each observation, clouds of gauge draws and of
# each model's draws. Table: per stat, % of draws each model ranks first (mean rank) among all models ranked.
plot_mc_worked <- function(gauge, sfx, joined_df, rank_df, mc_ranks, gauge_z, algo_z, n_show = 40, seed = 1) {
  r1 <- rank_df %>% filter(USGS_ID == gauge, rank == 1)
  a_nse  <- r1$algo[r1$metric == "NSE"]
  a_prec <- r1$algo[r1$metric == "Precision"]
  stopifnot(length(a_nse) == 1, length(a_prec) == 1, a_nse != a_prec)   # no ties at rank 1
  lab  <- function(x) paste("Model", algo_key$model[algo_key$algo == x])
  who_labels <- c(Gauge = "Gauge draws", nse = paste(lab(a_nse), "draws (best NSE)"), prec = paste(lab(a_prec), "draws (best Precision)"))
  cols <- setNames(c("black", "#2166ac", "#b2182b"), who_labels)

  df   <- joined_df %>% filter(USGS_ID == gauge) %>% arrange(time)
  ax   <- obs_axis(as.Date(df$time, tz = "UTC"))
  slot <- seq_len(nrow(df))
  set.seed(seed)
  clouds <- bind_rows(lapply(seq_len(n_show), function(k) bind_rows(
    tibble(slot = slot - 0.25, who = who_labels[["Gauge"]],
           Q = draw_linear(df[[paste0("Q_central_", sfx)]], df[[paste0("Q_lower_", sfx)]], df[[paste0("Q_upper_", sfx)]], gauge_z)),
    tibble(slot = slot,        who = who_labels[["nse"]],
           Q = draw_linear(df[[a_nse]], df[[paste0(a_nse, "_lower_CI")]], df[[paste0(a_nse, "_upper_CI")]], algo_z)),
    tibble(slot = slot + 0.25, who = who_labels[["prec"]],
           Q = draw_linear(df[[a_prec]], df[[paste0(a_prec, "_lower_CI")]], df[[paste0(a_prec, "_upper_CI")]], algo_z))))) %>%
    mutate(source = factor(who, levels = unname(who_labels)))

  hydro <- ggplot(clouds, aes(slot, Q, colour = source)) +
    geom_point(position = position_jitter(width = 0.05, height = 0, seed = 1), alpha = 0.25, size = 1.1, stroke = 0) +
    scale_colour_manual(values = cols, name = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    scale_x_continuous(breaks = ax$breaks, labels = ax$labels, expand = expansion(add = 0.6)) +
    labs(title = paste0(gauge, ", ", sfx, ": Monte Carlo, best NSE vs best Precision"),
         x = NULL, y = expression(Q~(m^3/s))) +
    theme_bw(base_size = 15) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top",
          legend.text = element_text(size = 13))

  stat_names <- c("KGE", "NSE", "r", "nBIAS", "Miss_bias", "Miss_magnitude")
  summ <- mc_ranks %>%
    filter(USGS_ID == gauge, metric %in% stat_names, algo %in% c(a_nse, a_prec)) %>%
    group_by(algo, metric) %>%
    summarise(cell = sprintf("%.0f%% (%.1f)", 100 * mean(rank == 1), mean(rank)), .groups = "drop") %>%
    pivot_wider(names_from = algo, values_from = cell) %>%
    mutate(metric = factor(metric, levels = stat_names)) %>% arrange(metric)
  tab <- data.frame(Stat = as.character(summ$metric), A = summ[[a_nse]], B = summ[[a_prec]], check.names = FALSE)
  names(tab)[2:3] <- c(lab(a_nse), lab(a_prec))
  tg <- gridExtra::tableGrob(tab, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 13,
          core = list(fg_params = list(hjust = 1, x = 0.95)), colhead = list(fg_params = list(fontface = "bold"))))

  patchwork::wrap_plots(hydro, patchwork::wrap_elements(full = tg), widths = c(3, 1.1))
}


# Gauges where one model alone is best on NSE and a different model alone is best on Precision (no ties
# at rank 1). Sorted clearest first: the NSE winner's Precision rank plus the Precision winner's NSE rank.
find_nse_precision_gauges <- function(rank_df) {
  nse_prec <- rank_df %>%
    filter(metric %in% c("NSE", "Precision")) %>%
    select(reach_id, USGS_ID, type, algo, metric, rank)

  sole_winner <- function(m) nse_prec %>%
    filter(metric == m, rank == 1) %>%
    group_by(reach_id, USGS_ID, type) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    select(reach_id, USGS_ID, type, algo)

  inner_join(rename(sole_winner("NSE"), NSE = algo), rename(sole_winner("Precision"), Precision = algo),
             by = c("reach_id", "USGS_ID", "type"), relationship = "one-to-one") %>%
    filter(NSE != Precision) %>%
    left_join(nse_prec %>% filter(metric == "Precision") %>% select(reach_id, type, algo, nse_best_precision_rank = rank),
              by = c("reach_id", "type", "NSE" = "algo"), relationship = "one-to-one") %>%
    left_join(nse_prec %>% filter(metric == "NSE") %>% select(reach_id, type, algo, precision_best_nse_rank = rank),
              by = c("reach_id", "type", "Precision" = "algo"), relationship = "one-to-one") %>%
    arrange(desc(nse_best_precision_rank + precision_best_nse_rank))
}

# Hydrograph of one gauge on a compressed observation axis: gauge band (grey bar) and value, the best-NSE
# model and the best-Precision model, each with its band. Table: both models' NSE and Precision and ranks.
plot_nse_vs_precision <- function(gauge, sfx, joined_df, val_df, rank_df) {
  r1 <- rank_df %>% filter(USGS_ID == gauge, rank == 1)
  a_nse  <- r1$algo[r1$metric == "NSE"]
  a_prec <- r1$algo[r1$metric == "Precision"]
  stopifnot(length(a_nse) == 1, length(a_prec) == 1, a_nse != a_prec)   # no ties at rank 1
  lab <- function(x) paste("Model", algo_key$model[algo_key$algo == x])
  labs2 <- c(lab(a_nse), lab(a_prec)); names(labs2) <- c("nse", "prec")
  cols <- setNames(c("#2166ac", "#b2182b"), paste0(labs2, c(" (best NSE)", " (best Precision)")))

  df <- joined_df %>% filter(USGS_ID == gauge) %>% arrange(time) %>% mutate(slot = row_number())
  ax <- obs_axis(as.Date(df$time, tz = "UTC"))
  mods <- bind_rows(
    tibble(slot = df$slot - 0.18, Q = df[[a_nse]],  lo = df[[paste0(a_nse,  "_lower_CI")]], up = df[[paste0(a_nse,  "_upper_CI")]], who = names(cols)[1]),
    tibble(slot = df$slot + 0.18, Q = df[[a_prec]], lo = df[[paste0(a_prec, "_lower_CI")]], up = df[[paste0(a_prec, "_upper_CI")]], who = names(cols)[2]))

  hydro <- ggplot() +
    geom_rect(data = df, aes(xmin = slot - 0.4, xmax = slot + 0.4, ymin = .data[[paste0("Q_lower_", sfx)]],
                             ymax = .data[[paste0("Q_upper_", sfx)]], fill = "Gauge")) +
    geom_segment(data = df, aes(x = slot - 0.4, xend = slot + 0.4, y = .data[[paste0("Q_central_", sfx)]],
                                yend = .data[[paste0("Q_central_", sfx)]], linetype = "Gauge"), linewidth = 0.8) +
    geom_linerange(data = mods, aes(slot, ymin = lo, ymax = up, colour = who), linewidth = 0.9) +
    geom_point(data = mods, aes(slot, Q, colour = who), size = 2.2) +
    scale_colour_manual(values = cols, name = NULL) +
    # fill and linetype share the one break "Gauge", so they merge into one key: line inside the band
    scale_fill_manual(values = c(Gauge = "grey80"), name = NULL) +
    scale_linetype_manual(values = c(Gauge = "solid"), name = NULL) +
    guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    scale_x_continuous(breaks = ax$breaks, labels = ax$labels, expand = expansion(add = 0.6)) +
    labs(x = NULL, y = expression(Q~(m^3/s)),
         title = paste0(gauge, ", ", sfx, ": best NSE is not best Precision")) +
    theme_bw(base_size = 15) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

  st <- c("NSE", "KGE", "Accuracy", "Precision")
  get <- function(x, m) val_df[[m]][val_df$USGS_ID == gauge & val_df$algo == x]
  rk  <- function(x, m) rank_df$rank[rank_df$USGS_ID == gauge & rank_df$algo == x & rank_df$metric == m]
  tab <- data.frame(Stat = st,
                    A = sapply(st, function(m) sprintf("%.2f (rank %d)", get(a_nse, m), rk(a_nse, m))),
                    B = sapply(st, function(m) sprintf("%.2f (rank %d)", get(a_prec, m), rk(a_prec, m))), check.names = FALSE)
  names(tab)[2:3] <- labs2
  tg <- gridExtra::tableGrob(tab, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 13,
          core = list(fg_params = list(hjust = 1, x = 0.95)), colhead = list(fg_params = list(fontface = "bold"))))
  patchwork::wrap_plots(hydro, patchwork::wrap_elements(full = tg), widths = c(3, 1.1))
}


# Single-comparison counterpart of plot_mc_best: rows = stats, columns = models, cell = % of gauges where
# the model ranks first (ties count for every tied model, so a row can exceed 100%). Same colours as the
# rank heatmaps. rank_df: rank_algos() output for one type (e.g. sos_rank_A15).
plot_rank_first <- function(rank_df, metrics = c("KGE", "NSE", "r", "nBIAS", "Miss_bias", "Miss_magnitude",
                                                  "Accuracy", "Precision")) {
  lv <- algo_key$model[algo_key$algo %in% unique(rank_df$algo)]

  ranked  <- rank_df %>% filter(metric %in% metrics, !is.na(rank))
  plot_df <- ranked %>%
    group_by(metric) %>%
    mutate(n_gauges = n_distinct(USGS_ID)) %>%
    group_by(metric, algo, n_gauges) %>%
    summarise(pct = 100 * sum(rank == 1) / first(n_gauges), .groups = "drop") %>%
    left_join(algo_key, by = "algo", relationship = "many-to-one") %>%
    mutate(model  = factor(model, levels = lv),
           metric = factor(metric, levels = rev(metrics)))

  ggplot(plot_df, aes(model, metric, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", pct), colour = ifelse(pct > 40, "white", "black")), size = 3.8) +
    # line above the stats that have no Monte Carlo counterpart
    geom_hline(yintercept = sum(metrics %in% c("Accuracy", "Precision")) + 0.5, linewidth = 0.8) +
    scale_fill_gradient(low = "white", high = "#08306b", name = "Gauges ranked\nfirst (%)") +
    scale_colour_identity() +
    labs(title = paste0("Share of gauges each model ranks first, ", paste(unique(rank_df$type), collapse = ", "),
                        " (n = ", n_distinct(ranked$USGS_ID), " gauges)"),
         x = "Model", y = NULL) +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank())
}