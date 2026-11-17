# Explainer figures: one synthetic gauge, eight model scenarios with every metric printed on every
# panel, and a Monte Carlo example (each draw: both models against the same gauge draw, ranked per stat).
# Defines functions only. Needs validate.R (trad_validate, batesian_validate) sourced first.
#   make_example_data(out_file): builds the synthetic gauge, scenarios and Monte Carlo draws, saves them as an rds
#   plot_example_figs(example_data): returns list(trad = <figure>, bates = <figure>, mc = <figure>)

make_example_data <- function(out_file, seed = 4, mc_seed = 7, curve_seed = 11, n_draw = 500, n_curve = 40) {
  set.seed(seed)

  hydro <- function(d) 100 + 25 * sin(2 * pi * (d - 60) / 365) +
    45 * exp(-((d - 130) / 14)^2) + 30 * exp(-((d - 260) / 18)^2)
  log_range  <- range(log(hydro(1:365)))
  half_width <- function(q) 0.05 + 0.25 * ((log(q) - log_range[1]) / diff(log_range))^1.5   # log half-width, widens with flow

  # gauge: continuous record and its band
  daily <- tibble(day = 1:365, Q = hydro(day),
                  Q_lower = Q * exp(-half_width(Q)), Q_upper = Q * exp(half_width(Q)))

  # overpasses every 6 days
  day   <- seq(3, 363, by = 6)
  n     <- length(day)
  truth <- hydro(day)
  hw    <- half_width(truth)
  upper <- truth * exp(hw)

  no_band <- function(m) list(model = m, model_lower = rep(NA_real_, n), model_upper = rep(NA_real_, n))
  band    <- function(m, w) list(model = m, model_lower = m * exp(-w * hw), model_upper = m * exp(w * hw))

  trad <- list(
    "High KGE, low NSE"           = no_band(truth * 1.15 + rnorm(n, 0, 6)),
    "High NSE, low KGE"           = no_band(mean(truth) + 0.5 * (truth - mean(truth))),
    "High bias, good correlation" = no_band(truth * 1.4 + rnorm(n, 0, 3)),
    "Low bias, poor correlation"  = no_band(mean(truth) + sample(truth - mean(truth)) * 0.9))

  bates <- list(
    "High miss bias, low miss magnitude" = no_band(upper * exp(runif(n, 0.01, 0.04))),
    "Low miss bias, low accuracy"        = no_band(truth * exp(sample(rep(c(-1, 1), length.out = n)) * (hw + runif(n, 0.02, 0.08)))),
    "Good precision"                     = band(truth * exp(rnorm(n, 0, 0.02)), 0.5),
    "Poor precision"                     = band(truth * exp(rnorm(n, 0, 0.02)), 3))

  scenario_df <- function(sc, group) bind_rows(lapply(names(sc), function(nm)
    tibble(group = group, scenario = nm, day = day, Q = truth,
           Q_lower = truth * exp(-hw), Q_upper = upper,
           model = sc[[nm]]$model, model_lower = sc[[nm]]$model_lower, model_upper = sc[[nm]]$model_upper)))

  # Monte Carlo: independent draws of the gauge (rating curve) and of each model ----------------
  sq <- 2 * (log(truth) - log_range[1]) / diff(log_range) - 1                 # scaled stage, -1..1
  # gauge draw: rating-curve error = offset + tilt with stage (shared across the record), sd ~ band / 1.96
  gauge_draw <- function() { z <- rnorm(2); truth * exp(hw / 1.96 * (0.8 * z[1] + 0.6 * z[2] * sq)) }
  # model draw: AR(1) noise (rho 0.8) around the model's central estimate, w x gauge half-width
  ar1 <- function(rho = 0.8) { e <- numeric(n); e[1] <- rnorm(1); for (i in 2:n) e[i] <- rho * e[i - 1] + sqrt(1 - rho^2) * rnorm(1); e }
  model_draw <- function(central, w) central * exp(w * hw / 1.96 * ar1())

  # c: central estimate; w: noise width of the draws (x gauge half-width)
  mc_scen <- list(
    "Model A clearly better"                   = list(A = list(c = truth * 1.03, w = 0.5), B = list(c = truth * 1.25, w = 2.0)),
    # A: unbiased noisy draws (centre scaled by exp(-sd^2 / 2) so the noise is unbiased in linear Q)
    # B: smooth draws offset 5.5% high. Same squared error, so the same NSE; B is almost always over
    "Same NSE, poor miss bias for B"           = list(A = list(c = truth * exp(-(1.0 * hw / 1.96)^2 / 2), w = 1.0),
                                                      B = list(c = truth * exp(0.055),                     w = 0.3)))

  # per draw: ONE gauge draw shared by both models, one draw of each model, stats against the gauge draw
  # (no band: Miss_bias = (over - under) / N, Miss_magnitude = mean |ln(model / gauge)|), then the models
  # are ranked on each stat within the draw. Only the ranks are kept.
  set.seed(mc_seed)
  mc_ranks <- bind_rows(lapply(names(mc_scen), function(nm) {
    sc <- mc_scen[[nm]]
    bind_rows(lapply(seq_len(n_draw), function(k) {
      g <- gauge_draw()
      bind_rows(lapply(c("A", "B"), function(mn) {
        m <- model_draw(sc[[mn]]$c, sc[[mn]]$w)
        bind_cols(tibble(reach_id = nm, USGS_ID = nm, algo = mn, type = "MC"),
                  trad_validate(g, m),
                  batesian_validate(day, g, m, g, g, rep(NA_real_, n), rep(NA_real_, n)) %>% select(Miss_bias, Miss_magnitude))
      })) %>%
        select(reach_id, USGS_ID, algo, type, KGE, NSE, r, nBIAS, Miss_bias, Miss_magnitude) %>%
        rank_algos(tie_tol = 0) %>%                      # exact ties only
        transmute(scenario = reach_id, draw = k, model = algo, metric, rank)
    }))
  }))

  # a smaller set of draws to plot
  mc_curves <- bind_rows(lapply(names(mc_scen), function(nm) {
    sc <- mc_scen[[nm]]
    set.seed(curve_seed)
    bind_rows(
      bind_rows(lapply(seq_len(n_curve), function(i) tibble(scenario = nm, day = day, Q = gauge_draw(), who = "Gauge", id = i))),
      bind_rows(lapply(c("A", "B"), function(mn) bind_rows(lapply(seq_len(n_curve), function(i)
        tibble(scenario = nm, day = day, Q = model_draw(sc[[mn]]$c, sc[[mn]]$w), who = mn, id = i))))))
  }))

  example_data <- list(daily     = daily,
                       overpass  = bind_rows(scenario_df(trad, "trad"), scenario_df(bates, "bates")),
                       mc_ranks  = mc_ranks,
                       mc_curves = mc_curves)
  saveRDS(example_data, out_file)
  example_data
}


plot_example_figs <- function(example_data) {
  daily <- example_data$daily
  op    <- example_data$overpass

  gauge <- list(geom_ribbon(data = daily, aes(day, ymin = Q_lower, ymax = Q_upper, fill = "Gauge")),
                geom_line(data = daily, aes(day, Q, colour = "Gauge"), linewidth = 0.6))
  th <- theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 12, lineheight = 1.1),
          panel.grid.minor = element_blank(), axis.title = element_blank())

  # legend: fill (bands) and colour (lines, points) share the breaks Gauge / Model and the same labels, so
  # they merge into one key per entry: the gauge line inside its band, the model point inside its band
  key_scales <- function(model_label, model_fill) list(
    scale_fill_manual(values = c(Gauge = alpha("#6baed6", 0.45), Model = model_fill), limits = c("Gauge", "Model"),
                      labels = c(Gauge = "Gauge and its band", Model = model_label), name = NULL),
    scale_colour_manual(values = c(Gauge = "#2171b5", Model = "black"), limits = c("Gauge", "Model"),
                        labels = c(Gauge = "Gauge and its band", Model = model_label), name = NULL),
    guides(colour = guide_legend(override.aes = list(linetype = c(1, 0), shape = c(NA, 16)))))
  legend_inside <- theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1),
                         legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
                         legend.text = element_text(size = 11))

  panel <- function(df, show_legend, grp) {
    t <- trad_validate(df$Q, df$model)
    b <- batesian_validate(df$day, df$Q, df$model, df$Q_lower, df$Q_upper, df$model_lower, df$model_upper)
    title <- sprintf("%s\nKGE %.2f | NSE %.2f | r %.2f | nBIAS %+.2f\nAccuracy %.2f | Miss bias %+.2f | Miss mag %.2f\nBcorr %.2f | Precision %s",
                     df$scenario[1], t$KGE, t$NSE, t$r, t$nBIAS,
                     b$Accuracy, b$Miss_bias, b$Miss_magnitude, b$Bcorr_rho0,
                     if (is.na(b$Precision)) "n/a" else sprintf("%.2f", b$Precision))

    # Batesian panels: the model band is drawn on every panel so the legend key exists; where the model has
    # no band it is zero width (lower = upper = model) and invisible
    bates <- grp == "bates"
    p <- ggplot() + gauge
    if (bates) {
      p <- p + geom_ribbon(data = mutate(df, lo = ifelse(is.na(model_lower), model, model_lower),
                                            up = ifelse(is.na(model_upper), model, model_upper)),
                           aes(day, ymin = lo, ymax = up, fill = "Model"))
    }
    p + geom_point(data = df, aes(day, model, colour = "Model"), size = 1.8) +
      key_scales(if (bates) "Model and its band" else "Model", if (bates) alpha("grey30", 0.35) else "transparent") +
      labs(title = title) + th +
      if (show_legend) legend_inside else theme(legend.position = "none")
  }

  figure <- function(grp, title) {
    df <- filter(op, group == grp)
    panels <- split(df, factor(df$scenario, levels = unique(df$scenario)))
    patchwork::wrap_plots(lapply(seq_along(panels), function(i) panel(panels[[i]], i == 1, grp)), ncol = 2) +
      patchwork::plot_annotation(title = title, theme = theme(plot.title = element_text(size = 16, face = "bold")))
  }

  # Monte Carlo figure: draws on top, A vs B test table under each scenario -------------------
  mc_cols  <- c(Gauge = "black", A = "#2166ac", B = "#b2182b")
  mc_stats <- c("KGE", "NSE", "r", "nBIAS", "Miss_bias", "Miss_magnitude")
  fmt_p    <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

  mc_curve_panel <- function(nm) {
    ggplot(filter(example_data$mc_curves, scenario == nm), aes(day, Q, group = interaction(who, id), colour = who)) +
      geom_line(alpha = 0.15, linewidth = 0.5) +
      scale_colour_manual(values = mc_cols, breaks = c("Gauge", "A", "B"),
                          labels = c("Gauge draws", "Model A draws", "Model B draws"), name = NULL) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.2))) +
      labs(title = nm, x = "Day", y = expression(Q~(m^3/s))) +
      theme_bw(base_size = 13) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
            legend.text = element_text(size = 14), legend.key.width = unit(1.6, "lines"),
            legend.spacing.x = unit(0.4, "lines"))
  }

  # share of draws in which each model ranks first on each stat (both rank 1 = tie)
  mc_table <- function(nm) {
    res <- filter(example_data$mc_ranks, scenario == nm)
    tab <- res %>%
      group_by(metric, draw) %>%
      summarise(winner = if (sum(rank == 1) > 1) "Tie" else model[rank == 1], .groups = "drop") %>%
      count(metric, winner) %>%
      complete(metric, winner = c("A", "B", "Tie"), fill = list(n = 0)) %>%
      group_by(metric) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
      select(-n) %>% pivot_wider(names_from = winner, values_from = pct) %>%
      mutate(metric = factor(metric, levels = mc_stats)) %>% arrange(metric)
    out <- data.frame(Stat = as.character(tab$metric),
                      `A ranks 1st (%)` = sprintf("%.0f", tab$A), `B ranks 1st (%)` = sprintf("%.0f", tab$B),
                      `Tie (%)` = sprintf("%.0f", tab$Tie), check.names = FALSE)
    patchwork::wrap_elements(full = gridExtra::tableGrob(out, rows = NULL,
      theme = gridExtra::ttheme_minimal(base_size = 12, colhead = list(fg_params = list(fontface = "bold")))))
  }

  mc_names <- unique(example_data$mc_ranks$scenario)
  n_draw   <- max(example_data$mc_ranks$draw)
  n_curve  <- max(example_data$mc_curves$id)
  mc_fig <- patchwork::wrap_plots(mc_curve_panel(mc_names[1]), mc_curve_panel(mc_names[2]),
                                  mc_table(mc_names[1]), mc_table(mc_names[2]),
                                  ncol = 2, heights = c(1, 0.7)) +
    patchwork::plot_annotation(
      title = paste0("Monte Carlo ranking: ", n_draw, " draws, each model drawn against the same gauge draw"),
      theme = theme(plot.title = element_text(size = 16, face = "bold")))

  list(trad  = figure("trad", "Traditional metrics"),
       bates = figure("bates", "Batesian metrics"),
       mc    = mc_fig)
}