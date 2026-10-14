source("schematic_data.R")

# gauge log-space standard error at each overpass, from its own (flow-dependent) band
s  <- half_width(gauge) / 1.96
s1 <- head(s, -1)
s2 <- tail(s, -1)
ci_half <- function(rho) 1.96 * sqrt(s1^2 + s2^2 - 2 * rho * s1 * s2)

steps <- data.frame(mid_day = (head(day, -1) + tail(day, -1)) / 2,
                    d_gauge = diff(log(gauge)))

bcorr <- function(d_model, rho) mean(abs(d_model - steps$d_gauge) <= ci_half(rho))

y_lim_d <- range(c(steps$d_gauge - ci_half(0), steps$d_gauge + ci_half(0),
                   unlist(lapply(scenarios, function(m) diff(log(m))))))

hit_levels <- c("Model change, hit (rho = 0.8)", "Model change, miss (rho = 0.8)")

ci_layers <- list(
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.4),
  geom_rect(data = steps, aes(xmin = mid_day - 3, xmax = mid_day + 3,
                              ymin = d_gauge - ci_half(0), ymax = d_gauge + ci_half(0),
                              fill = "Gauge change 95% CI, rho = 0")),
  geom_rect(data = steps, aes(xmin = mid_day - 3, xmax = mid_day + 3,
                              ymin = d_gauge - ci_half(0.8), ymax = d_gauge + ci_half(0.8),
                              fill = "Gauge change 95% CI, rho = 0.8")),
  scale_fill_manual(values = c("Gauge change 95% CI, rho = 0"   = "grey75",
                               "Gauge change 95% CI, rho = 0.8" = "black")),
  scale_y_continuous(limits = y_lim_d),
  scale_x_continuous(limits = range(daily$day))
)

truth_steps <- ggplot() + ci_layers +
  labs(title = "Gauge change between\nconsecutive overpasses") +
  base_theme

columns <- lapply(names(scenarios), function(name) {
  d_model <- diff(log(scenarios[[name]]))
  hit     <- abs(d_model - steps$d_gauge) <= ci_half(0.8)

  top <- model_hydro_panel(name, sprintf("%s\nPearson r %.2f  |  B(0) %.2f  |  B(0.8) %.2f",
                                         name, cor(scenarios[[name]], gauge),
                                         bcorr(d_model, 0), bcorr(d_model, 0.8)))

  bottom <- ggplot() + ci_layers +
    geom_point(data = data.frame(steps, d_model,
                                 hit = factor(ifelse(hit, hit_levels[1], hit_levels[2]), levels = hit_levels)),
               aes(x = mid_day, y = d_model, shape = hit), colour = model_blue, size = 2.6, stroke = 1.1) +
    scale_shape_manual(values = setNames(c(16, 1), hit_levels), drop = FALSE) +
    guides(shape = guide_legend(override.aes = list(colour = model_blue))) +
    base_theme

  top / bottom
})

fig <- wrap_plots(c(list(truth_hydro_panel() / truth_steps), columns), design = "ABC\nADE") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical")


                                 plot(fig)
# ggsave("bcorr_schematic.png", fig, width = 20, height = 9, dpi = 150, bg = "white")
