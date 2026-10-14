# Shared synthetic data for both schematic figures: same hydrograph, band, overpasses and models.
library(ggplot2)
library(patchwork)

set.seed(11)

hydro <- function(day) {
  50 + 30 * sin(2 * pi * (day - 60) / 365) +
    120 * exp(-((day - 120) / 12)^2) + 80 * exp(-((day - 250) / 15)^2)
}

# 95% band half-width in log space grows with flow: 0.25 at the lowest flow, 0.65 at the highest
log_range <- range(log(hydro(1:365)))
half_width <- function(gauge) 0.25 + 0.40 * (log(gauge) - log_range[1]) / diff(log_range)

daily <- data.frame(day = 1:365)
daily$gauge   <- hydro(daily$day)
daily$LowerCI <- daily$gauge * exp(-half_width(daily$gauge))
daily$UpperCI <- daily$gauge * exp( half_width(daily$gauge))

# overpasses
n       <- 26
day     <- sort(sample(seq(5, 360, by = 7), n))
gauge   <- hydro(day)
LowerCI <- gauge * exp(-half_width(gauge))
UpperCI <- gauge * exp( half_width(gauge))

log_centre <- mean(log(gauge))
scenarios <- list(
  "Within the band"   = gauge * exp(rnorm(n, 0,    0.08)),
  "Biased high"       = gauge * 2.2 * exp(rnorm(n, 0, 0.05)),
  "Unbiased, noisy"   = gauge * exp(rnorm(n, 0,    0.45)),
  "Too smooth"        = exp(log_centre + 0.2 * (log(gauge) - log_centre))
)

y_lim_h <- range(c(daily$LowerCI, daily$UpperCI, unlist(scenarios)))
model_blue <- "#3a7dc9"

base_theme <- theme_void(base_size = 14) +
  theme(panel.border    = element_rect(colour = "grey70", fill = NA),
        plot.title      = element_text(face = "bold", hjust = 0.5, size = 12.5, lineheight = 1.1),
        plot.margin     = margin(4, 6, 4, 6),
        legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 12))

# identical hydrograph panels for both figures
truth_hydro_panel <- function() {
  ggplot(daily, aes(day)) +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = "Gauge 95% band (outline in model panels)")) +
    scale_fill_manual(values = c("Gauge 95% band (outline in model panels)" = "black")) +
    scale_y_log10(limits = y_lim_h) +
    labs(title = "Truth\n(gauge 95% band)") +
    base_theme
}

model_hydro_panel <- function(name, title) {
  ggplot(daily, aes(day)) +
    geom_line(aes(y = LowerCI), linewidth = 0.5) +
    geom_line(aes(y = UpperCI), linewidth = 0.5) +
    geom_point(data = data.frame(day, model = scenarios[[name]]),
               aes(y = model, colour = "Model prediction at overpass"), size = 2.4) +
    scale_colour_manual(values = c("Model prediction at overpass" = model_blue)) +
    scale_y_log10(limits = y_lim_h) +
    labs(title = title) +
    base_theme
}
