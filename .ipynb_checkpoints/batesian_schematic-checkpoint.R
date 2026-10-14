source("schematic_data.R")

batesian_stats <- function(model) {
  is_OP <- model > UpperCI
  is_UP <- model < LowerCI
  miss  <- is_OP | is_UP
  # Miss magnitude: log distance outside the band edge crossed; symmetric for over/under
  c(Accuracy       = mean(!miss),
    Miss_bias      = (sum(is_OP) - sum(is_UP)) / n,
    Miss_magnitude = if (!any(miss)) 0 else
                       mean(c(log(model[is_OP] / UpperCI[is_OP]), log(LowerCI[is_UP] / model[is_UP]))))
}

model_panels <- lapply(names(scenarios), function(name) {
  s <- batesian_stats(scenarios[[name]])
  model_hydro_panel(name, sprintf("%s\nAccuracy %.2f  |  Miss bias %+.2f\nMiss magnitude %.2f",
                                  name, s["Accuracy"], s["Miss_bias"], s["Miss_magnitude"]))
})
fig <- wrap_plots(c(list(truth_hydro_panel()), model_panels), design = "ABC\nADE") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical")


plot(fig)

# ggsave("batesian_schematic.png", fig, width = 20, height = 5.2, dpi = 150, bg = "white")
