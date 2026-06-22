# ============================================================
# Extract model outputs
# ============================================================

if (!compareMultipleOutputs) {

  outputs <- extract_model_outputs(path_to_results)

} else {

  paths_df <- build_model_paths(
    output_root = output_root,
    version_label = version_label,
    ecosystemType = ecosystemType,
    existing_only = TRUE,
    include_single_species = c(TRUE, FALSE),
    include_spatial_model = c(TRUE, FALSE)
  )

  outputs <- extract_model_outputs(paths_df$path_to_results)
}

# ============================================================
# Save extracted outputs
# ============================================================

save(
  outputs,
  file = file.path(path_to_postProcessing, "model_outputs.RData")
)

print(outputs)

# ============================================================
# Common plotting theme
# ============================================================

theme_model <- function() {

  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.2)
    )
}

# ============================================================
# Create plots
# ============================================================

p_occ_det <- plot(
  outputs,
  type = "error",
  components = c("det_effects", "occ_effects")
) +
  ggplot2::labs(
    title = "Detection and occupancy effects",
    subtitle = ecosystemType
  ) +
  theme_model()

p_bayes <- plot(
  outputs,
  type = "error",
  components = c("bayes_P")
) +
  ggplot2::labs(
    title = "Posterior predictive checks",
    subtitle = ecosystemType,
    x = "Bayesian p-value"
  ) +
  theme_model()

p_comm <- plot(
  outputs,
  type = "error",
  components = c("beta_summary", "alpha_summary")
) +
  ggplot2::labs(
    title = "Community-level parameter summaries",
    subtitle = ecosystemType,
    x = "Posterior mean (95% credible interval)"
  ) +
  theme_model()

# ============================================================
# Save helper
# ============================================================

save_plot <- function(plot_obj,
                      filename,
                      width = 14,
                      height = 10) {

  ggplot2::ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(filename, ".png")
    ),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )

  ggplot2::ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(filename, ".pdf")
    ),
    plot = plot_obj,
    width = width,
    height = height
  )
}

# ============================================================
# Save figures
# ============================================================

save_plot(
  p_occ_det,
  "occupancy_detection_effects"
)

save_plot(
  p_bayes,
  "posterior_predictive_checks"
)

save_plot(
  p_comm,
  "community_parameter_summaries"
)

# ============================================================
# Display plots in current session
# ============================================================

print(p_occ_det)
print(p_bayes)
print(p_comm)
