library(glmmTMB)
library(seqwrap)
library(splines)
library(marginaleffects)
library(dplyr)
library(ggplot2)
library(broom.mixed)
library(purrr)
library(tidyr)

# Load metadata
metadata <- readRDS("data/Trainome_metadata.RDS")



# Load the beta_binomial model

beta_binom_model <- readRDS("data/splined_beta_binomial_model.RDS")

# Filter out the NULL models while preserving the names of valid ones
valid_models <- compact(beta_binom_model@models)

# test training effecs across all 100 timepoints


df_predictions <- map_dfr(
  valid_models,
  function(mod) {
    tryCatch({
      # Force standard error calculation using vcov = TRUE
      predictions(
        mod, 
        type = "response", 
        re_formula = NA,
        vcov = TRUE, 
        newdata = datagrid(scaled_age = seq(0, 1, length.out = 100),
                           time = c("PreExc", "PostExc")) # Forces predictions for BOTH states
      )
    }, error = function(e) {
      # If a model fails the matrix math, return NULL and continue the loop
      return(NULL) 
    })
  },
  .id = "target"
)

# saveRDS(df_predictions, "data/splined_UPDATED_beta_binomial_predictions.RDS")

# calculate slope of age across grid
# Rinstead of globally

# COMPARE WIH FIRST VERSION

df_age_slopes <- map_dfr(
  valid_models,
  function(mod) {
    tryCatch({
      # Calculate age slopes explicitly at both time points
      avg_slopes(
        mod, 
        variables = "scaled_age", 
        by = "time",
        newdata = datagrid(scaled_age = seq(0, 1, by = 0.25)), # Slopes at milestones
        type = "response", 
        re_formula = NA
      )
    }, error = function(e) return(NULL))
  },
  .id = "target"
)

# saveRDS(df_age_slopes, "data/splined_UPDATED_beta_binomial_age_slopes.RDS")



# The Marginal Effect of Time Across Age
df_time_effects <- map_dfr(
  valid_models,
  function(mod) {
    tryCatch({
      # Calculate the discrete change of time (post vs pre) across age steps
      avg_comparisons(
        mod,
        variables = "time",
        by = "scaled_age",
        # Test at specific age anchors: 0, 0.25, 0.5, 0.75, and 1.0
        newdata = datagrid(scaled_age = seq(0, 1, by = 0.25)),
        type = "response",
        re_formula = NA
      )
    }, error = function(e) return(NULL))
  },
  .id = "target"
)


# df_time_effects <- readRDS("data/splined_beta_binomial_time_effect.RDS")
# saveRDS(df_time_effects, "data/splined_beta_binomial_time_effect.RDS")

# REPEAT FOR BINOMIAL MODEL

binomial_model <- readRDS("data/splined_binom_results.RDS")


valid_models1 <- compact(binomial_model@models) 



df_predictions1 <- map_dfr(
  valid_models1,
  function(mod) {
    tryCatch({
      # Force standard error calculation using vcov = TRUE
      predictions(
        mod, 
        type = "response", 
        re_formula = NA,
        vcov = TRUE, 
        newdata = datagrid(scaled_age = seq(0, 1, length.out = 100),
                           time = c("PreExc", "PostExc")) # Forces predictions for BOTH states
      )
    }, error = function(e) {
      # If a model fails the matrix math, return NULL and continue the loop
      return(NULL) 
    })
  },
  .id = "target"
)

# saveRDS(df_predictions1, "data/splined_UPDATED_binomial_predictions.RDS")





df_age_slopes1 <- map_dfr(
  valid_models1,
  function(mod) {
    tryCatch({
      # Calculate age slopes explicitly at both time points
      avg_slopes(
        mod, 
        variables = "scaled_age", 
        by = "time",
        newdata = datagrid(scaled_age = seq(0, 1, by = 0.25)), # Slopes at milestones
        type = "response", 
        re_formula = NA
      )
    }, error = function(e) return(NULL))
  },
  .id = "target"
)

# saveRDS(df_age_slopes1, "data/splined_UPDATED_binomial_age_slopes.RDS")



# VISUALISATION

# extract introns that are significatly active by time

# 1. Find introns where the exercise training effect is significantly active at ANY age milestone
fdr_significant_targets <- df_time_effects %>%
  group_by(target) %>%
  # Adjust p-values across all models and milestones globally using Benjamini-Hochberg
  mutate(adj.p = p.adjust(p.value, method = "BH")) %>% 
   filter(adj.p < 0.05) # %>%
  # pull(target) %>%
  # unique()


# Visualise the non-linear relatiohsip between aging and intron splicing efficiency


# Filter predictions down to your robust FDR targets
df_sig_predictions <- df_predictions %>%
  filter(target %in% fdr_significant_targets$target)

ggplot(df_sig_predictions, aes(x = scaled_age, y = estimate, color = time, fill = time)) +
  # 1. Plot all individual intron trajectories as ultra-faint lines to visualize the density/spread
  geom_line(aes(group = interaction(target, time)), alpha = 0.01, linewidth = 0.8) +
  # 2. Overlay a powerful non-linear generalized additive model (GAM) smoother for the consensus trend
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, alpha = 0.2) +
  scale_color_manual(values = c("PreExc" = "#e41a1c", "PostExc" = "#377eb8")) +
  scale_fill_manual(values = c("PreExc" = "#e41a1c", "PostExc" = "#377eb8")) +
  coord_cartesian(ylim = c(0.9,1)) +
  labs(
    title = "Global Non-Linear Splicing Dynamics Across Aging",
    subtitle = paste("Consensus trends for", length(unique(fdr_significant_targets$target)), "FDR-significant introns"),
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Splicing Efficiency (Response Scale)",
    color = "Exercise State",
    fill = "Exercise State"
  ) +
  theme_minimal()




# quantifying how aging changes SE both at baseline and post exercise


# 1. Calculate the change in splicing efficiency between age steps
df_derived_slopes <- df_predictions %>%
  filter(target %in% fdr_significant_targets$target) %>%
  group_by(target, time) %>%
  arrange(scaled_age) %>%
  # Change in efficiency divided by change in age = Slope
  mutate(slope = (estimate - lag(estimate)) / (scaled_age - lag(scaled_age))) %>%
  filter(!is.na(slope))

# 2. Plot the distribution of these trajectories over the continuous age scale
ggplot(df_derived_slopes, aes(x = scaled_age, y = slope, color = time, fill = time)) +
  # Faint background lines for individual introns
  geom_line(aes(group = interaction(target, time)), alpha = 0.01, linewidth = 0.5) +
  # Bold global consensus non-linear spline trend
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("PreExc" = "#e41a1c", "PostExc" = "#377eb8")) +
  scale_fill_manual(values = c("PreExc" = "#e41a1c", "PostExc" = "#377eb8")) +
  labs(
    title = "Splicing Decay Velocity Across the Aging Lifecycle",
    subtitle = "Values below 0 indicate active splicing degradation at that specific age",
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Splicing Efficiency Slope (Rate of Change)",
    color = "Exercise",
    fill = "Exercise"
  ) +
  theme_minimal()




# Define the boundaries of age sequence from the data
min_age <- min(df_predictions$scaled_age)
max_age <- max(df_predictions$scaled_age)

# Step 1: Calculate the absolute baseline aging shift per intron
fdr_age_active_targets <- df_predictions %>%
  # Focus strictly on the natural aging state
  filter(time == "PreExc") %>%
  # Grab the absolute youngest and oldest prediction points
  filter(scaled_age == min_age | scaled_age == max_age) %>%
  mutate(age_status = if_else(scaled_age == min_age, "Young", "Old")) %>%
  select(target, age_status, estimate) %>%
  # Pivot side-by-side to find the lifespan delta
  pivot_wider(names_from = age_status, values_from = estimate) %>%
  mutate(lifespan_change = Old - Young) %>%
  # Keep introns in the top 10% of changes (severe decline or extreme increase)
  filter(lifespan_change < quantile(lifespan_change, 0.10, na.rm = TRUE) | 
           lifespan_change > quantile(lifespan_change, 0.90, na.rm = TRUE)) %>%
  pull(target)

# Step 2: Subset predictions to only your age-active intron list
df_plot_data <- df_predictions %>%
  filter(target %in% fdr_age_active_targets)

# Step 3: Plot the interaction between non-linear aging and exercise
ggplot(df_plot_data, aes(x = scaled_age, y = estimate, color = time, fill = time)) +
  # 1. Background: Draw all trajectories ultra-faintly to show global density
  geom_line(aes(group = interaction(target, time)), alpha = 0.01, linewidth = 0.4) +
  # 2. Foreground: Bold GAM smoother showing the absolute biological consensus trend
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, alpha = 0.2) +
  # Use highly readable, publication-standard colors
  scale_color_manual(values = c("PreExc" = "#d95f02", "PostExc" = "#1b9e77")) +
  scale_fill_manual(values = c("PreExc" = "#d95f02", "PostExc" = "#1b9e77")) +
  coord_cartesian(ylim = c(0.95, 1)) +
  labs(
    title = "Global Non-Linear Splicing Trajectories for Age-Active Introns",
    subtitle = paste("Showing consensus trends for", length(fdr_age_active_targets), "highly responsive aging targets"),
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Predicted Splicing Efficiency (Response Scale)",
    color = "Training State",
    fill = "Training State"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


# Visualise genes significant in old age
# Find targets significant ONLY in advanced age
old_age_responders <- df_time_effects %>%
  # 1. Correct p-values globally within each specific age group
  group_by(scaled_age) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  # 2. Filter for robust significance (FDR < 0.05) strictly at advanced milestones
  filter(fdr < 0.05 & scaled_age >= 0.75) %>%
  pull(target) %>%
  unique()

# Find targets significant in youth to subtract them (ensuring age-specificity)
young_age_responders <- df_time_effects %>%
  group_by(scaled_age) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  filter(fdr < 0.05 & scaled_age <= 0.25) %>%
  pull(target) %>%
  unique()

# Final Filter: Keep genes significant in old age, but NOT in youth
strict_old_age_targets <- setdiff(old_age_responders, young_age_responders)





df_specific_plot <- df_predictions %>%
  filter(target %in% strict_old_age_targets)

ggplot(df_specific_plot, aes(x = scaled_age, y = estimate, color = time, fill = time)) +
  # 1. Individual trajectories of your uniquely old-age responsive genes
  geom_line(aes(group = interaction(target, time)), alpha = 0.02, linewidth = 0.4) +
  # 2. Bold universal trend curve
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, alpha = 0.2) +
  scale_color_manual(values = c("PreExc" = "#7570b3", "PostExc" = "#e7298a")) +
  scale_fill_manual(values = c("PreExc" = "#7570b3", "PostExc" = "#e7298a")) +
  labs(
    title = "Exercise Rescue Effect is Limited Exclusively to Older Ages",
    subtitle = paste("Consensus trajectories for", length(strict_old_age_targets), "uniquely old-age sensitive introns"),
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Predicted Splicing Efficiency",
    color = "Training State",
    fill = "Training State"
  ) +
  theme_minimal()





min_age <- min(df_predictions$scaled_age)
max_age <- max(df_predictions$scaled_age)

# 1. Identify targets with a large lifespan change under baseline conditions
natural_aging_metrics <- df_predictions %>%
  filter(time == "PreExc") %>%
  filter(scaled_age == min_age | scaled_age == max_age) %>%
  mutate(age_status = if_else(scaled_age == min_age, "Young", "Old")) %>%
  select(target, age_status, estimate) %>%
  pivot_wider(names_from = age_status, values_from = estimate) %>%
  mutate(lifespan_change = Old - Young)

# 2. Select the top 10% most severe natural aging targets (both up and down)
severe_aging_targets <- natural_aging_metrics %>%
  filter(lifespan_change < quantile(lifespan_change, 0.10, na.rm = TRUE) | 
           lifespan_change > quantile(lifespan_change, 0.90, na.rm = TRUE)) %>%
  pull(target)

# 3. Calculate the global exercise effect (Post - Pre) across all age steps
exercise_effects_summary <- df_predictions %>%
  filter(target %in% severe_aging_targets) %>%
  select(target, scaled_age, time, estimate) %>%
  pivot_wider(names_from = time, values_from = estimate) %>%
  mutate(exercise_diff = abs(PostExc - PreExc)) %>%
  group_by(target) %>%
  # Find the maximum exercise response across the whole lifespan
  summarise(max_exercise_response = max(exercise_diff, na.rm = TRUE))

# 4. Strict Filter: Keep severe aging targets with very little exercise response 
# (e.g., exercise effect is in the bottom 30% of your aging subset)
age_only_targets <- exercise_effects_summary %>%
  filter(max_exercise_response < quantile(max_exercise_response, 0.30, na.rm = TRUE)) %>%
  pull(target)



# 1. Directional mapping for faceting
direction_map <- natural_aging_metrics %>%
  filter(target %in% age_only_targets) %>%
  mutate(aging_direction = if_else(lifespan_change < 0, "Splicing Decreases with Age", "Splicing Increases with Age")) %>%
  select(target, aging_direction)

# 2. Prepare the plotting data
df_age_only_plot <- df_predictions %>%
  filter(target %in% age_only_targets) %>%
  inner_join(direction_map, by = "target")

# 3. Plot the parallel aging tracks
ggplot(df_age_only_plot, aes(x = scaled_age, y = estimate, color = time, fill = time)) +
  # Background individual trajectories
  geom_line(aes(group = interaction(target, time)), alpha = 0.01, linewidth = 0.4) +
  # Bold consensus GAM curves
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, alpha = 0.2) +
  # Split the plot side-by-side by direction of aging
  facet_wrap(~aging_direction, scales = "free_y") +
  scale_color_manual(values = c("PreExc" = "#33a02c", "PostExc" = "#b2df8a")) +
  scale_fill_manual(values = c("PreExc" = "#33a02c", "PostExc" = "#b2df8a")) +
  labs(
    title = "Irreversible Splicing Trajectories Driven Strictly by Aging",
    subtitle = paste("Showing consensus trajectories for", length(age_only_targets), "exercise-insensitive aging markers"),
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Predicted Splicing Efficiency",
    color = "Training State",
    fill = "Training State"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12), # Formats facet headers
    legend.position = "bottom"
  )


library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Isolate the subset of introns that are statistically robust aging targets
fdr_aging_introns <- df_age_slopes %>%
  # Apply Benjamini-Hochberg FDR correction globally across all 5,000 models
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  # Keep only the statistically definitive aging relationships
  filter(fdr < 0.05) %>%
  pull(target) %>%
  unique()

message("Statistically confirmed ", length(fdr_aging_introns), " FDR-significant aging introns.")

# 2. Of these significant introns, find which ones have parallel curves (insensitive to exercise)
min_age <- min(df_predictions$scaled_age)
max_age <- max(df_predictions$scaled_age)

exercise_insensitivity_matrix <- df_predictions %>%
  # Keep only the statistically verified aging introns
  filter(target %in% fdr_aging_introns) %>%
  select(target, scaled_age, time, estimate) %>%
  pivot_wider(names_from = time, values_from = estimate) %>%
  # Calculate the absolute difference between Pre and Post curves across the whole lifecycle
  mutate(exercise_diff = abs(PostExc - PreExc)) %>%
  group_by(target) %>%
  summarise(max_exercise_response = max(exercise_diff, na.rm = TRUE))

# 3. Final Filter: Keep FDR aging introns where the exercise response is minimal
# (e.g., the bottom 30% of the significant aging pool)
strict_fdr_age_only_targets <- exercise_insensitivity_matrix %>%
  filter(max_exercise_response < quantile(max_exercise_response, 0.30, na.rm = TRUE)) %>%
  pull(target)

message("Isolated ", length(strict_fdr_age_only_targets), " strictly FDR-significant, irreversible aging markers.")



aging_direction_map <- df_predictions %>%
  filter(target %in% strict_fdr_age_only_targets & time == "PreExc") %>%
  filter(scaled_age == min_age | scaled_age == max_age) %>%
  mutate(age_status = if_else(scaled_age == min_age, "Young", "Old")) %>%
  select(target, age_status, estimate) %>%
  pivot_wider(names_from = age_status, values_from = estimate) %>%
  mutate(aging_direction = if_else(Old - Young < 0, "Splicing Decreases with Age (FDR Sig)", "Splicing Increases with Age (FDR Sig)")) %>%
  select(target, aging_direction)


df_final_plot_data <- df_predictions %>%
  filter(target %in% strict_fdr_age_only_targets) %>%
  inner_join(aging_direction_map, by = "target")

ggplot(df_final_plot_data, aes(x = scaled_age, y = estimate, color = time, fill = time)) +
  # 1. Background density: Faint lines for true FDR-significant aging trajectories
  geom_line(aes(group = interaction(target, time)), alpha = 0.02, linewidth = 0.4) +
  # 2. Foreground consensus: Bold non-linear GAM curves showing overlapping tracks
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, alpha = 0.2) +
  facet_wrap(~aging_direction, scales = "free_y") +
  # Use complementary colors to show overlapping, parallel tracks
  scale_color_manual(values = c("PreExc" = "#4daf4a", "PostExc" = "#984ea3")) +
  scale_fill_manual(values = c("PreExc" = "#4daf4a", "PostExc" = "#984ea3")) +
  labs(
    title = "Irreversible Splicing Trajectories Driven Strictly by Aging",
    subtitle = paste("Consensus curves for", length(strict_fdr_age_only_targets), "FDR-significant, exercise-insensitive introns"),
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Predicted Splicing Efficiency (Response Scale)",
    color = "Training State",
    fill = "Training State"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )



library(dplyr)
library(ggplot2)

# 1. Calculate the moving velocity (slope) across continuous age steps for your FDR aging targets
df_splicing_velocity <- df_predictions %>%
  # Filter for your FDR-significant aging targets (from your previous filtering step)
  filter(target %in% strict_fdr_age_only_targets & time == "PreExc") %>%
  group_by(target) %>%
  arrange(scaled_age) %>%
  # Slope = change in splicing efficiency divided by change in age
  mutate(
    velocity = (estimate - lag(estimate)) / (scaled_age - lag(scaled_age))
  ) %>%
  filter(!is.na(velocity))

# 2. Plot the velocity curve to find the exact global "dip" point
ggplot(df_splicing_velocity, aes(x = scaled_age, y = velocity)) +
  # Faint background lines for all introns
  geom_line(aes(group = target), alpha = 0.01, color = "gray") +
  # Bold global consensus curve using a GAM smoother
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red", linewidth = 1.5) +
  # Reference line at 0 (above 0 = increasing, below 0 = decreasing)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Splicing Decay Velocity Across the Aging Lifespan",
    subtitle = "The lowest point on the red curve marks the exact age where the 'dip' is steepest",
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Splicing Acceleration/Decay Velocity"
  ) +
  theme_minimal()


# 1. Find the exact age of the absolute steepest drop for each individual intron
df_steepest_dips <- df_splicing_velocity %>%
  group_by(target) %>%
  # Filter for the single age step where the negative velocity was at its absolute maximum
  filter(velocity == min(velocity)) %>%
  select(target, age_of_steepest_dip = scaled_age, worst_velocity = velocity)

# 2. Plot a density curve to see where these dips cluster together
ggplot(df_steepest_dips, aes(x = age_of_steepest_dip)) +
  geom_density(fill = "darkblue", alpha = 0.4, color = "darkblue", linewidth = 1) +
  labs(
    title = "Distribution of Splicing Inflection Points Across Age",
    subtitle = "Peaks indicate specific age milestones where a massive wave of introns experience a dip",
    x = "Scaled Age at Which Splicing Dips",
    y = "Density of Dipping Introns"
  ) +
  theme_minimal()


library(dplyr)

# 1. Calculate the average global splicing trajectory across your selected introns
global_consensus_profile <- df_splicing_velocity %>%
  group_by(scaled_age) %>%
  summarise(
    avg_splicing_efficiency = mean(estimate, na.rm = TRUE),
    avg_velocity = mean(velocity, na.rm = TRUE)
  )

# 2. Isolate the exact bottom of the mid-life dip (searching around the 0.3 - 0.6 window)
dip_data <- global_consensus_profile %>%
  filter(scaled_age >= 0.30 & scaled_age <= 0.60) %>%
  filter(avg_splicing_efficiency == min(avg_splicing_efficiency))

exact_dip_age <- dip_data$scaled_age[1]

# 3. Isolate the exact peak of the late-life rise (searching around the 0.6 - 0.9 window)
rise_data <- global_consensus_profile %>%
  filter(scaled_age >= 0.60 & scaled_age <= 0.90) %>%
  filter(avg_splicing_efficiency == max(avg_splicing_efficiency))

exact_rise_age <- rise_data$scaled_age[1]

# 4. Print the exact metrics out
cat(paste0(
  "--- CRITICAL INFLECTION POINTS ---\n",
  "Exact Age of the Bottom of the Dip: ", round(exact_dip_age, 4), "\n",
  "Exact Age of the Peak of the Rise:  ", round(exact_rise_age, 4), "\n"
))


raw_min_age <- min(metadata$age, na.rm = TRUE) 
raw_max_age <- max(metadata$age, na.rm = TRUE)

biological_dip_age  <- (exact_dip_age * (raw_max_age - raw_min_age)) + raw_min_age
biological_rise_age <- (exact_rise_age * (raw_max_age - raw_min_age)) + raw_min_age

cat(paste0(
  "Biological Age at bottom of Dip: ", round(biological_dip_age, 1), " years old\n",
  "Biological Age at peak of Rise:   ", round(biological_rise_age, 1), " years old\n"
))

# Use your previous plotting script as the base, and add these layers:
ggplot(df_plot_data, aes(x = scaled_age, y = estimate, color = time, fill = time)) +
  geom_line(aes(group = interaction(target, time)), alpha = 0.01, linewidth = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, alpha = 0.2) +
  
  # Add vertical line for the exact bottom of the dip
  geom_vline(xintercept = exact_dip_age, linetype = "dotted", color = "darkred", linewidth = 1) +
  # Add vertical line for the exact top of the rebound
  geom_vline(xintercept = exact_rise_age, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  
  # Text annotations on the chart
  annotate("text", x = exact_dip_age, y = max(df_plot_data$estimate)*0.9, 
           label = paste("Dip:", round(biological_dip_age, 1), "yrs"), color = "darkred", hjust = 1.1, fontface = "bold") +
  annotate("text", x = exact_rise_age, y = max(df_plot_data$estimate)*0.9, 
           label = paste("Rise:", round(biological_rise_age, 1), "yrs"), color = "darkgreen", hjust = -0.1, fontface = "bold") +
  
  scale_color_manual(values = c("PreExc" = "#d95f02", "PostExc" = "#1b9e77")) +
  scale_fill_manual(values = c("PreExc" = "#d95f02", "PostExc" = "#1b9e77")) +
  theme_minimal()



library(dplyr)
library(tidyr)

# 1. Isolate the upward velocity window and pivot it side-by-side
df_paired_velocity_wide <- df_predictions %>%
  filter(target %in% strict_fdr_age_only_targets) %>%
  group_by(target, time) %>%
  arrange(scaled_age) %>%
  # Velocity = change in splicing efficiency per step of age
  mutate(velocity = (estimate - lag(estimate)) / (scaled_age - lag(scaled_age))) %>%
  filter(scaled_age > 0.50 & scaled_age <= 0.70) %>%
  summarise(mean_rise_velocity = mean(velocity, na.rm = TRUE), .groups = "drop") %>%
  # FIX: Pivot the time variable into separate columns so they pair perfectly by target
  pivot_wider(names_from = time, values_from = mean_rise_velocity)

# 2. Run the paired t-test using the vector inputs
paired_velocity_test <- t.test(
  df_paired_velocity_wide$PostExc, 
  df_paired_velocity_wide$PreExc, 
  paired = TRUE
)

print(paired_velocity_test)


# 1. Isolate the exact peak point (scaled_age 0.70) and pivot wide
df_paired_height_wide <- df_predictions %>%
  filter(target %in% strict_fdr_age_only_targets) %>%
  filter(abs(scaled_age - 0.70) == min(abs(scaled_age - 0.70))) %>%
  select(target, time, estimate) %>%
  # FIX: Pivot the time variable into separate columns to pair by intron target
  pivot_wider(names_from = time, values_from = estimate)

# 2. Run the paired t-test on absolute efficiency at the peak
paired_height_test <- t.test(
  df_paired_height_wide$PostExc, 
  df_paired_height_wide$PreExc, 
  paired = TRUE
)

print(paired_height_test)


# 1. Calculate the moving velocity (slope) across continuous age steps for your FDR aging targets
df_splicing_velocity <- df_predictions %>%
  # Filter for your FDR-significant aging targets (from your previous filtering step)
  filter(target %in% strict_fdr_age_only_targets & time == "PreExc") %>%
  group_by(target) %>%
  arrange(scaled_age) %>%
  # Slope = change in splicing efficiency divided by change in age
  mutate(
    velocity = (estimate - lag(estimate)) / (scaled_age - lag(scaled_age))
  ) %>%
  filter(!is.na(velocity))

# 2. Plot the velocity curve to find the exact global "dip" point
ggplot(df_splicing_velocity, aes(x = scaled_age, y = velocity)) +
  # Faint background lines for all introns
  geom_line(aes(group = target), alpha = 0.01, color = "gray") +
  # Bold global consensus curve using a GAM smoother
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red", linewidth = 1.5) +
  # Reference line at 0 (above 0 = increasing, below 0 = decreasing)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Splicing Decay Velocity Across the Aging Lifespan",
    subtitle = "The lowest point on the red curve marks the exact age where the 'dip' is steepest",
    x = "Scaled Age (0 = Youngest, 1 = Oldest)",
    y = "Splicing Acceleration/Decay Velocity"
  ) +
  theme_minimal()

