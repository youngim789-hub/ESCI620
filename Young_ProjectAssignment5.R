# ============================================================
# FACILITY POLLUTION PROFILE CLUSTERING
# Approach: K-means
# Unit of analysis: One facility, identified by DEC ID
# Features: Average reported annual emissions for 9 pollutants
# Validation: Silhouette scores and elbow method
# ============================================================

# Run once if needed:
# install.packages(c("tidyverse", "cluster"))

library(tidyverse)
library(cluster)

set.seed(123)

# ------------------------------------------------------------
# 1. SELECT A MODELING APPROACH
# ------------------------------------------------------------

# Which facilities have similar average annual pollution profiles?
#
# K-means uses only emissions variables.
# Facility names, IDs, locations, years, and quality-control flags
# are not clustering features.

pollutants <- c(
  "VOC (tons)",
  "NOx (tons)",
  "CO (tons)",
  "CO2 (tons)",
  "Particulates (tons)",
  "PM10 (tons)",
  "PM2.5 (tons)",
  "HAPS (tons)",
  "SO2 (tons)"
)

raw <- read_csv(
  "EMISSIONS_CLEANED_111126.csv",
  col_types = cols(.default = col_character()),
  na = c("", "NA", "N/A", "NULL"),
  show_col_types = FALSE
)

required_columns <- c("DEC ID", "Facility Name", "Year", pollutants)

missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. PREPROCESS FOR MODELING
# ------------------------------------------------------------

# Preserve IDs as character strings and convert emissions to numeric.
# Negative or non-finite emissions are treated as missing.

clean <- raw %>%
  mutate(
    `DEC ID` = trimws(`DEC ID`),
    Year = suppressWarnings(as.integer(Year)),
    across(
      all_of(pollutants),
      ~ {
        x <- suppressWarnings(
          as.numeric(gsub(",", "", .x, fixed = TRUE))
        )
        x[!is.finite(x) | x < 0] <- NA_real_
        x
      }
    )
  ) %>%
  filter(
    !is.na(`DEC ID`),
    `DEC ID` != "",
    !is.na(Year)
  ) %>%
  select(all_of(required_columns)) %>%
  distinct()

# Conflicting records for the same facility and year should be reviewed rather than silently summed or averaged.
duplicate_facility_years <- clean %>%
  count(`DEC ID`, Year, name = "records") %>%
  filter(records > 1)

if (nrow(duplicate_facility_years) > 0) {
  write_csv(
    duplicate_facility_years,
    "facility_year_duplicates_to_review.csv"
  )
  
  stop(
    "Multiple non-identical records remain for some facility-years. ",
    "Review facility_year_duplicates_to_review.csv before modeling."
  )
}

# Avoid NaN when every value in a group is missing.
mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

first_nonmissing <- function(x) {
  valid <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(valid) == 0) NA_character_ else valid[1]
}

# One observation per facility.
# Means use available reported values, not an assumption that missing emissions equal zero.
facilities <- clean %>%
  group_by(`DEC ID`) %>%
  summarise(
    `Facility Name` = first_nonmissing(`Facility Name`),
    first_year = min(Year),
    last_year = max(Year),
    years_observed = n_distinct(Year),
    across(all_of(pollutants), mean_or_na),
    .groups = "drop"
  ) %>%
  mutate(
    missing_pollutants = rowSums(
      is.na(as.matrix(pick(all_of(pollutants))))
    )
  )

# Retain excluded facilities as an audit file.
excluded_facilities <- facilities %>%
  filter(missing_pollutants == length(pollutants))

write_csv(excluded_facilities, "facilities_excluded_all_missing.csv")

facilities <- facilities %>%
  filter(missing_pollutants < length(pollutants))

if (nrow(facilities) < 10) {
  stop("Too few usable facilities for this train/test workflow.")
}

# Split at the FACILITY level: approximately 80% training, 20% test.
train_index <- sample(
  seq_len(nrow(facilities)),
  size = floor(0.80 * nrow(facilities))
)

train <- facilities[train_index, ]
test <- facilities[-train_index, ]

x_train <- as.matrix(select(train, all_of(pollutants)))
x_test <- as.matrix(select(test, all_of(pollutants)))

# Learn median imputation values from training facilities only.
train_medians <- apply(x_train, 2, median, na.rm = TRUE)

if (any(!is.finite(train_medians))) {
  stop("At least one pollutant is entirely missing in training data.")
}

impute_matrix <- function(x, medians) {
  for (j in seq_len(ncol(x))) {
    x[is.na(x[, j]), j] <- medians[j]
  }
  x
}

x_train_imputed <- impute_matrix(x_train, train_medians)
x_test_imputed <- impute_matrix(x_test, train_medians)

# Log-transform to reduce the influence of very large emissions.
# log1p(x) = log(1 + x), so zero emissions remain valid.
x_train_log <- log1p(x_train_imputed)
x_test_log <- log1p(x_test_imputed)

# Remove features with no variation in the training data.
training_sd <- apply(x_train_log, 2, sd)
keep <- is.finite(training_sd) & training_sd > 0

if (!all(keep)) {
  message(
    "Dropping constant training features: ",
    paste(names(training_sd)[!keep], collapse = ", ")
  )
}

if (!any(keep)) {
  stop("No varying pollution features remain.")
}

used_pollutants <- pollutants[keep]

x_train_log <- x_train_log[, keep, drop = FALSE]
x_test_log <- x_test_log[, keep, drop = FALSE]

# Estimate scaling parameters on training data only.
training_mean <- colMeans(x_train_log)
training_sd <- apply(x_train_log, 2, sd)

standardize <- function(x, means, sds) {
  sweep(sweep(x, 2, means, "-"), 2, sds, "/")
}

x_train_scaled <- standardize(
  x_train_log, training_mean, training_sd
)

x_test_scaled <- standardize(
  x_test_log, training_mean, training_sd
)

# No one-hot encoding is needed:
# all clustering features are numeric.
# Facility Name and DEC ID are retained only as identifiers.

# ------------------------------------------------------------
# 3. TRAIN AND VALIDATE THE MODEL
# ------------------------------------------------------------

# Evaluate candidate k values using training data only.
# Cap k to avoid requesting more clusters than distinct profiles.
max_k <- min(
  10L,
  nrow(x_train_scaled) - 1L,
  nrow(unique(as.data.frame(x_train_scaled)))
)

if (max_k < 2) {
  stop("Not enough distinct training profiles for clustering.")
}

k_values <- 2:max_k

# Silhouette calculations require pairwise distances.
# Use a fixed sample to limit memory consumption.
silhouette_index <- sample(
  seq_len(nrow(x_train_scaled)),
  min(2000L, nrow(x_train_scaled))
)

silhouette_distances <- dist(
  x_train_scaled[silhouette_index, , drop = FALSE]
)

safe_silhouette <- function(labels, distances) {
  n_groups <- length(unique(labels))
  
  if (n_groups < 2 || n_groups >= length(labels)) {
    return(NA_real_)
  }
  
  mean(
    cluster::silhouette(as.integer(labels), distances)[, "sil_width"]
  )
}

models <- vector("list", length(k_values))

validation <- tibble(
  k = k_values,
  total_within_ss = NA_real_,
  mean_silhouette = NA_real_
)

for (i in seq_along(k_values)) {
  set.seed(123 + k_values[i])
  
  models[[i]] <- kmeans(
    x_train_scaled,
    centers = k_values[i],
    nstart = 50,
    iter.max = 100
  )
  
  validation$total_within_ss[i] <- models[[i]]$tot.withinss
  
  validation$mean_silhouette[i] <- safe_silhouette(
    models[[i]]$cluster[silhouette_index],
    silhouette_distances
  )
}

if (all(is.na(validation$mean_silhouette))) {
  stop("Could not calculate valid training silhouette scores.")
}

best_index <- which.max(
  replace(
    validation$mean_silhouette,
    is.na(validation$mean_silhouette),
    -Inf
  )
)

best_k <- validation$k[best_index]
final_model <- models[[best_index]]

print(validation)
cat("\nSelected number of clusters:", best_k, "\n")

# Elbow plot.
elbow_plot <- ggplot(
  validation,
  aes(x = k, y = total_within_ss)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = k_values) +
  labs(
    title = "Elbow Method",
    x = "Number of clusters (k)",
    y = "Total within-cluster sum of squares"
  ) +
  theme_minimal()

# Silhouette selection plot.
silhouette_plot <- ggplot(
  validation,
  aes(x = k, y = mean_silhouette)
) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(xintercept = best_k, linetype = "dashed") +
  scale_x_continuous(breaks = k_values) +
  labs(
    title = "Training Silhouette Scores",
    subtitle = "Higher values indicate better-separated clusters",
    x = "Number of clusters (k)",
    y = "Mean silhouette width"
  ) +
  theme_minimal()

print(elbow_plot)
print(silhouette_plot)

# Assign held-out facilities to the nearest TRAINING centroid.
assign_nearest <- function(x, centers) {
  squared_distances <- vapply(
    seq_len(nrow(centers)),
    function(j) {
      rowSums(sweep(x, 2, centers[j, ], "-")^2)
    },
    numeric(nrow(x))
  )
  
  max.col(-squared_distances, ties.method = "first")
}

train_cluster <- final_model$cluster
test_cluster <- assign_nearest(x_test_scaled, final_model$centers)

# Test-set silhouette: descriptive separation among held-out
# facilities, using assignments from the training model.
test_silhouette_index <- sample(
  seq_len(nrow(x_test_scaled)),
  min(2000L, nrow(x_test_scaled))
)

test_silhouette <- safe_silhouette(
  test_cluster[test_silhouette_index],
  dist(x_test_scaled[test_silhouette_index, , drop = FALSE])
)

# Average squared distance to assigned training centroid.
test_centroid_mse <- mean(
  rowSums(
    (
      x_test_scaled -
        final_model$centers[test_cluster, , drop = FALSE]
    )^2
  )
)

metrics <- tibble(
  selected_k = best_k,
  training_facilities = nrow(train),
  test_facilities = nrow(test),
  training_mean_silhouette =
    validation$mean_silhouette[best_index],
  test_mean_silhouette = test_silhouette,
  training_mean_squared_centroid_distance =
    final_model$tot.withinss / nrow(train),
  test_mean_squared_centroid_distance = test_centroid_mse
)

print(metrics)

# ------------------------------------------------------------
# 4. INTERPRET AND EXPORT RESULTS
# ------------------------------------------------------------

assignments <- bind_rows(
  train %>% mutate(dataset = "Train", cluster = train_cluster),
  test %>% mutate(dataset = "Test", cluster = test_cluster)
) %>%
  arrange(cluster, `DEC ID`)

# Descriptive cluster summaries in ORIGINAL tons.
# These use the observed facility averages, not imputed values.
cluster_profiles <- assignments %>%
  group_by(cluster) %>%
  summarise(
    facilities = n(),
    across(
      all_of(pollutants),
      list(
        mean = ~ mean_or_na(.x),
        median = ~ if (all(is.na(.x))) {
          NA_real_
        } else {
          median(.x, na.rm = TRUE)
        }
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

print(cluster_profiles, width = Inf)

# Heatmap of model centers in standardized log-emissions space.
center_data <- as_tibble(
  final_model$centers,
  .name_repair = "minimal"
) %>%
  mutate(cluster = factor(seq_len(nrow(final_model$centers)))) %>%
  pivot_longer(
    cols = all_of(used_pollutants),
    names_to = "pollutant",
    values_to = "standardized_log_emissions"
  )

profile_plot <- ggplot(
  center_data,
  aes(x = pollutant, y = cluster,
      fill = standardized_log_emissions)
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0
  ) +
  labs(
    title = "Pollution Profiles of Facility Clusters",
    x = NULL,
    y = "Cluster",
    fill = "Scaled\nlog emissions"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(profile_plot)

write_csv(validation, "cluster_k_validation.csv")
write_csv(metrics, "clustering_validation_metrics.csv")
write_csv(assignments, "facility_cluster_assignments.csv")
write_csv(cluster_profiles, "cluster_pollution_profiles.csv")

ggsave("elbow_plot.png", elbow_plot, width = 7, height = 5)
ggsave("silhouette_plot.png", silhouette_plot, width = 7, height = 5)
ggsave("cluster_profiles.png", profile_plot, width = 10, height = 5)

# Save the fitted model and preprocessing parameters together.
saveRDS(
  list(
    model = final_model,
    all_pollutants = pollutants,
    used_pollutants = used_pollutants,
    imputation_medians = train_medians,
    training_mean = training_mean,
    training_sd = training_sd,
    selected_k = best_k,
    validation = metrics
  ),
  "facility_pollution_clustering_model.rds"
)