library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(ggbeeswarm)
library(ggcorrplot)

# -----------------------------
# Paths
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/proportion_output"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/summary"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Sample directories
# -----------------------------
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV_RppH", sample_dirs)]

get_timepoint <- function(name) {
  if (grepl("RppH_5h", name)) return("5h")
  if (grepl("RppH_", name)) return("24h")
  return("24h") # Defaulting to 24h as per your previous logic
}

# -----------------------------
# Read all replicate data
# -----------------------------
all_data <- bind_rows(lapply(sample_dirs, function(sample_dir) {
  files <- list.files(sample_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  folder_id <- basename(sample_dir)
  
  bind_rows(lapply(files, function(f) {
    df <- read.delim(f, check.names = FALSE)
    # Unique ID per file
    df$Replicate <- paste(folder_id, basename(f), sep = "_")
    df$Timepoint <- get_timepoint(folder_id)
    df
  }))
}))

# -----------------------------
# Add Phase & CPM
# -----------------------------
all_data <- all_data %>%
  mutate(
    Phase = sub(".*_", "", Gene_Timepoint),
    Phase = factor(Phase, levels = c("early", "intermediate", "late"))
  ) %>%
  filter(!is.na(Phase)) %>%
  mutate(across(c(Natural, Unnatural, polyA, micropolyA), ~ replace_na(., 0))) %>%
  rowwise() %>%
  mutate(
    TotalReads = sum(Natural, Unnatural, polyA, micropolyA),
    polyA_CPM = ifelse(TotalReads > 0, polyA / TotalReads * 1e6, 0)
  ) %>%
  ungroup()

# -----------------------------
# 1. Standardized Correlation Heatmaps (Fixed 7x7 Grid)
# -----------------------------
# Create a Master list of all 7 replicates so the order never changes
master_replicate_list <- sort(unique(all_data$Replicate))

unique_phases <- levels(all_data$Phase)
pdf(file.path(plot_base, "Replicate_Correlation_by_Phase_Standardized.pdf"), width = 11, height = 9)

for (ph in unique_phases) {
  
  # Filter for phase
  phase_data <- all_data %>%
    filter(Phase == ph) %>%
    group_by(Gene_Timepoint, Replicate) %>%
    summarise(polyA_CPM = mean(polyA_CPM), .groups = "drop")
  
  # Pivot and FORCE all master replicates to be present
  phase_wide <- phase_data %>%
    pivot_wider(names_from = Replicate, values_from = polyA_CPM, values_fill = 0)
  
  # Add columns for any of the 7 replicates missing in this specific phase
  missing_reps <- setdiff(master_replicate_list, colnames(phase_wide))
  for(m in missing_reps) { phase_wide[[m]] <- 0 }
  
  # Select columns in FIXED order and remove Gene names for matrix
  cor_input <- phase_wide %>%
    select(all_of(master_replicate_list))
  
  # Compute Spearman
  cor_matrix_ph <- cor(cor_input, method = "spearman")
  cor_matrix_ph[is.na(cor_matrix_ph)] <- 0 # Handle zero-count samples
  
  # Plot with hc.order = FALSE to keep our fixed order
  p <- ggcorrplot(cor_matrix_ph, 
                  hc.order = FALSE, 
                  type = "lower", 
                  lab = TRUE, 
                  lab_size = 2.5,
                  colors = c("#6D9ECB", "white", "#E46726"),
                  title = paste("Spearman Correlation:", tools::toTitleCase(ph), "Phase")) +
    labs(subtitle = paste("Fixed Sample Order | Genes analyzed:", nrow(phase_wide))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))
  
  print(p)
}
dev.off()

# -----------------------------
# 2. Variability (CV) Boxplot
# -----------------------------
cv_data <- all_data %>%
  group_by(Gene_Timepoint, Phase, Timepoint) %>%
  summarise(
    mean_cpm = mean(polyA_CPM),
    sd_cpm = sd(polyA_CPM),
    .groups = "drop"
  ) %>%
  mutate(CV = ifelse(mean_cpm > 0, sd_cpm / mean_cpm, 0))

pdf(file.path(plot_base, "Expression_Variability_CV.pdf"), width = 8, height = 6)
p_cv <- ggplot(cv_data, aes(x = Phase, y = CV, fill = Phase)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_quasirandom(size = 1, alpha = 0.6) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  theme_bw() +
  labs(title = "Gene Expression Variability by Phase", 
       y = "Coefficient of Variation (SD/Mean)", x = "Viral Phase")
print(p_cv)
dev.off()





library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(ggbeeswarm)
library(ggcorrplot)

# -----------------------------
# Paths
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/proportion_output"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/summary"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV_RppH", sample_dirs)]

get_timepoint <- function(name) {
  if (grepl("RppH_5h", name)) return("5h")
  return("24h") 
}

# -----------------------------
# Read all replicate data
# -----------------------------
all_data <- bind_rows(lapply(sample_dirs, function(sample_dir) {
  files <- list.files(sample_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  folder_id <- basename(sample_dir)
  
  bind_rows(lapply(files, function(f) {
    df <- read.delim(f, check.names = FALSE)
    df$Replicate <- paste(folder_id, basename(f), sep = "_")
    df$Timepoint <- get_timepoint(folder_id)
    df
  }))
}))

# -----------------------------
# Add Phase & micropolyA CPM (UPDATED)
# -----------------------------
all_data <- all_data %>%
  mutate(
    Phase = sub(".*_", "", Gene_Timepoint),
    Phase = factor(Phase, levels = c("early", "intermediate", "late"))
  ) %>%
  filter(!is.na(Phase)) %>%
  mutate(across(c(Natural, Unnatural, polyA, micropolyA), ~ replace_na(., 0))) %>%
  rowwise() %>%
  mutate(
    TotalReads = sum(Natural, Unnatural, polyA, micropolyA),
    # FOCUS ON MICROPOLYA HERE
    micropolyA_CPM = ifelse(TotalReads > 0, micropolyA / TotalReads * 1e6, 0)
  ) %>%
  ungroup()

# -----------------------------
# 1. Standardized Correlation Heatmaps (micropolyA)
# -----------------------------
master_replicate_list <- sort(unique(all_data$Replicate))
unique_phases <- levels(all_data$Phase)

pdf(file.path(plot_base, "Replicate_Correlation_by_Phase_micropolyA.pdf"), width = 11, height = 9)

for (ph in unique_phases) {
  
  phase_data <- all_data %>%
    filter(Phase == ph) %>%
    group_by(Gene_Timepoint, Replicate) %>%
    summarise(micropolyA_CPM = mean(micropolyA_CPM), .groups = "drop")
  
  phase_wide <- phase_data %>%
    pivot_wider(names_from = Replicate, values_from = micropolyA_CPM, values_fill = 0)
  
  missing_reps <- setdiff(master_replicate_list, colnames(phase_wide))
  for(m in missing_reps) { phase_wide[[m]] <- 0 }
  
  cor_input <- phase_wide %>%
    select(all_of(master_replicate_list))
  
  # Spearman is excellent for micropolyA which can be sparse/low count
  cor_matrix_ph <- cor(cor_input, method = "spearman")
  cor_matrix_ph[is.na(cor_matrix_ph)] <- 0 
  
  p <- ggcorrplot(cor_matrix_ph, 
                  hc.order = FALSE, 
                  type = "lower", 
                  lab = TRUE, 
                  lab_size = 2.5,
                  colors = c("#6D9ECB", "white", "#E46726"),
                  title = paste("Spearman Correlation (micropolyA):", tools::toTitleCase(ph), "Phase")) +
    labs(subtitle = paste("Fixed Sample Order | Genes analyzed:", nrow(phase_wide))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))
  
  print(p)
}
dev.off()

# -----------------------------
# 2. Variability (CV) Boxplot (micropolyA)
# -----------------------------
cv_data <- all_data %>%
  group_by(Gene_Timepoint, Phase, Timepoint) %>%
  summarise(
    mean_cpm = mean(micropolyA_CPM),
    sd_cpm = sd(micropolyA_CPM),
    .groups = "drop"
  ) %>%
  mutate(CV = ifelse(mean_cpm > 0, sd_cpm / mean_cpm, 0))

pdf(file.path(plot_base, "Expression_Variability_CV_micropolyA.pdf"), width = 8, height = 6)
p_cv <- ggplot(cv_data, aes(x = Phase, y = CV, fill = Phase)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_quasirandom(size = 1, alpha = 0.6) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  theme_bw() +
  labs(title = "micropolyA Expression Variability by Phase", 
       y = "Coefficient of Variation (SD/Mean)", x = "Viral Phase")
print(p_cv)
dev.off()
