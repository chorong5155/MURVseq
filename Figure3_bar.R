library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(ggbeeswarm)

# -----------------------------
# Paths
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/proportion_output"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV_RppH", sample_dirs)]

get_timepoint <- function(name) {
  if (grepl("5h", name)) return("5h")
  else {return("24h")}
}

# -----------------------------
# Read all replicate data
# -----------------------------
all_data <- bind_rows(lapply(sample_dirs, function(sample_dir) {
  files <- list.files(sample_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  bind_rows(lapply(files, function(f) {
    df <- read.delim(f, check.names = FALSE)
    df$Replicate <- basename(f)
    df$Timepoint <- get_timepoint(basename(sample_dir))
    df
  }))
}))

# -----------------------------
# Add Phase & clean counts
# -----------------------------
all_data <- all_data %>%
  mutate(
    Phase = sub(".*_", "", Gene_Timepoint),
    Phase = factor(Phase, levels = c("early", "intermediate", "late"))
  ) %>%
  filter(!is.na(Phase)) %>%
  mutate(across(
    c(Natural, Unnatural, polyA, micropolyA),
    ~ replace_na(., 0)
  ))

# Calculate CPM relative to total sample reads
all_data <- all_data %>%
  group_by(Replicate, Timepoint) %>%  
  mutate(
    TotalReads_sample = sum(Natural + Unnatural + polyA + micropolyA, na.rm = TRUE),
    Natural_CPM = ifelse(TotalReads_sample > 0, Natural / TotalReads_sample * 1e6, 0),
    # Note: Unnatural_CPM already includes micropolyA based on your logic below
    Unnatural_CPM = ifelse(TotalReads_sample > 0, (Unnatural + micropolyA) / TotalReads_sample * 1e6, 0),
    polyA_CPM = ifelse(TotalReads_sample > 0, polyA / TotalReads_sample * 1e6, 0)
  ) %>%
  ungroup()

genes_of_interest <- c(
  "VACWR025_late",
  "VACWR122_intermediate",
  "VACWR028_early",
  "VACWR103_early",
  "VACWR056_late",
  "VACWR070_intermediate",
  "VACWR150_intermediate"
)

# -----------------------------
# Calculate Percentages & Level Timepoints
# -----------------------------
all_data_bar <- all_data %>%
  filter(Gene_Timepoint %in% genes_of_interest) %>%
  # Explicitly level the Timepoint
  mutate(Timepoint = factor(Timepoint, levels = c("5h", "24h"))) %>%
  # Group by gene and replicate to find the total signal for that specific instance
  group_by(Gene_Timepoint, Replicate, Timepoint) %>%
  mutate(
    Gene_Total_CPM = Natural_CPM + Unnatural_CPM + polyA_CPM,
    Natural_Pct   = ifelse(Gene_Total_CPM > 0, (Natural_CPM / Gene_Total_CPM) * 100, 0),
    Unnatural_Pct = ifelse(Gene_Total_CPM > 0, (Unnatural_CPM / Gene_Total_CPM) * 100, 0),
    polyA_Pct     = ifelse(Gene_Total_CPM > 0, (polyA_CPM / Gene_Total_CPM) * 100, 0)
  ) %>%
  ungroup() %>%
  select(Gene_Timepoint, Phase, Timepoint, Replicate, Natural_Pct, Unnatural_Pct, polyA_Pct) %>%
  pivot_longer(
    cols = c(Natural_Pct, Unnatural_Pct, polyA_Pct),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    Category = factor(Category, 
                      levels = c("Natural_Pct", "Unnatural_Pct", "polyA_Pct"),
                      labels = c("Natural", "Unnatural", "polyA"))
  )

# -----------------------------
# Plotting
# -----------------------------
pdf(file.path(plot_base, "SelectedGenes_Percentage_boxplot.pdf"), width = 6, height = 5)

for (gene in genes_of_interest) {
  
  df_gene <- all_data_bar %>%
    filter(Gene_Timepoint == gene)
  
  p <- ggplot(df_gene, aes(x = Category, y = Percentage, fill = Category)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7,
      width = 0.6
    ) +
    # Switched to geom_quasirandom for better spacing of replicate points
    geom_quasirandom(
      aes(color = Category),
      width = 0.15,
      size = 2,
      alpha = 0.8
    ) +
    facet_wrap(~ Timepoint) +
    scale_fill_manual(values = wes_palette("Darjeeling1", 3)) +
    scale_color_manual(values = wes_palette("Darjeeling1", 3)) +
    theme_classic(base_size = 12) +
    labs(
      title = gene,
      subtitle = "Relative distribution of transcript types",
      x = NULL,
      y = "Percentage of Total Signal (%)",
      fill = "Category",
      color = "Category"
    ) +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    ylim(0, 100)
  
  print(p)
}

dev.off()



# -----------------------------
# Updated Genes of Interest for only poluA and natural
# -----------------------------
genes_of_interest <- c(
  "VACWR056_late",
  "VACWR070_intermediate",
  "VACWR129_late",
  "VACWR150_intermediate"
)

# -----------------------------
# Calculate Normalized Percentages (Natural + polyA only)
# -----------------------------
all_data_bar <- all_data %>%
  filter(Gene_Timepoint %in% genes_of_interest) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("5h", "24h"))) %>%
  group_by(Gene_Timepoint, Replicate, Timepoint) %>%
  mutate(
    # Create a new denominator excluding Unnatural
    Natural_polyA_Sum = Natural_CPM + polyA_CPM,
    # Normalize so these two equal 100%
    Natural_Normalized_Pct = ifelse(Natural_polyA_Sum > 0, (Natural_CPM / Natural_polyA_Sum) * 100, 0),
    polyA_Normalized_Pct   = ifelse(Natural_polyA_Sum > 0, (polyA_CPM / Natural_polyA_Sum) * 100, 0)
  ) %>%
  ungroup() %>%
  select(Gene_Timepoint, Phase, Timepoint, Replicate, Natural_Normalized_Pct, polyA_Normalized_Pct) %>%
  pivot_longer(
    cols = c(Natural_Normalized_Pct, polyA_Normalized_Pct),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    Category = factor(Category, 
                      levels = c("Natural_Normalized_Pct", "polyA_Normalized_Pct"),
                      labels = c("Natural", "polyA"))
  )

# -----------------------------
# 1. Normalized Calculation
# -----------------------------
all_data_bar <- all_data %>%
  filter(Gene_Timepoint %in% genes_of_interest) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("5h", "24h"))) %>%
  group_by(Gene_Timepoint, Replicate, Timepoint) %>%
  mutate(
    # New 100% total excluding Unnatural
    Natural_polyA_Sum = Natural_CPM + polyA_CPM,
    Natural_Normalized_Pct = ifelse(Natural_polyA_Sum > 0, (Natural_CPM / Natural_polyA_Sum) * 100, 0),
    polyA_Normalized_Pct   = ifelse(Natural_polyA_Sum > 0, (polyA_CPM / Natural_polyA_Sum) * 100, 0)
  ) %>%
  ungroup() %>%
  select(Gene_Timepoint, Timepoint, Replicate, Natural_Normalized_Pct, polyA_Normalized_Pct) %>%
  pivot_longer(
    cols = c(Natural_Normalized_Pct, polyA_Normalized_Pct),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    # Natural = 1 (Green), polyA = 2 (Yellow) in Darjeeling1
    Category = factor(Category, 
                      levels = c("Natural_Normalized_Pct", "polyA_Normalized_Pct"),
                      labels = c("Natural", "polyA"))
  )

# -----------------------------
# 2. Plotting with Yellow polyA
# -----------------------------
# WesAnderson Darjeeling1 colors: 1=#FF0000(Red), 2=#00A08A(Teal), 3=#F2AD00(Yellow), 4=#F98400(Orange), 5=#5BBCD6(Blue)
# We pick colors 2 (Teal/Green) and 3 (Yellow) for the 2-category comparison
my_colors <- wes_palette("Darjeeling1", 5)[c(1, 3)] 

pdf(file.path(plot_base, "Normalized_Yellow_polyA_boxplot.pdf"), width = 6, height = 5)

for (gene in genes_of_interest) {
  
  df_gene <- all_data_bar %>%
    filter(Gene_Timepoint == gene)
  
  p <- ggplot(df_gene, aes(x = Category, y = Percentage, fill = Category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_quasirandom(aes(color = Category), width = 0.15, size = 2, alpha = 0.8) +
    facet_wrap(~ Timepoint) +
    # This applies the colors: Natural gets Greenish-Teal, polyA gets Yellow
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    theme_classic(base_size = 12) +
    labs(
      title = gene,
      subtitle = "Natural vs polyA (100% Total)",
      x = NULL,
      y = "Normalized Percentage (%)"
    ) +
    theme(legend.position = "top", strip.background = element_blank()) +
    ylim(0, 100)
  
  print(p)
}

dev.off()

