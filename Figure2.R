#Figure2A&sup
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)

# -----------------------------
# Paths
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/proportion_output"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/proportion"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Samples
# -----------------------------
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
samples <- basename(sample_dirs)

get_timepoint <- function(name) {
  if (grepl("RppH_5h", name)) return("5h")
  if (grepl("RppH_", name)) return("24h")
  return("unknown")
}

# -----------------------------
# Read all replicate data
# -----------------------------
all_data <- bind_rows(lapply(sample_dirs, function(sample_dir) {
  files <- list.files(sample_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  bind_rows(lapply(files, function(f) {
    d <- read.delim(f, check.names = FALSE)
    d$Replicate <- basename(f)
    d$Timepoint <- get_timepoint(basename(sample_dir))
    d
  }))
}))

# -----------------------------
# Merge micropolyA into Unnatural
# -----------------------------
all_data <- all_data %>%
  mutate(Unnatural = Unnatural + micropolyA) %>%
  select(-micropolyA)


convert_to_cpm <- function(df, count_cols = c("Natural","Unnatural","polyA")) {
  df <- df %>%
    mutate(TotalReads = rowSums(select(., all_of(count_cols)), na.rm = TRUE)) %>%
    mutate(across(all_of(count_cols), ~ . / sum(.) * 1e6))
  return(df)
}

# -----------------------------
# Function to plot CPM stacked bars ordered by Phase and CPM
# -----------------------------
plot_cpm_combined <- function(df, timepoint, output_dir) {
  
  # Filter by timepoint
  df_tp <- df %>% filter(Timepoint == timepoint)
  
  # Sum replicates per gene
  df_sum <- df_tp %>%
    group_by(Gene_Timepoint, Phase = sub(".*_", "", Gene_Timepoint)) %>%
    summarise(
      Natural = sum(Natural, na.rm = TRUE),
      Unnatural = sum(Unnatural, na.rm = TRUE),
      polyA = sum(polyA, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_sum$Phase <- factor(df_sum$Phase, levels = c("early","intermediate","late","unknown"))
  
  # Convert to CPM
  df_sum <- df_sum %>%
    mutate(TotalGeneCounts = Natural + Unnatural + polyA) %>%
    mutate(across(c(Natural, Unnatural, polyA), ~ . / TotalGeneCounts * 1e6))
  
  # Order genes within phase by total CPM
  #df_sum <- df_sum %>%
  #arrange(Phase, desc(Natural + Unnatural + polyA)) %>%
  #mutate(Gene_Timepoint = factor(Gene_Timepoint, levels = unique(Gene_Timepoint)))
  
  # Pivot longer for ggplot
  df_long <- df_sum %>%
    pivot_longer(cols = c(Natural, Unnatural, polyA),
                 names_to = "Category", values_to = "CPM") %>%
    mutate(Category = factor(Category, levels = c("Natural","Unnatural","polyA")))
  
  # Plot
  p <- ggplot(df_long, aes(x = Gene_Timepoint, y = CPM, fill = Category)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = wes_palette("Darjeeling1", 3, type = "discrete")) +
    theme_classic(base_size = 10) +
    labs(title = paste0("CPM per Gene (", timepoint, ")"),
         x = "Gene & Timepoint",
         y = "Counts Per Million (CPM)",
         fill = "UTR Type") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          legend.position = "top")
  
  # Save
  ggsave(file.path(output_dir, paste0("CPM_stacked_", timepoint, ".pdf")),
         plot = p, width = 20, height = 4)
  
  print(p)
}

# -----------------------------
# Generate plots for 5h and 24h
# -----------------------------
plot_cpm_combined(all_data, "5h", plot_base)
plot_cpm_combined(all_data, "24h", plot_base)



plot_CPM_combined <- function(df, timepoint, output_dir) {
  
  df_tp <- df %>% filter(Timepoint == timepoint)
  
  summary_df <- df_tp %>%
    group_by(Gene_Timepoint) %>%
    summarise(
      Natural   = sum(Natural, na.rm = TRUE),
      Unnatural = sum(Unnatural, na.rm = TRUE),
      polyA     = sum(polyA, na.rm = TRUE),
      .groups = "drop"
    )
  
  total_counts <- sum(summary_df$Natural +
                        summary_df$Unnatural +
                        summary_df$polyA)
  
  summary_df <- summary_df %>%
    mutate(
      across(c(Natural, Unnatural, polyA), ~ . / total_counts * 1e6),
      Phase = sub(".*_", "", Gene_Timepoint),
      Phase = factor(Phase, levels = c("early", "intermediate", "late")),
      total_CPM = Natural + Unnatural + polyA
    ) %>%
    arrange(Phase, desc(total_CPM))
  
  data_long <- summary_df %>%
    pivot_longer(
      cols = c(Natural, Unnatural, polyA),
      names_to = "Category",
      values_to = "CPM"
    ) %>%
    mutate(
      Category = factor(Category, levels = c("Natural","Unnatural","polyA")),
      Gene_Timepoint = factor(Gene_Timepoint, levels = unique(Gene_Timepoint))
    )
  
  p <- ggplot(data_long, aes(x = Gene_Timepoint, y = CPM, fill = Category)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = wesanderson::wes_palette("Darjeeling1", 3)) +
    theme_classic(base_size = 10) +
    labs(
      title = paste0("Combined UTR Counts per Gene in CPM (", timepoint, ")"),
      x = "Gene",
      y = "Counts Per Million (CPM)",
      fill = "UTR Type"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      legend.position = "top"
    )
  
  ggsave(
    file.path(output_dir, paste0("UTR_CPM_ordered_", timepoint, ".pdf")),
    plot = p, width = 20, height = 4
  )
  
  print(p)
}

# -----------------------------
# Generate CPM plots for 5h and 24h
# -----------------------------
plot_CPM_combined(all_data, "5h", plot_base)
plot_CPM_combined(all_data, "24h", plot_base)


#Figure2B
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(wesanderson)

# -----------------------------
# Paths and settings
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/proportion_output"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/summary"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

phase_colors <- c(early = "skyblue", intermediate = "orange", late = "salmon")

get_timepoint <- function(name) {
  if (grepl("VacV_RppH_5h", name)) return("5h")
  if (grepl("VacV_RppH", name)) return("24h")
  return(NA)
}

sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV_RppH", sample_dirs)]

# Pairwise comparisons for Wilcoxon test (same as before)
comparisons <- list(
  c("early", "intermediate"),
  c("early", "late"),
  c("intermediate", "late")
)

# -----------------------------
# Read and combine all replicates
# -----------------------------
all_data <- lapply(sample_dirs, function(sample_dir) {
  replicate_name <- basename(sample_dir)
  timepoint <- get_timepoint(replicate_name)
  if (is.na(timepoint)) return(NULL)
  
  f <- file.path(sample_dir, "combined_wide.txt")
  if (!file.exists(f)) return(NULL)
  
  df <- read.delim(f, check.names = FALSE)
  if (nrow(df) == 0) return(NULL)
  
  df$Replicate <- replicate_name
  df$Timepoint <- timepoint
  
  # Assign Phase
  df$Phase <- ifelse(grepl("early", df$Gene_Timepoint), "early",
                     ifelse(grepl("intermediate", df$Gene_Timepoint), "intermediate",
                            ifelse(grepl("late", df$Gene_Timepoint), "late", NA)))
  df$Phase <- factor(df$Phase, levels = c("early", "intermediate", "late"))
  
  df <- df %>% filter(!is.na(Phase)) %>% droplevels()
  
  df <- df %>% mutate(across(c(Natural, Unnatural, polyA, micropolyA), ~replace_na(.,0)))
  
  df <- df %>%
    rowwise() %>%
    mutate(
      TotalReads = sum(Natural, Unnatural, polyA, micropolyA),
      polyA_CPM = ifelse(TotalReads > 0, polyA / TotalReads * 1e6, 0),
      micropolyA_CPM = ifelse(TotalReads > 0, micropolyA / TotalReads * 1e6, 0)
    ) %>%
    ungroup()
  
  df_long_CPM <- df %>%
    select(Replicate, Phase, Timepoint, Gene_Timepoint, polyA_CPM, micropolyA_CPM) %>%
    pivot_longer(cols = c(polyA_CPM, micropolyA_CPM),
                 names_to = "Category",
                 values_to = "Value") %>%
    mutate(Category = factor(Category,
                             levels = c("polyA_CPM", "micropolyA_CPM"),
                             labels = c("polyA","micropolyA")),
           Type = "CPM")
  
  df_long_pct <- df_long_CPM %>%
    group_by(Replicate, Phase, Gene_Timepoint) %>%
    mutate(Value = Value / sum(Value, na.rm = TRUE) * 100) %>%
    ungroup() %>%
    mutate(Type = "Percent")
  
  bind_rows(df_long_pct)
}) %>% bind_rows()

# -----------------------------
# Compute median per Phase × Category × Type
# -----------------------------
df_median <- all_data %>%
  group_by(Phase, Category, Type) %>%
  summarize(Median = median(Value, na.rm = TRUE), .groups = "drop")
write.csv(df_median, file.path(plot_base, "medians_combined.csv"), row.names = FALSE)

# -----------------------------
# Plot all genes with median + Wilcoxon test
# -----------------------------
p <- ggplot(all_data, aes(x = Phase, y = Value, fill = Phase)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE) +
  geom_quasirandom(aes(color = Phase), width = 0.25, size = 1.5, alpha = 0.7) +
  geom_point(data = df_median, aes(x = Phase, y = Median), color = "black", size = 2) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    label.y.npc = 0.95,
    size = 3
  ) +
  facet_grid(Type ~ Category + Timepoint, scales = "free_y") +
  scale_fill_manual(values = phase_colors) +
  scale_color_manual(values = phase_colors) +
  theme_classic(base_size = 12) +
  labs(
    title = "polyA vs micropolyA – Combined replicates",
    x = "Phase",
    y = "Percent per phase",
    fill = "Phase",
    color = "Phase"
  ) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 11)
  )

pdf_file <- file.path(plot_base, "polyA_vs_micropolyA_combined_all_replicates.pdf")
ggsave(filename = pdf_file, plot = p, width = 8, height = 4.5)

message("Saved combined plot for all replicates")


library(dplyr)
library(tidyr)

# Pairwise phase comparisons
phase_comparisons <- list(
  c("early", "intermediate"),
  c("early", "late"),
  c("intermediate", "late")
)

# Initialize empty table
p_values_phase <- data.frame()

# Loop over all Timepoint × Category × Type
for(tp in unique(all_data$Timepoint)) {
  for(cat in unique(all_data$Category)) {
    for(tp_type in unique(all_data$Type)) {
      
      df_sub <- all_data %>%
        filter(Timepoint == tp, Category == cat, Type == tp_type)
      
      for(comp in phase_comparisons) {
        ph1 <- df_sub %>% filter(Phase == comp[1]) %>% pull(Value)
        ph2 <- df_sub %>% filter(Phase == comp[2]) %>% pull(Value)
        
        pval <- if(length(ph1) > 0 & length(ph2) > 0) {
          wilcox.test(ph1, ph2)$p.value
        } else {
          NA
        }
        
        p_values_phase <- rbind(p_values_phase, data.frame(
          Timepoint = tp,
          Category = cat,
          Type = tp_type,
          Comparison = paste(comp, collapse = "_vs_"),
          p_value = pval
        ))
      }
    }
  }
}

# View and save
print(p_values_phase)
write.csv(p_values_phase, file.path(plot_base, "wilcox_pvalues_phase_combined.csv"), row.names = FALSE)


library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(purrr)
library(scales)

# -----------------------------
# Paths
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/micropolyA_processed"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/summary"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Phases and colors
# -----------------------------
phases <- c("early", "intermediate", "late")
phase_colors <- c(
  early = "skyblue",
  intermediate = "orange",
  late = "salmon"
)


# -----------------------------
# Timepoint extractor & Leveling
# -----------------------------
get_timepoint <- function(name) {
  if (grepl("VacV_RppH_5h", name)) return("5h")
  if (grepl("VacV_RppH", name)) return("24h")
  return(NA)
}

# -----------------------------
# Read all data
# -----------------------------
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)

all_data <- map_dfr(sample_dirs, function(sample_dir) {
  sample_name <- basename(sample_dir)
  timepoint <- get_timepoint(sample_name)
  if (is.na(timepoint)) return(NULL)
  
  gene_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
  
  map_dfr(gene_dirs, function(gene_dir) {
    gene_name <- basename(gene_dir)
    tsv_file <- list.files(gene_dir, pattern = "_full\\.tsv$", full.names = TRUE)
    if (length(tsv_file) == 0) return(NULL)
    
    df <- read_tsv(tsv_file, col_types = cols(.default = "c"), show_col_types = FALSE)
    
    df %>%
      mutate(
        Length = as.numeric(polyA_length),
        Gene = gene_name,
        Phase = case_when(
          grepl("early", gene_name, ignore.case = TRUE) ~ "early",
          grepl("intermediate", gene_name, ignore.case = TRUE) ~ "intermediate",
          grepl("late", gene_name, ignore.case = TRUE) ~ "late",
          TRUE ~ NA_character_
        ),
        Sample = sample_name,
        Timepoint = timepoint
      ) %>%
      filter(!is.na(Length), !is.na(Phase)) %>%
      group_by(Sample, Timepoint, Gene, Phase, Length) %>%
      summarise(Count = n(), .groups = "drop")
  })
})

# -----------------------------
# PER-REPLICATE NORMALIZATION (The Fix)
# -----------------------------
all_data <- all_data %>%
  group_by(Sample) %>% # Grouping by Sample ensures CPM is relative to each library
  mutate(CPM = Count / sum(Count) * 1e6) %>%
  ungroup() %>%
  mutate(
    Phase = factor(Phase, levels = phases),
    Timepoint = factor(Timepoint, levels = c("5h", "24h")) # Leveling timepoints
  )

# -----------------------------
# Summarise CPM per sample & gene (keeps Replicates separate)
# -----------------------------
gene_sample_summary <- all_data %>%
  group_by(Sample, Timepoint, Gene, Phase) %>%
  summarise(
    Total_CPM = sum(CPM),
    Dominant_Length = Length[which.max(Count)],
    .groups = "drop"
  )

# -----------------------------
# Calculate Stats for Ordering
# -----------------------------
gene_timepoint_summary <- gene_sample_summary %>%
  group_by(Timepoint, Gene, Phase) %>%
  summarise(
    Median_CPM = median(Total_CPM),
    .groups = "drop"
  )

# -----------------------------
# Plotting Loop
# -----------------------------
topN <- 10

# Use levels(gene_timepoint_summary$Timepoint) to ensure 5h is processed first
for (tp in levels(gene_timepoint_summary$Timepoint)) {
  
  # Filter for current timepoint
  tp_data <- gene_timepoint_summary %>% filter(Timepoint == tp)
  
  if(nrow(tp_data) == 0) next
  
  # Select Top Genes by Median
  top_genes <- tp_data %>%
    arrange(desc(Median_CPM)) %>%
    slice_head(n = topN) %>%
    mutate(Gene_Phase = paste(Gene, Phase, sep = "_"))
  
  gene_order <- rev(top_genes$Gene_Phase)
  
  # Prepare points for this specific timepoint
  plot_points <- gene_sample_summary %>%
    filter(Timepoint == tp) %>%
    semi_join(top_genes, by = c("Gene", "Phase")) %>%
    mutate(
      Gene_Phase = paste(Gene, Phase, sep = "_"),
      Gene_Phase = factor(Gene_Phase, levels = gene_order)
    )
  
  # Summary for labels
  plot_summary <- plot_points %>%
    group_by(Gene_Phase, Phase) %>%
    summarise(
      Median_CPM = median(Total_CPM),
      # Mode of dominant length across replicates
      Dominant_Length = names(sort(table(Dominant_Length), decreasing = TRUE))[1],
      .groups = "drop"
    )
  
  p <- ggplot(plot_points, aes(x = Gene_Phase, y = Total_CPM, fill = Phase)) +
    geom_boxplot(width = 0.7, alpha = 0.5, outlier.shape = NA, color = "black") +
    # Use ggbeeswarm for better replicate visualization
    ggbeeswarm::geom_quasirandom(width = 0.15, size = 2, shape = 21, color = "black", alpha = 0.9) +
    geom_text(
      data = plot_summary,
      aes(y = Median_CPM, label = Dominant_Length),
      hjust = -0.5, # Move text to the right of the median line
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(values = phase_colors) +
    coord_flip() +
    labs(
      title = paste0("Top ", topN, " Genes by Median CPM (", tp, ")"),
      subtitle = "Each point represents an individual replicate",
      x = "Gene (Phase)",
      y = "CPM (Normalized per Sample)",
      fill = "Viral Phase"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  # Save PDF
  pdf(file.path(plot_base, paste0("nonpure_Top", topN, "_", tp, ".pdf")), width = 8, height = 6)
  print(p)
  dev.off()
  
  message("Saved plot for time point: ", tp)
}


library(dplyr)
library(ggplot2)
library(readr)

# -----------------------------
# Paths
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/summary"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Phases and colors
# -----------------------------
phases <- c("early", "intermediate", "late")
phase_colors <- c(early = "skyblue", intermediate = "orange", late = "salmon")

# -----------------------------
# Timepoint extractor
# -----------------------------
# ... [Keep initial paths and color setups] ...

# -----------------------------
# Timepoint extractor & Leveling
# -----------------------------
get_timepoint <- function(name) {
  if (grepl("VacV_RppH_5h", name)) return("5h")
  if (grepl("VacV_RppH", name)) return("24h")
  return(NA)
}

# -----------------------------
# Read all data
# -----------------------------
all_data <- data.frame()
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  timepoint <- get_timepoint(sample_name)
  if (is.na(timepoint)) next
  
  polyA_folder <- file.path(sample_dir, "polyA_length_per_gene")
  if (!dir.exists(polyA_folder)) next
  
  files <- list.files(polyA_folder, pattern = "_polyA_lengths.txt$", full.names = TRUE)
  if (length(files) == 0) next
  
  for (f in files) {
    lines <- readLines(f)
    idx <- grep("leader length counts:", lines, ignore.case = TRUE)
    if (length(idx) == 0 || idx >= length(lines)) next
    data_lines <- lines[(idx + 1):length(lines)]
    data_lines <- data_lines[nzchar(data_lines)]
    
    df <- tryCatch(read.table(text = data_lines, header = FALSE, col.names = c("Count", "Length")),
                   error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    
    header <- lines[1]
    gene_match <- regexec("Poly-A leader length summary for (VACWR[0-9]+) \\((early|intermediate|late)\\)", header)
    match_data <- regmatches(header, gene_match)[[1]]
    if (length(match_data) != 3) next
    
    df$Gene <- match_data[2]
    df$Phase <- match_data[3]
    df$Sample <- sample_name
    df$Timepoint <- timepoint
    
    all_data <- bind_rows(all_data, df)
  }
}

# -----------------------------
# PER-REPLICATE NORMALIZATION
# -----------------------------
all_data <- all_data %>%
  group_by(Sample) %>% 
  mutate(CPM = Count / sum(Count, na.rm = TRUE) * 1e6) %>%
  ungroup() %>%
  mutate(
    Phase = factor(Phase, levels = phases),
    Timepoint = factor(Timepoint, levels = c("5h", "24h")) # Correct leveling
  )

# -----------------------------
# Summarise CPM per Sample (keeps Replicates separate)
# -----------------------------
gene_sample_summary <- all_data %>%
  group_by(Sample, Timepoint, Gene, Phase) %>%
  summarise(
    Total_CPM = sum(CPM),
    Dominant_Length = Length[which.max(Count)],
    .groups = "drop"
  )

# -----------------------------
# Stats for Ordering Only
# -----------------------------
gene_timepoint_summary <- gene_sample_summary %>%
  group_by(Timepoint, Gene, Phase) %>%
  summarise(
    Median_CPM = median(Total_CPM),
    .groups = "drop"
  )

# -----------------------------
# Plotting Loop
# -----------------------------
topN <- 10

# Iterating through factor levels ensures 5h is processed first
for (tp in levels(gene_timepoint_summary$Timepoint)) {
  
  tp_data <- gene_timepoint_summary %>% filter(Timepoint == tp)
  if(nrow(tp_data) == 0) next
  
  # Select Top Genes by Median for this timepoint
  top_genes <- tp_data %>%
    arrange(desc(Median_CPM)) %>%
    slice_head(n = topN) %>%
    mutate(Gene_Phase = paste(Gene, Phase, sep = "_"))
  
  gene_order <- rev(top_genes$Gene_Phase)
  
  # Prepare plotting data
  plot_points <- gene_sample_summary %>%
    filter(Timepoint == tp) %>%
    semi_join(top_genes, by = c("Gene", "Phase")) %>%
    mutate(
      Gene_Phase = paste(Gene, Phase, sep = "_"),
      Gene_Phase = factor(Gene_Phase, levels = gene_order)
    )
  
  # Summary for labels (Median across replicates)
  plot_summary <- plot_points %>%
    group_by(Gene_Phase, Phase) %>%
    summarise(
      Median_Value = median(Total_CPM),
      # Calculate the most common dominant length across replicates
      Dominant_Length = names(sort(table(Dominant_Length), decreasing = TRUE))[1],
      .groups = "drop"
    )
  
  p <- ggplot(plot_points, aes(x = Gene_Phase, y = Total_CPM, fill = Phase)) +
    geom_boxplot(width = 0.7, alpha = 0.5, outlier.shape = NA, color = "black") +
    # Use quasirandom dots for a professional look
    ggbeeswarm::geom_quasirandom(width = 0.15, size = 2.5, shape = 21, color = "black", alpha = 0.9) +
    # Add Dominant Length label near the median
    geom_text(
      data = plot_summary, 
      aes(y = Median_Value, label = Dominant_Length),
      hjust = -0.6, 
      size = 3.5, 
      fontface = "bold"
    ) +
    scale_fill_manual(values = phase_colors) +
    coord_flip() +
    labs(
      title = paste0("Top ", topN, " Genes by PolyA Leader CPM (", tp, ")"),
      subtitle = "Normalized per replicate; labels indicate dominant leader length",
      x = "Gene (Phase)",
      y = "CPM (Normalized per Sample)",
      fill = "Viral Phase"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  # Save PDF
  pdf(file.path(plot_base, paste0("polyA_Top", topN, "_Box_Dots_", tp, ".pdf")), width = 8, height = 6)
  print(p)
  dev.off()
  
  message("Saved PolyA box+dot plot for time point: ", tp)
}


