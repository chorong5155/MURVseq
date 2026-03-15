#Figure7B
library(ggplot2)
library(dplyr)
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/combined_TSS"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/TSS_combined"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

TSS_TYPES <- c("natural_TSS", "classicslippage_TSS", "slippage_TSS")

sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV_RppH", basename(sample_dirs))]

# Helper to get timepoint or replicate from sample name
get_timepoint <- function(name) {
  if (grepl("5h", name)) return("5h")
  else {return("24h")
  }
}
# Initialize all_data **before** sample loop
all_data <- list()

# =========================
# Process each sample
# =========================
for (sample_dir in sample_dirs) {
  
  sample_name <- basename(sample_dir)
  timepoint <- get_timepoint(sample_name)
  if (is.na(timepoint)) next
  
  gene_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
  
  for (gene_path in gene_dirs) {
    gene_name <- basename(gene_path)
    bed_files <- list.files(gene_path, full.names = TRUE)
    
    for (bed in bed_files) {
      TSS <- case_when(
        grepl("natural_TSS_counts", bed) ~ "natural_TSS",
        grepl("classicslippage_TSS_counts", bed) ~ "classicslippage_TSS",
        grepl("micro_TSS_counts", bed) ~ "slippage_TSS",
        TRUE ~ NA_character_
      )
      if (is.na(TSS)) next
      
      df <- tryCatch(
        read.table(bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1:4],
        error = function(e) NULL
      )
      if (is.null(df) || nrow(df) == 0) next
      colnames(df) <- c("chr", "start", "end", "count")
      
      df <- df[grepl("VacV", df$chr, ignore.case = TRUE), ]
      if (nrow(df) == 0) next
      
      df <- df %>%
        mutate(
          count = as.numeric(count),
          pos   = as.numeric(start),
          gene  = gene_name,
          type  = factor(TSS, levels = TSS_TYPES),
          sample = sample_name,
          timepoint = factor(timepoint, levels = c("5h", "24h"))
        ) %>%
        filter(count >= 100)  # filter low counts
      
      # Add to global list
      all_data[[length(all_data) + 1]] <- df
    }
  }
}

# =========================
# Combine all samples
# =========================
df_all <- bind_rows(all_data)
if (nrow(df_all) == 0) stop("No TSS data found after filtering!")

# =========================
# Plot per gene
# =========================
genes <- unique(df_all$gene)
for (gene_name in genes) {
  df_gene <- df_all %>% filter(gene == gene_name)
  if (nrow(df_gene) == 0) next
  df_gene <- df_gene %>% mutate(count_log = count + 1)
  p <- ggplot(df_gene, aes(x = pos, y = count_log, color = type, size = count_log)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(
      values = c(
        natural_TSS = "steelblue",
        classicslippage_TSS = "darkgreen",
        slippage_TSS = "firebrick"
      ),
      drop = FALSE
    ) +
    scale_size_continuous(range = c(1, 5)) +
    scale_y_log10() +
    labs(
      title = paste(gene_name),
      x = "Genomic position",
      y = "Read count (log10)",
      color = "TSS type",
      size = "Read count (log10)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.margin = margin(b = 20),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    facet_wrap(~timepoint, ncol = 1)
  
  out_pdf <- file.path(plot_base, paste0(gene_name, "_TSS_by_timepoint.pdf"))
  ggsave(out_pdf, plot = p, width = 8, height = 8 , units = "in")
  message("Saved: ", out_pdf)
}
#7D
library(tidyverse)

# -----------------------------
# Directories
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/uAUG"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/uAUG2"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# uAUG categories
# -----------------------------
category_colors <- c(
  "Category1_Productive_Upstream"    = "#0072B2",
  "Category2_Upstream_Stopped"       = "#D55E00",
  "Category3_Overlapping_OutOfFrame" = "#009E73",
  "Category4_No_uAUG"                = "#CC79A7"
)

# -----------------------------
# Only VacV samples
# -----------------------------
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV", sample_dirs)]

if (length(sample_dirs) == 0) stop("⚠️ No VacV samples found!")

# -----------------------------
# Stage patterns
# -----------------------------
stage_patterns <- c(
  Early        = "VACWR[0-9]+_early_aligned_categorized\\.tsv$",
  Intermediate = "VACWR[0-9]+_intermediate_aligned_categorized\\.tsv$",
  Late         = "VACWR[0-9]+_late_aligned_uncategorized\\.tsv$"
)

for (stage in names(stage_patterns)) {
  
  cat("\n🔹 Processing stage:", stage, "\n")
  
  # -----------------------------
  # Read all files matching stage
  # -----------------------------
  stage_df <- map_df(sample_dirs, function(sample_dir) {
    
    # Extract timepoint from folder name
    timepoint <- ifelse(grepl("5h", basename(sample_dir)), "5h", "24h")
    
    files <- list.files(
      sample_dir,
      pattern = stage_patterns[[stage]],
      full.names = TRUE
    )
    
    if (length(files) == 0) return(NULL)
    
    map_df(files, function(f) {
      
      gene <- sub("_aligned_categorized\\.tsv$", "", basename(f))
      
      df <- suppressWarnings(read_table(f, col_names = FALSE))
      
      tibble(
        gene = gene,
        category = df[[ncol(df)]],
        timepoint = timepoint
      ) %>%
        filter(!is.na(category), category %in% names(category_colors))
    })
  }) 
  
  if (nrow(stage_df) == 0) {
    cat("  ⚠️ No data for stage:", stage, "\n")
    next
  }
  #use the loop like above, this is for overall proportion.
  summary_df <- stage_df %>%
    group_by(timepoint, category) %>%
    summarise(total_n = n(), .groups = "drop") %>%
    group_by(timepoint) %>%
    mutate(fraction = total_n / sum(total_n)) %>%
    ungroup()
  
  # -----------------------------
  # Plot Proportion of Category
  # -----------------------------
  p <- ggplot(summary_df, aes(x = timepoint, y = fraction, fill = category)) +
    geom_col(color = "black", width = 0.7) +
    # Add percentage labels
    geom_text(
      aes(label = scales::percent(fraction, accuracy = 0.1)),
      position = position_stack(vjust = 0.5),
      size = 5,
      fontface = "bold",
      color = "white"
    ) +
    scale_fill_manual(values = category_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = paste("Overall", stage, "uAUG Category Distribution"),
      subtitle = "Aggregate of all genes in this stage",
      x = "Timepoint",
      y = "Proportion of Total Reads",
      fill = "uAUG Category"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  # -----------------------------
  # Save
  # -----------------------------
  out_file <- file.path(plot_base, paste0(stage, "_overall_category_proportions.pdf"))
  ggsave(out_file, plot = p, width = 7, height = 6)
  
  cat("  Saved summary plot:", out_file, "\n")
}

#7E
library(tidyverse)
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/uAUG"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/uAUG_combined"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

stage_patterns <- c(
  Early        = "VACWR[0-9]+_early_aligned_categorized\\.tsv$",
  Intermediate = "VACWR[0-9]+_intermediate_aligned_categorized\\.tsv$",
  Late         = "VACWR[0-9]+_late_aligned_categorized\\.tsv$"
)
category_colors <- c(
  "Category1_Productive_Upstream"    = "#0072B2",
  "Category2_Upstream_Stopped"       = "#D55E00",
  "Category3_Overlapping_OutOfFrame" = "#009E73",
  "Category4_No_uAUG"                = "#CC79A7"
)

sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV", sample_dirs)]
# -----------------------------
# 1. COLLECT ALL DATA FROM ALL STAGES
# -----------------------------
all_stages_list <- list()

for (stage_name in names(stage_patterns)) {
  cat("Reading data for:", stage_name, "\n")
  
  for (sample_dir in sample_dirs) {
    timepoint <- ifelse(grepl("5h", basename(sample_dir)), "5h", "24h")
    files <- list.files(sample_dir, pattern = stage_patterns[[stage_name]], full.names = TRUE)
    
    if (length(files) == 0) next
    
    for (f in files) {
      gene_name <- sub("_aligned_categorized\\.tsv$", "", basename(f))
      lines <- readLines(f)
      if (length(lines) == 0) next
      
      categories <- str_extract(lines, "Category[1-4]_[A-Za-z_]+")
      categories <- categories[!is.na(categories)]
      
      if (length(categories) > 0) {
        all_stages_list[[length(all_stages_list) + 1]] <- data.frame(
          gene = as.character(gene_name), 
          category = as.character(categories), 
          timepoint = as.character(timepoint),
          stage = stage_name, # Track which stage this gene belongs to
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# Force into a single master data frame
master_raw_df <- do.call(rbind, all_stages_list)
rownames(master_raw_df) <- NULL

# -----------------------------
# 2. FILTER: 24h Category 4 < 20%
# -----------------------------
genes_to_keep <- master_raw_df %>%
  filter(timepoint == "24h") %>%
  group_by(gene, category) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup() %>%
  filter(category == "Category4_No_uAUG" & fraction < 0.20) %>%
  pull(gene) %>% 
  unique()

# -----------------------------
# 3. PREPARE PLOT DATA
# -----------------------------
grand_plot_data <- master_raw_df %>%
  filter(gene %in% genes_to_keep) %>%
  group_by(stage, timepoint, gene, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(timepoint, gene) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup() %>%
  # Order stages correctly
  mutate(stage = factor(stage, levels = c("Early", "Intermediate", "Late")))

# -----------------------------
# 4. THE GRAND PLOT
# -----------------------------
p_grand <- ggplot(grand_plot_data, aes(x = gene, y = fraction, fill = category)) +
  geom_col(color = "black", linewidth = 0.1) +
  scale_fill_manual(values = category_colors) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  # Facet by Stage (Horizontal) and Timepoint (Vertical)
  facet_grid(timepoint ~ stage, scales = "free_x", space = "free_x") +
  labs(
    title = "Vaccinia Transcriptome: High uAUG Usage Genes",
    subtitle = "Genes filtered for < 20% Category 4 (No uAUG) at 24hpi",
    x = "Gene", y = "Proportion of Reads", fill = "Category"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# -----------------------------
# 5. SAVE
# -----------------------------
# Calculate width based on total unique genes to keep it readable
total_genes <- length(unique(grand_plot_data$gene))
out_grand <- file.path(plot_base, "Summary_24h_Filtered_uAUG.pdf")

ggsave(out_grand, plot = p_grand, width = max(10, total_genes * 0.3), height = 10)

cat("\nGrand Summary Plot saved with", total_genes, "total genes.\n")



#sup8
library(tidyverse)

# -----------------------------
# Directories
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/uAUG"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/uAUG_combined"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# uAUG categories
# -----------------------------
category_colors <- c(
  "Category1_Productive_Upstream"    = "#0072B2",
  "Category2_Upstream_Stopped"       = "#D55E00",
  "Category3_Overlapping_OutOfFrame" = "#009E73",
  "Category4_No_uAUG"                = "#CC79A7"
)

# -----------------------------
# Only VacV samples
# -----------------------------
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV", sample_dirs)]
if (length(sample_dirs) == 0) stop("⚠️ No VacV samples found!")

# -----------------------------
# Stage patterns
# -----------------------------
stage_patterns <- c(
  Early        = "VACWR[0-9]+_early_aligned_categorized\\.tsv$",
  Intermediate = "VACWR[0-9]+_intermediate_aligned_categorized\\.tsv$",
  Late         = "VACWR[0-9]+_late_aligned_categorized\\.tsv$"
)

# -----------------------------
# MAIN LOOP: per stage
# -----------------------------
for (stage in names(stage_patterns)) {
  
  cat("\n🔹 Processing stage:", stage, "\n")
  
  # -----------------------------
  # Read all files matching stage
  # -----------------------------
  stage_df <- map_df(sample_dirs, function(sample_dir) {
    
    # Extract timepoint from folder name
    timepoint <- ifelse(grepl("5h", basename(sample_dir)), "5h", "24h")
    
    files <- list.files(
      sample_dir,
      pattern = stage_patterns[[stage]],
      full.names = TRUE
    )
    
    if (length(files) == 0) return(NULL)
    
    map_df(files, function(f) {
      
      gene <- sub("_aligned_categorized\\.tsv$", "", basename(f))
      
      df <- suppressWarnings(read_table(f, col_names = FALSE))
      
      tibble(
        gene = gene,
        category = df[[ncol(df)]],
        timepoint = timepoint
      ) %>%
        filter(!is.na(category), category %in% names(category_colors))
    })
  }) 
  
  if (nrow(stage_df) == 0) {
    cat("  ⚠️ No data for stage:", stage, "\n")
    next
  }
  
  # -----------------------------
  # Count and fraction per gene x category x timepoint
  # -----------------------------
  stage_df <- stage_df %>%
    count(timepoint, gene, category) %>%
    group_by(timepoint, gene) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup() %>%
    complete(
      timepoint, gene, category = names(category_colors),
      fill = list(n = 0, frac = 0)
    )
  
  # -----------------------------
  # Plot stacked bar
  # -----------------------------
  p <- ggplot(stage_df, aes(x = gene, y = frac, fill = category)) +
    geom_col(color = "black", linewidth = 0.2) +
    scale_fill_manual(values = category_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    facet_wrap(~timepoint, scales = "free_x") +
    labs(
      title = paste(stage, "uAUG distribution by timepoint"),
      x = "Gene",
      y = "Fraction of reads",
      fill = "uAUG category"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # -----------------------------
  # Save
  # -----------------------------
  out_file <- file.path(plot_base, paste0(stage, "_genes_stacked_by_timepoint.pdf"))
  
  ggsave(
    filename = out_file,
    plot = p,
    width = max(6, length(unique(stage_df$gene)) * 0.6),
    height = 4,
    dpi = 300,
    limitsize = FALSE
  )
  
  cat("  Saved:", out_file, "\n")
}

cat("\nAll combined stacked bar plots saved by stage and timepoint\n")



library(tidyverse)

# -----------------------------
# Directories
# -----------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/uAUG"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/uAUG2"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# uAUG categories
# -----------------------------
category_colors <- c(
  "Category1_Productive_Upstream"    = "#0072B2",
  "Category2_Upstream_Stopped"       = "#D55E00",
  "Category3_Overlapping_OutOfFrame" = "#009E73",
  "Category4_No_uAUG"                = "#CC79A7"
)

# -----------------------------
# Only VacV samples
# -----------------------------
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("VacV", sample_dirs)]

if (length(sample_dirs) == 0) stop("⚠️ No VacV samples found!")

# -----------------------------
# Stage patterns
# -----------------------------
stage_patterns <- c(
  Early        = "VACWR[0-9]+_early_unaligned_categorized\\.tsv$",
  Intermediate = "VACWR[0-9]+_intermediate_unaligned_categorized\\.tsv$",
  Late         = "VACWR[0-9]+_late_unaligned_uncategorized\\.tsv$"
)

for (stage in names(stage_patterns)) {
  
  cat("\n🔹 Processing stage:", stage, "\n")
  
  # -----------------------------
  # Read all files matching stage
  # -----------------------------
  stage_df <- map_df(sample_dirs, function(sample_dir) {
    
    # Extract timepoint from folder name
    timepoint <- ifelse(grepl("5h", basename(sample_dir)), "5h", "24h")
    
    files <- list.files(
      sample_dir,
      pattern = stage_patterns[[stage]],
      full.names = TRUE
    )
    
    if (length(files) == 0) return(NULL)
    
    map_df(files, function(f) {
      
      gene <- sub("_unaligned_categorized\\.tsv$", "", basename(f))
      
      df <- suppressWarnings(read_table(f, col_names = FALSE))
      
      tibble(
        gene = gene,
        category = df[[ncol(df)]],
        timepoint = timepoint
      ) %>%
        filter(!is.na(category), category %in% names(category_colors))
    })
  }) 
  
  if (nrow(stage_df) == 0) {
    cat("  ⚠️ No data for stage:", stage, "\n")
    next
  }
  
  
  summary_df <- stage_df %>%
    group_by(timepoint, category) %>%
    summarise(total_n = n(), .groups = "drop") %>%
    group_by(timepoint) %>%
    mutate(fraction = total_n / sum(total_n)) %>%
    ungroup()
  
  # -----------------------------
  # Plot Proportion of Category
  # -----------------------------
  p <- ggplot(summary_df, aes(x = timepoint, y = fraction, fill = category)) +
    geom_col(color = "black", width = 0.7) +
    # Add percentage labels
    geom_text(
      aes(label = scales::percent(fraction, accuracy = 0.1)),
      position = position_stack(vjust = 0.5),
      size = 5,
      fontface = "bold",
      color = "white"
    ) +
    scale_fill_manual(values = category_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = paste("Overall", stage, "uAUG Category Distribution"),
      subtitle = "Aggregate of all genes in this stage",
      x = "Timepoint",
      y = "Proportion of Total Reads",
      fill = "uAUG Category"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  # -----------------------------
  # Save
  # -----------------------------
  out_file <- file.path(plot_base, paste0(stage, "_unaligned_category_proportions.pdf"))
  ggsave(out_file, plot = p, width = 7, height = 6)
  
  cat("  Saved summary plot:", out_file, "\n")
}


