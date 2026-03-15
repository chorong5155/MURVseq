library(readr)
library(dplyr)
library(ggplot2)

root_dir <- "/Users/chorongpark/Desktop/MUSVseq/micropolyA_processed"
outdir <- "/Users/chorongpark/Desktop/MUSVseq/micropolyA_processed/restseq_plots_selected"
if (!dir.exists(outdir)) dir.create(outdir)

# Genes you want to test
genes_to_test <- c(
  "VACWR028","VACWR025", "VACWR122", "VACWR103"
)

# Get all first-level sample folders containing "VacV_RppH_"
sample_dirs <- list.dirs(root_dir, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("VacV_RppH_", sample_dirs, ignore.case = TRUE)]

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  cat("Processing sample:", sample_name, "\n")
  
  # Create output folder for this sample
  sample_outdir <- file.path(outdir, sample_name)
  if (!dir.exists(sample_outdir)) dir.create(sample_outdir)
  
  # Get gene subfolders and filter to only genes_to_test
  gene_dirs <- list.dirs(sample_dir, recursive = FALSE)
  gene_dirs <- gene_dirs[sapply(gene_dirs, function(g) any(sapply(genes_to_test, function(x) grepl(x, g, ignore.case = TRUE))))]
  
  for (gene_dir in gene_dirs) {
    gene_name <- basename(gene_dir)
    
    files <- list.files(gene_dir, pattern = "subtracted_17nt.tsv$", full.names = TRUE)
    
    for (file in files) {
      cat("  Processing file:", basename(file), "\n")
      
      df <- read_tsv(file, show_col_types = FALSE) %>%
        mutate(polyA_length = as.numeric(as.character(polyA_length)))
      
      if(nrow(df) == 0 || any(is.na(df$polyA_length)) || !("rest_sequence" %in% colnames(df))) {
        cat("    Skipping", basename(file), "- no data or missing polyA_length/rest_sequence\n")
        next
      }
      
      # Summarize counts
      df_summary <- df %>%
        group_by(polyA_length, rest_sequence) %>%
        summarise(Count = n(), .groups = "drop")
      
      # Select top 30 most frequent rest sequences
      top30_rest <- df_summary %>%
        group_by(rest_sequence) %>%
        summarise(Total = sum(Count), .groups = "drop") %>%
        arrange(desc(Total)) %>%
        slice_head(n = 30)
      
      df_top30 <- df_summary %>%
        filter(rest_sequence %in% top30_rest$rest_sequence)
      
      # Plot
      p <- ggplot(df_top30, aes(x = reorder(rest_sequence, Count), y = polyA_length)) +
        geom_tile(aes(fill = Count), color = "white") +
        scale_fill_gradient(low = "lightblue", high = "darkblue") +
        labs(
          title = paste0("Top 30 Rest Sequences vs PolyA Length: ", gene_name),
          x = "Rest Sequence", y = "PolyA Length"
        ) +
        coord_flip() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
              axis.text.y = element_text(size = 5),
              axis.title = element_text(size = 6),
              plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
              plot.margin = margin(10, 10, 10, 10))
      
      outfile <- file.path(sample_outdir, paste0(tools::file_path_sans_ext(basename(file)), ".png"))
      ggsave(filename = outfile, plot = p, width = 4, height = 6)
      cat("    Saved plot to", outfile, "\n")
    }
  }
}

cat("All done!\n")

library(tidyverse)
library(Hmisc)

# -------------------------------
# Base setup (keep as is)
# -------------------------------
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/micropolyA_processed"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/micropolyA"
dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

phases <- c("early", "intermediate", "late")
phase_colors <- c(
  early = "skyblue",
  intermediate = "orange",
  late = "salmon"
)

get_timepoint <- function(name) {
  if (grepl("VacV_RppH_5h", name)) return("5h")
  if (grepl("VacV_RppH", name)) return("24h")
  return(NA)
}

# -------------------------------
# Target gene
# -------------------------------
target_gene <- <- c(
  "VACWR025_late","VACWR028_early", "VACWR122_intermediate", "VACWR103_early"
)

# -------------------------------
# Get top-level sample folders
# -------------------------------
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
# Keep only RppH samples
sample_dirs <- sample_dirs[grepl("VacV_RppH", sample_dirs)]
samples <- basename(sample_dirs)

# -------------------------------
# Loop through each sample
# -------------------------------
for (sample_dir in sample_dirs) {
  
  sample_name <- basename(sample_dir)
  sample_plot_dir <- file.path(plot_base, sample_name)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get gene folders inside the sample
  gene_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
  
  # Only process the target gene
  gene_dirs <- gene_dirs[basename(gene_dirs) == target_gene]
  if (length(gene_dirs) == 0) next  # skip if target_gene not in this sample
  
  # Read all polyA data for the target gene
  sample_polyA <- map_dfr(gene_dirs, function(gene_dir) {
    tsv_file <- list.files(gene_dir, pattern = "_full.tsv$", full.names = TRUE)
    if (length(tsv_file) == 0) return(NULL)
    
    read_tsv(tsv_file, col_types = cols(
      header = col_character(),
      polyA_seq = col_character(),
      polyA_length = col_double(),
      rest_sequence = col_character()
    ), show_col_types = FALSE) %>%
      mutate(
        polyA_length = as.numeric(polyA_length),
        gene = basename(gene_dir),
        Timepoint = get_timepoint(sample_name),
        phase_color = case_when(
          grepl("early", gene, ignore.case = TRUE) ~ "early",
          grepl("intermediate", gene, ignore.case = TRUE) ~ "intermediate",
          grepl("late", gene, ignore.case = TRUE) ~ "late",
          TRUE ~ "other"
        )
      ) %>%
      filter(!is.na(polyA_length))
  })
  
  if (nrow(sample_polyA) == 0) next  # skip if no data
  
  # -------------------------------
  # Plot histogram per sample
  # -------------------------------
  p <- ggplot(sample_polyA, aes(x = polyA_length, y = ..count.., fill = phase_color)) +
    geom_histogram(binwidth = 1, color = "black", size = 0.3) +
    scale_fill_manual(values = phase_colors) +
    scale_x_continuous(breaks = c(3, 6, 9, 15, 30, 40)) +
    scale_y_log10() +
    labs(
      title = paste0(sample_name, " Poly(A) Tail Length Distribution: ", target_gene),
      x = "Poly(A) Tail Length (nt)",
      y = "Read Count",
      fill = "Phase"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      strip.text = element_text(size = 6, face = "bold"),
      axis.text = element_text(size = 6)
    )
  
  # -------------------------------
  # Save PDF
  # -------------------------------
  out_file <- file.path(sample_plot_dir, paste0(target_gene, "_polyA_histogram.pdf"))
  ggsave(out_file, plot = p, width = 8, height = 4, dpi = 300)
  
  message("Saved plot for sample ", sample_name)
}

