library(tidyverse)
library(patchwork)

# ============================================================
# DIRECTORIES
# ============================================================

indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11/uAUGtest"

plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11/rplot/uAUG_combined_test"

dir.create(
  plot_base,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# uAUG CATEGORIES
# ============================================================

category_colors <- c(
  "Category1_Productive_Upstream"     = "#0072B2",
  "Category2_Upstream_Stopped"        = "#D55E00",
  "Category3_Overlapping_OutOfFrame"  = "#009E73",
  "Category4_No_uAUG"                 = "#CC79A7",
  "Category5_No_CDS_Stop"             = "#999999"
)

category_order <- names(category_colors)


# ============================================================
# SAMPLE DIRECTORIES
#
# INCLUDE:
#   RppH samples
#
# EXCLUDE:
#   5ppp samples
#
# TIMEPOINT:
#   contains "5h" = 5h
#   all other RppH = 24h
# ============================================================

sample_dirs <- list.dirs(
  indir,
  recursive = FALSE,
  full.names = TRUE
)

sample_dirs <- sample_dirs[
  grepl(
    "RppH",
    basename(sample_dirs),
    ignore.case = TRUE
  ) &
    !grepl(
      "5ppp",
      basename(sample_dirs),
      ignore.case = TRUE
    )
]

if (length(sample_dirs) == 0) {
  stop(
    "No RppH sample directories found under:\n",
    indir
  )
}


cat("\n============================================\n")
cat("SAMPLES INCLUDED\n")
cat("============================================\n")

for (x in sample_dirs) {
  
  sample_name <- basename(x)
  
  timepoint <- ifelse(
    grepl(
      "5h",
      sample_name,
      ignore.case = TRUE
    ),
    "5h",
    "24h"
  )
  
  cat(
    sample_name,
    " -> ",
    timepoint,
    "\n",
    sep = ""
  )
}


# ============================================================
# FILE PATTERNS
# ============================================================

stage_patterns_aligned <- c(
  Early =
    "^VACWR[0-9]+_early_utr_highscore_aligned_seqs_categorized\\.tsv$",
  Intermediate =
    "^VACWR[0-9]+_intermediate_utr_highscore_aligned_seqs_categorized\\.tsv$",
  Late =
    "^VACWR[0-9]+_late_utr_highscore_aligned_seqs_categorized\\.tsv$"
)

stage_patterns_unaligned <- c(
  Early =
    "^VACWR[0-9]+_early_utr_unaligned_categorized\\.tsv$",
  Intermediate =
    "^VACWR[0-9]+_intermediate_utr_unaligned_categorized\\.tsv$",
  Late =
    "^VACWR[0-9]+_late_utr_unaligned_categorized\\.tsv$"
)


# ============================================================
# FUNCTION: READ uAUG DATA
# ============================================================

read_uaug_data <- function(stage_patterns, file_type) {
  
  all_list <- list()
  
  for (stage_name in names(stage_patterns)) {
    
    cat(
      "\nReading",
      file_type,
      ":",
      stage_name,
      "\n"
    )
    
    for (sample_dir in sample_dirs) {
      
      sample_name <- basename(sample_dir)
      
      # --------------------------------------------------------
      # TIMEPOINT
      # --------------------------------------------------------
      
      timepoint <- ifelse(
        grepl(
          "5h",
          sample_name,
          ignore.case = TRUE
        ),
        "5h",
        "24h"
      )
      
      # --------------------------------------------------------
      # FIND FILES
      # --------------------------------------------------------
      
      files <- list.files(
        sample_dir,
        pattern = stage_patterns[[stage_name]],
        full.names = TRUE
      )
      
      cat(
        "  ",
        sample_name,
        " [",
        timepoint,
        "] : ",
        length(files),
        " files\n",
        sep = ""
      )
      
      if (length(files) == 0) {
        next
      }
      
      # --------------------------------------------------------
      # READ EACH GENE FILE
      # --------------------------------------------------------
      
      for (f in files) {
        
        filename <- basename(f)
        
        gene_name <- str_extract(
          filename,
          "^VACWR[0-9]+"
        )
        
        if (is.na(gene_name)) {
          warning(
            "Could not extract gene from: ",
            filename
          )
          next
        }
        
        dat <- suppressWarnings(
          read_tsv(
            f,
            col_names = TRUE,
            show_col_types = FALSE
          )
        )
        
        if (
          nrow(dat) == 0 ||
          !"category" %in% colnames(dat)
        ) {
          next
        }
        
        tmp <- tibble(
          gene = gene_name,
          category = as.character(dat$category),
          sample = sample_name,
          timepoint = timepoint,
          stage = stage_name,
          file_type = file_type
        ) %>%
          filter(
            !is.na(category),
            category %in% category_order
          )
        
        if (nrow(tmp) > 0) {
          all_list[[length(all_list) + 1]] <- tmp
        }
      }
    }
  }
  
  if (length(all_list) == 0) {
    stop(
      paste(
        "No usable",
        file_type,
        "data found!"
      )
    )
  }
  
  bind_rows(all_list)
}


# ============================================================
# READ ALIGNED + UNALIGNED
# ============================================================

aligned_df <- read_uaug_data(
  stage_patterns_aligned,
  "Aligned"
)

unaligned_df <- read_uaug_data(
  stage_patterns_unaligned,
  "Unaligned"
)


cat("\n============================================\n")
cat("READ COUNTS\n")
cat("============================================\n")

cat("Aligned:", nrow(aligned_df), "\n")
cat("Unaligned:", nrow(unaligned_df), "\n")


# ============================================================
# CATEGORY COUNTS
# ============================================================

cat("\n============================================\n")
cat("ALIGNED CATEGORY COUNTS\n")
cat("============================================\n")

print(
  aligned_df %>%
    count(
      timepoint,
      category,
      name = "reads"
    ) %>%
    arrange(
      timepoint,
      category
    )
)


cat("\n============================================\n")
cat("UNALIGNED CATEGORY COUNTS\n")
cat("============================================\n")

print(
  unaligned_df %>%
    count(
      timepoint,
      category,
      name = "reads"
    ) %>%
    arrange(
      timepoint,
      category
    )
)


# ============================================================
# FILTER GENES
#
# Remove gene/timepoint combinations where:
# Category4_No_uAUG >= 20%
#
# Category 1-5 are all included in denominator.
# ============================================================

get_genes_to_keep <- function(df) {
  
  category_fraction <- df %>%
    
    count(
      timepoint,
      gene,
      category,
      name = "n"
    ) %>%
    
    group_by(
      timepoint,
      gene
    ) %>%
    
    mutate(
      fraction = n / sum(n)
    ) %>%
    
    ungroup()
  
  
  cat4_remove <- category_fraction %>%
    
    filter(
      category == "Category4_No_uAUG",
      fraction >= 0.20
    ) %>%
    
    distinct(
      timepoint,
      gene
    )
  
  
  df %>%
    
    distinct(
      timepoint,
      gene
    ) %>%
    
    anti_join(
      cat4_remove,
      by = c(
        "timepoint",
        "gene"
      )
    )
}


aligned_genes_to_keep <- get_genes_to_keep(
  aligned_df
)

unaligned_genes_to_keep <- get_genes_to_keep(
  unaligned_df
)


# ============================================================
# PREPARE PLOT DATA
# ============================================================

prepare_plot_data <- function(
    df,
    genes_to_keep
) {
  
  counts <- df %>%
    
    semi_join(
      genes_to_keep,
      by = c(
        "timepoint",
        "gene"
      )
    ) %>%
    
    count(
      timepoint,
      stage,
      gene,
      category,
      name = "n"
    )
  
  
  gene_info <- counts %>%
    
    distinct(
      timepoint,
      stage,
      gene
    )
  
  
  complete_grid <- gene_info %>%
    
    crossing(
      category = category_order
    )
  
  
  plot_data <- complete_grid %>%
    
    left_join(
      counts,
      by = c(
        "timepoint",
        "stage",
        "gene",
        "category"
      )
    ) %>%
    
    mutate(
      n = replace_na(
        n,
        0L
      )
    ) %>%
    
    group_by(
      timepoint,
      stage,
      gene
    ) %>%
    
    mutate(
      total_reads = sum(n),
      fraction = n / total_reads
    ) %>%
    
    ungroup() %>%
    
    mutate(
      stage = factor(
        stage,
        levels = c(
          "Early",
          "Intermediate",
          "Late"
        )
      ),
      
      timepoint = factor(
        timepoint,
        levels = c(
          "5h",
          "24h"
        )
      ),
      
      category = factor(
        category,
        levels = category_order
      )
    )
  
  plot_data
}


aligned_plot_data <- prepare_plot_data(
  aligned_df,
  aligned_genes_to_keep
)

unaligned_plot_data <- prepare_plot_data(
  unaligned_df,
  unaligned_genes_to_keep
)


# ============================================================
# GENE ORDER
# ============================================================

get_gene_order <- function(dat) {
  
  dat %>%
    
    distinct(gene) %>%
    
    mutate(
      gene_num = as.numeric(
        str_extract(
          gene,
          "[0-9]+"
        )
      )
    ) %>%
    
    arrange(
      gene_num
    ) %>%
    
    pull(gene)
}


# ============================================================
# FUNCTION: GENE-LEVEL TIMEPOINT PLOT
# ============================================================

make_timepoint_plot <- function(
    plot_data,
    tp,
    data_label
) {
  
  dat <- plot_data %>%
    filter(
      timepoint == tp
    )
  
  if (nrow(dat) == 0) {
    
    return(
      ggplot() +
        theme_void() +
        labs(
          title = paste(
            data_label,
            tp,
            "- No data"
          )
        )
    )
  }
  
  
  gene_order <- get_gene_order(dat)
  
  dat <- dat %>%
    
    mutate(
      gene = factor(
        gene,
        levels = gene_order
      )
    )
  
  
  ggplot(
    dat,
    aes(
      x = gene,
      y = fraction,
      fill = category
    )
  ) +
    
    geom_col(
      color = "black",
      linewidth = 0.15,
      width = 0.9
    ) +
    
    scale_fill_manual(
      values = category_colors,
      breaks = category_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      limits = c(
        0,
        1
      ),
      expand = expansion(
        mult = c(
          0,
          0.01
        )
      )
    ) +
    
    facet_grid(
      . ~ stage,
      scales = "free_x",
      space = "free_x"
    ) +
    
    labs(
      title = paste0(
        "Vaccinia Transcriptome: ",
        data_label,
        " ",
        tp
      ),
      subtitle =
        "Gene/timepoint combinations with >=20% Category 4 excluded",
      x = "Gene",
      y = "Proportion of Reads",
      fill = "Category"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 10
      ),
      
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      
      strip.background = element_rect(
        fill = "grey90"
      ),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      legend.position = "bottom",
      
      panel.spacing.x = unit(
        0.15,
        "lines"
      ),
      
      plot.margin = margin(
        5,
        5,
        5,
        5
      )
    )
}


# ============================================================
# GENE-LEVEL PLOTS
# ============================================================

aligned_5h <- make_timepoint_plot(
  aligned_plot_data,
  "5h",
  "Aligned"
)

aligned_24h <- make_timepoint_plot(
  aligned_plot_data,
  "24h",
  "Aligned"
)

unaligned_5h <- make_timepoint_plot(
  unaligned_plot_data,
  "5h",
  "Unaligned"
)

unaligned_24h <- make_timepoint_plot(
  unaligned_plot_data,
  "24h",
  "Unaligned"
)


# ============================================================
# COMBINE GENE-LEVEL PLOTS
# ============================================================

combined_gene_plot <-
  
  (
    aligned_5h /
      aligned_24h /
      unaligned_5h /
      unaligned_24h
  ) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1,
      1,
      1
    ),
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )


# ============================================================
# SAVE GENE-LEVEL PDF
# ============================================================

out_file_gene <- file.path(
  plot_base,
  "Summary_uAUG_Aligned_Unaligned_Filtered.pdf"
)


max_genes <- max(
  n_distinct(
    aligned_plot_data$gene
  ),
  n_distinct(
    unaligned_plot_data$gene
  )
)


ggsave(
  filename = out_file_gene,
  plot = combined_gene_plot,
  width = max(
    10,
    max_genes * 0.30
  ),
  height = 15,
  limitsize = FALSE
)


cat(
  "\nSaved gene-level plot:\n",
  out_file_gene,
  "\n"
)


# ============================================================
# STAGE × TIMEPOINT SUMMARY
# ============================================================

make_stage_timepoint_plot <- function(
    plot_data,
    data_label
) {
  
  # ----------------------------------------------------------
  # Sum reads over retained genes
  # ----------------------------------------------------------
  
  dat <- plot_data %>%
    
    group_by(
      stage,
      timepoint,
      category
    ) %>%
    
    summarise(
      n = sum(n),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # Existing stage/timepoint combinations
  # ----------------------------------------------------------
  
  stage_tp <- dat %>%
    
    distinct(
      stage,
      timepoint
    )
  
  
  # ----------------------------------------------------------
  # Add zeroes for missing categories
  # ----------------------------------------------------------
  
  dat <- stage_tp %>%
    
    crossing(
      category = category_order
    ) %>%
    
    left_join(
      dat,
      by = c(
        "stage",
        "timepoint",
        "category"
      )
    ) %>%
    
    mutate(
      n = replace_na(
        n,
        0
      )
    ) %>%
    
    group_by(
      stage,
      timepoint
    ) %>%
    
    # IMPORTANT:
    # Do not use if_else() here.
    # n is a vector while sum(n) is one value.
    mutate(
      total_n = sum(n),
      fraction = n / total_n
    ) %>%
    
    ungroup() %>%
    
    mutate(
      stage = factor(
        stage,
        levels = c(
          "Early",
          "Intermediate",
          "Late"
        )
      ),
      
      timepoint = factor(
        timepoint,
        levels = c(
          "5h",
          "24h"
        )
      ),
      
      category = factor(
        category,
        levels = category_order
      )
    )
  
  
  # ----------------------------------------------------------
  # PLOT
  # ----------------------------------------------------------
  
  ggplot(
    dat,
    aes(
      x = timepoint,
      y = fraction,
      fill = category
    )
  ) +
    
    geom_col(
      color = "black",
      linewidth = 0.2,
      width = 0.7
    ) +
    
    geom_text(
      aes(
        label = case_when(
          fraction >= 0.03 ~
            scales::percent(
              fraction,
              accuracy = 1
            ),
          TRUE ~ ""
        )
      ),
      position = position_stack(
        vjust = 0.5
      ),
      size = 3
    ) +
    
    scale_fill_manual(
      values = category_colors,
      breaks = category_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(
        0,
        1
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    facet_grid(
      . ~ stage
    ) +
    
    labs(
      title = paste0(
        "Vaccinia Transcriptome: ",
        data_label
      ),
      subtitle =
        "All retained genes combined within each stage",
      x = NULL,
      y = "Proportion of Reads",
      fill = "Category"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      
      axis.text.x = element_text(
        size = 11
      ),
      
      strip.background = element_rect(
        fill = "grey90"
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      
      legend.position = "bottom",
      
      panel.spacing.x = unit(
        0.3,
        "lines"
      ),
      
      plot.margin = margin(
        5,
        5,
        5,
        5
      )
    )
}


# ============================================================
# MAKE STAGE SUMMARY PLOTS
# ============================================================

aligned_stage_plot <- make_stage_timepoint_plot(
  aligned_plot_data,
  "Aligned"
)

unaligned_stage_plot <- make_stage_timepoint_plot(
  unaligned_plot_data,
  "Unaligned"
)


# ============================================================
# COMBINE STAGE SUMMARY PLOTS
# ============================================================

combined_stage_plot <-
  
  (
    aligned_stage_plot /
      unaligned_stage_plot
  ) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1
    ),
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )


# ============================================================
# SAVE STAGE SUMMARY PDF
# ============================================================

out_file_stage <- file.path(
  plot_base,
  "Summary_uAUG_Stage_5h_24h.pdf"
)


ggsave(
  filename = out_file_stage,
  plot = combined_stage_plot,
  width = 7,
  height = 8,
  limitsize = FALSE
)


cat(
  "\nSaved stage summary:\n",
  out_file_stage,
  "\n"
)


# ============================================================
# SAVE ALIGNED PLOT DATA
# ============================================================

write_tsv(
  aligned_plot_data %>%
    
    mutate(
      stage = as.character(stage),
      timepoint = as.character(timepoint),
      category = as.character(category),
      gene = as.character(gene)
    ),
  
  file.path(
    plot_base,
    "Aligned_uAUG_plot_data.tsv"
  )
)


# ============================================================
# SAVE UNALIGNED PLOT DATA
# ============================================================

write_tsv(
  unaligned_plot_data %>%
    
    mutate(
      stage = as.character(stage),
      timepoint = as.character(timepoint),
      category = as.character(category),
      gene = as.character(gene)
    ),
  
  file.path(
    plot_base,
    "Unaligned_uAUG_plot_data.tsv"
  )
)


# ============================================================
# DONE
# ============================================================

cat("\n============================================\n")
cat("uAUG ANALYSIS COMPLETE\n")
cat("============================================\n")

cat("Included: RppH samples only\n")
cat("Excluded: 5ppp samples\n")
cat("5h: sample name contains '5h'\n")
cat("24h: all remaining RppH samples\n")
cat("Categories included: Category 1-5\n")

cat(
  "\nOutput directory:\n",
  plot_base,
  "\n"
)

cat("============================================\n")








library(tidyverse)
library(patchwork)

# ============================================================
# DIRECTORIES
# ============================================================

indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11/uAUGtest"

plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11/rplot/uAUG_combined_test"

dir.create(
  plot_base,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# uAUG CATEGORIES
# ============================================================

category_colors <- c(
  "Category1_Productive_Upstream"    = "#0072B2",
  "Category2_Upstream_Stopped"       = "#D55E00",
  "Category3_Overlapping_OutOfFrame" = "#009E73",
  "Category4_No_uAUG"                = "#CC79A7",
  "Category5_No_CDS_Stop"            = "#999999"
)

category_order <- names(category_colors)


# ============================================================
# SAMPLE DIRECTORIES
#
# INCLUDE:
#   RppH samples
#
# EXCLUDE:
#   5ppp samples
#
# TIMEPOINT:
#   contains "5h" = 5h
#   all other RppH = 24h
# ============================================================

sample_dirs <- list.dirs(
  indir,
  recursive = FALSE,
  full.names = TRUE
)

sample_dirs <- sample_dirs[
  grepl(
    "RppH",
    basename(sample_dirs),
    ignore.case = TRUE
  ) &
    !grepl(
      "5ppp",
      basename(sample_dirs),
      ignore.case = TRUE
    )
]

if (length(sample_dirs) == 0) {
  stop(
    "No RppH sample directories found under:\n",
    indir
  )
}


cat("\n============================================\n")
cat("SAMPLES INCLUDED\n")
cat("============================================\n")

for (x in sample_dirs) {
  
  sample_name <- basename(x)
  
  timepoint <- ifelse(
    grepl(
      "5h",
      sample_name,
      ignore.case = TRUE
    ),
    "5h",
    "24h"
  )
  
  cat(
    sample_name,
    " -> ",
    timepoint,
    "\n",
    sep = ""
  )
}


# ============================================================
# FILE PATTERNS
# ============================================================

stage_patterns_aligned <- c(
  Early =
    "^VACWR[0-9]+_early_utr_highscore_aligned_seqs_categorized\\.tsv$",
  Intermediate =
    "^VACWR[0-9]+_intermediate_utr_highscore_aligned_seqs_categorized\\.tsv$",
  Late =
    "^VACWR[0-9]+_late_utr_highscore_aligned_seqs_categorized\\.tsv$"
)

stage_patterns_unaligned <- c(
  Early =
    "^VACWR[0-9]+_early_utr_unaligned_categorized\\.tsv$",
  Intermediate =
    "^VACWR[0-9]+_intermediate_utr_unaligned_categorized\\.tsv$",
  Late =
    "^VACWR[0-9]+_late_utr_unaligned_categorized\\.tsv$"
)


# ============================================================
# FUNCTION: READ uAUG DATA
# ============================================================

read_uaug_data <- function(stage_patterns, file_type) {
  
  all_list <- list()
  
  for (stage_name in names(stage_patterns)) {
    
    cat(
      "\nReading",
      file_type,
      ":",
      stage_name,
      "\n"
    )
    
    for (sample_dir in sample_dirs) {
      
      sample_name <- basename(sample_dir)
      
      # --------------------------------------------------------
      # TIMEPOINT
      # --------------------------------------------------------
      
      timepoint <- ifelse(
        grepl(
          "5h",
          sample_name,
          ignore.case = TRUE
        ),
        "5h",
        "24h"
      )
      
      # --------------------------------------------------------
      # FIND FILES
      # --------------------------------------------------------
      
      files <- list.files(
        sample_dir,
        pattern = stage_patterns[[stage_name]],
        full.names = TRUE
      )
      
      cat(
        "  ",
        sample_name,
        " [",
        timepoint,
        "] : ",
        length(files),
        " files\n",
        sep = ""
      )
      
      if (length(files) == 0) {
        next
      }
      
      # --------------------------------------------------------
      # READ EACH GENE FILE
      # --------------------------------------------------------
      
      for (f in files) {
        
        filename <- basename(f)
        
        gene_name <- str_extract(
          filename,
          "^VACWR[0-9]+"
        )
        
        if (is.na(gene_name)) {
          warning(
            "Could not extract gene from: ",
            filename
          )
          next
        }
        
        dat <- suppressWarnings(
          read_tsv(
            f,
            col_names = TRUE,
            show_col_types = FALSE
          )
        )
        
        if (
          nrow(dat) == 0 ||
          !"category" %in% colnames(dat)
        ) {
          next
        }
        
        tmp <- tibble(
          gene = gene_name,
          category = as.character(dat$category),
          sample = sample_name,
          timepoint = timepoint,
          stage = stage_name,
          file_type = file_type
        ) %>%
          filter(
            !is.na(category),
            category %in% category_order
          )
        
        if (nrow(tmp) > 0) {
          all_list[[length(all_list) + 1]] <- tmp
        }
      }
    }
  }
  
  if (length(all_list) == 0) {
    stop(
      paste(
        "No usable",
        file_type,
        "data found!"
      )
    )
  }
  
  bind_rows(all_list)
}


# ============================================================
# READ ALIGNED + UNALIGNED
# ============================================================

aligned_df <- read_uaug_data(
  stage_patterns_aligned,
  "Aligned"
)

unaligned_df <- read_uaug_data(
  stage_patterns_unaligned,
  "Unaligned"
)


cat("\n============================================\n")
cat("READ COUNTS\n")
cat("============================================\n")

cat("Aligned:", nrow(aligned_df), "\n")
cat("Unaligned:", nrow(unaligned_df), "\n")


# ============================================================
# CATEGORY COUNTS
# ============================================================

cat("\n============================================\n")
cat("ALIGNED CATEGORY COUNTS\n")
cat("============================================\n")

print(
  aligned_df %>%
    count(
      timepoint,
      category,
      name = "reads"
    ) %>%
    arrange(
      timepoint,
      category
    )
)


cat("\n============================================\n")
cat("UNALIGNED CATEGORY COUNTS\n")
cat("============================================\n")

print(
  unaligned_df %>%
    count(
      timepoint,
      category,
      name = "reads"
    ) %>%
    arrange(
      timepoint,
      category
    )
)


# ============================================================
# PREPARE PLOT DATA
#
# NO GENE FILTERING
# All genes are retained.
# ============================================================

prepare_plot_data <- function(df) {
  
  counts <- df %>%
    count(
      timepoint,
      stage,
      gene,
      category,
      name = "n"
    )
  
  gene_info <- counts %>%
    distinct(
      timepoint,
      stage,
      gene
    )
  
  complete_grid <- gene_info %>%
    crossing(
      category = category_order
    )
  
  plot_data <- complete_grid %>%
    left_join(
      counts,
      by = c(
        "timepoint",
        "stage",
        "gene",
        "category"
      )
    ) %>%
    mutate(
      n = replace_na(n, 0L)
    ) %>%
    group_by(
      timepoint,
      stage,
      gene
    ) %>%
    mutate(
      total_reads = sum(n),
      fraction = n / total_reads
    ) %>%
    ungroup() %>%
    mutate(
      stage = factor(
        stage,
        levels = c(
          "Early",
          "Intermediate",
          "Late"
        )
      ),
      timepoint = factor(
        timepoint,
        levels = c(
          "5h",
          "24h"
        )
      ),
      category = factor(
        category,
        levels = category_order
      )
    )
  
  plot_data
}


aligned_plot_data <- prepare_plot_data(
  aligned_df
)

unaligned_plot_data <- prepare_plot_data(
  unaligned_df
)


# ============================================================
# GENE ORDER
# ============================================================

get_gene_order <- function(dat) {
  
  dat %>%
    distinct(gene) %>%
    mutate(
      gene_num = as.numeric(
        str_extract(
          gene,
          "[0-9]+"
        )
      )
    ) %>%
    arrange(
      gene_num
    ) %>%
    pull(gene)
}


# ============================================================
# FUNCTION: GENE-LEVEL TIMEPOINT PLOT
# ============================================================

make_timepoint_plot <- function(
    plot_data,
    tp,
    data_label
) {
  
  dat <- plot_data %>%
    filter(
      timepoint == tp
    )
  
  if (nrow(dat) == 0) {
    
    return(
      ggplot() +
        theme_void() +
        labs(
          title = paste(
            data_label,
            tp,
            "- No data"
          )
        )
    )
  }
  
  gene_order <- get_gene_order(dat)
  
  dat <- dat %>%
    mutate(
      gene = factor(
        gene,
        levels = gene_order
      )
    )
  
  ggplot(
    dat,
    aes(
      x = gene,
      y = fraction,
      fill = category
    )
  ) +
    
    geom_col(
      color = "black",
      linewidth = 0.15,
      width = 0.9
    ) +
    
    scale_fill_manual(
      values = category_colors,
      breaks = category_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      limits = c(
        0,
        1
      ),
      expand = expansion(
        mult = c(
          0,
          0.01
        )
      )
    ) +
    
    facet_grid(
      . ~ stage,
      scales = "free_x",
      space = "free_x"
    ) +
    
    labs(
      title = paste0(
        "Vaccinia Transcriptome: ",
        data_label,
        " ",
        tp
      ),
      subtitle = "All genes included",
      x = "Gene",
      y = "Proportion of Reads",
      fill = "Category"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 10
      ),
      
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      
      strip.background = element_rect(
        fill = "grey90"
      ),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      legend.position = "bottom",
      
      panel.spacing.x = unit(
        0.15,
        "lines"
      ),
      
      plot.margin = margin(
        5,
        5,
        5,
        5
      )
    )
}


# ============================================================
# GENE-LEVEL PLOTS
# ============================================================

aligned_5h <- make_timepoint_plot(
  aligned_plot_data,
  "5h",
  "Aligned"
)

aligned_24h <- make_timepoint_plot(
  aligned_plot_data,
  "24h",
  "Aligned"
)

unaligned_5h <- make_timepoint_plot(
  unaligned_plot_data,
  "5h",
  "Unaligned"
)

unaligned_24h <- make_timepoint_plot(
  unaligned_plot_data,
  "24h",
  "Unaligned"
)


# ============================================================
# COMBINE GENE-LEVEL PLOTS
# ============================================================

combined_gene_plot <-
  (
    aligned_5h /
      aligned_24h /
      unaligned_5h /
      unaligned_24h
  ) +
  patchwork::plot_layout(
    heights = c(
      1,
      1,
      1,
      1
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )


# ============================================================
# SAVE GENE-LEVEL PDF
# ============================================================

out_file_gene <- file.path(
  plot_base,
  "Summary_uAUG_Aligned_Unaligned_AllGenes.pdf"
)

max_genes <- max(
  n_distinct(
    aligned_plot_data$gene
  ),
  n_distinct(
    unaligned_plot_data$gene
  )
)

ggsave(
  filename = out_file_gene,
  plot = combined_gene_plot,
  width = max(
    10,
    max_genes * 0.30
  ),
  height = 15,
  limitsize = FALSE
)

cat(
  "\nSaved gene-level plot:\n",
  out_file_gene,
  "\n"
)


# ============================================================
# STAGE × TIMEPOINT SUMMARY
# ============================================================

make_stage_timepoint_plot <- function(
    plot_data,
    data_label
) {
  
  dat <- plot_data %>%
    group_by(
      stage,
      timepoint,
      category
    ) %>%
    summarise(
      n = sum(n),
      .groups = "drop"
    )
  
  stage_tp <- dat %>%
    distinct(
      stage,
      timepoint
    )
  
  dat <- stage_tp %>%
    crossing(
      category = category_order
    ) %>%
    left_join(
      dat,
      by = c(
        "stage",
        "timepoint",
        "category"
      )
    ) %>%
    mutate(
      n = replace_na(
        n,
        0
      )
    ) %>%
    group_by(
      stage,
      timepoint
    ) %>%
    mutate(
      total_n = sum(n),
      fraction = n / total_n
    ) %>%
    ungroup() %>%
    mutate(
      stage = factor(
        stage,
        levels = c(
          "Early",
          "Intermediate",
          "Late"
        )
      ),
      timepoint = factor(
        timepoint,
        levels = c(
          "5h",
          "24h"
        )
      ),
      category = factor(
        category,
        levels = category_order
      )
    )
  
  ggplot(
    dat,
    aes(
      x = timepoint,
      y = fraction,
      fill = category
    )
  ) +
    
    geom_col(
      color = "black",
      linewidth = 0.2,
      width = 0.7
    ) +
    
    geom_text(
      aes(
        label = case_when(
          fraction >= 0.03 ~
            scales::percent(
              fraction,
              accuracy = 1
            ),
          TRUE ~ ""
        )
      ),
      position = position_stack(
        vjust = 0.5
      ),
      size = 3
    ) +
    
    scale_fill_manual(
      values = category_colors,
      breaks = category_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(
        0,
        1
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    facet_grid(
      . ~ stage
    ) +
    
    labs(
      title = paste0(
        "Vaccinia Transcriptome: ",
        data_label
      ),
      subtitle = "All genes combined within each stage",
      x = NULL,
      y = "Proportion of Reads",
      fill = "Category"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      
      axis.text.x = element_text(
        size = 11
      ),
      
      strip.background = element_rect(
        fill = "grey90"
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      
      legend.position = "bottom",
      
      panel.spacing.x = unit(
        0.3,
        "lines"
      ),
      
      plot.margin = margin(
        5,
        5,
        5,
        5
      )
    )
}


# ============================================================
# MAKE STAGE SUMMARY PLOTS
# ============================================================

aligned_stage_plot <- make_stage_timepoint_plot(
  aligned_plot_data,
  "Aligned"
)

unaligned_stage_plot <- make_stage_timepoint_plot(
  unaligned_plot_data,
  "Unaligned"
)


# ============================================================
# COMBINE STAGE SUMMARY
# ============================================================

combined_stage_plot <-
  (
    aligned_stage_plot /
      unaligned_stage_plot
  ) +
  patchwork::plot_layout(
    heights = c(
      1,
      1
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )


# ============================================================
# SAVE STAGE SUMMARY PDF
# ============================================================

out_file_stage <- file.path(
  plot_base,
  "Summary_uAUG_Stage_5h_24h_AllGenes.pdf"
)

ggsave(
  filename = out_file_stage,
  plot = combined_stage_plot,
  width = 7,
  height = 8,
  limitsize = FALSE
)

cat(
  "\nSaved stage summary:\n",
  out_file_stage,
  "\n"
)


# ============================================================
# SAVE ALIGNED PLOT DATA
# ============================================================

write_tsv(
  aligned_plot_data %>%
    mutate(
      stage = as.character(stage),
      timepoint = as.character(timepoint),
      category = as.character(category),
      gene = as.character(gene)
    ),
  file.path(
    plot_base,
    "Aligned_uAUG_plot_data_AllGenes.tsv"
  )
)


# ============================================================
# SAVE UNALIGNED PLOT DATA
# ============================================================

write_tsv(
  unaligned_plot_data %>%
    mutate(
      stage = as.character(stage),
      timepoint = as.character(timepoint),
      category = as.character(category),
      gene = as.character(gene)
    ),
  file.path(
    plot_base,
    "Unaligned_uAUG_plot_data_AllGenes.tsv"
  )
)


# ============================================================
# DONE
# ============================================================

cat("\n============================================\n")
cat("uAUG ANALYSIS COMPLETE\n")
cat("============================================\n")

cat("Included: RppH samples only\n")
cat("Excluded: 5ppp samples\n")
cat("5h: sample name contains '5h'\n")
cat("24h: all remaining RppH samples\n")
cat("Categories included: Category 1-5\n")
cat("Gene filtering: NONE\n")

cat(
  "\nOutput directory:\n",
  plot_base,
  "\n"
)

cat("============================================\n")
