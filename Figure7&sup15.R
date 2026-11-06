library(tidyverse)

# ============================================================
# SETTINGS
# ============================================================

indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_3/combined_TSS"

plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_3/rplot/TSS_combined"

dir.create(
  plot_base,
  recursive = TRUE,
  showWarnings = FALSE)

# EXACT SAMPLE NAMES
sample_5h <- "HAP1_VacV_RppH_5h_19215FL-12-03_S3_L003"
sample_24h <- "HAP1_VacV_RppH_4_19215FL-12-04_S4_L003"


MIN_READS <- 100
N_ITERATIONS <- 100
SWITCH_THRESHOLD <- 0.50

set.seed(12345)



TSS_TYPES <- c(
  "natural_TSS",
  "classicslippage_TSS",
  "slippage_TSS")


load_sample <- function(sample_name, timepoint) {
  
  sample_dir <- file.path(
    indir,
    sample_name  )
  
  if (!dir.exists(sample_dir)) {
    
    stop(
      "\nSample directory not found:\n",
      sample_dir,
      "\n"
    ) }
  
  cat("\n====================================================\n")
  
  cat(
    "Loading:",
    sample_name,
    "\n")
  
  cat(
    "Timepoint:",
    timepoint,
    "\n")
  
  cat( "====================================================\n")
  
  
  gene_dirs <- list.dirs(
    sample_dir,
    recursive = FALSE,
    full.names = TRUE
  )
  
  all_data <- list()
  
  
  for (gene_path in gene_dirs) {
    
    gene_name <- basename(gene_path)
    
    bed_files <- list.files(
      gene_path,
      pattern = "\\.bed$",
      full.names = TRUE  )
    
    
    for (bed in bed_files) {
      
      filename <- basename(bed)
      
      
    
      TSS <- case_when(
        
        grepl(
          "natural_TSS_counts",
          filename,
          ignore.case = TRUE
        ) ~ "natural_TSS",
        
        grepl(
          "classicslippage_TSS_counts",
          filename,
          ignore.case = TRUE
        ) ~ "classicslippage_TSS",
        
        grepl(
          "micro_TSS_counts",
          filename,
          ignore.case = TRUE
        ) ~ "slippage_TSS",
        
        TRUE ~ NA_character_
      )
      
      
      if (is.na(TSS)) {
        next }
      
      
      # ------------------------------------------------------
      # Read BED
      # ------------------------------------------------------
      
      dat <- tryCatch(
        
        read.table(
          bed,
          sep = "\t",
          header = FALSE,
          stringsAsFactors = FALSE
        )[, 1:4],
        
        error = function(e) NULL
      )
      
      
      if (is.null(dat)) {
        next
      }
      
      if (nrow(dat) == 0) {
        next
      }
      
      
      colnames(dat) <- c(
        "chr",
        "start",
        "end",
        "count"
      )
      
      
  
      dat <- dat %>%
        filter(
          grepl(
            "VacV",
            chr,
            ignore.case = TRUE  )
        )
      
      
      if (nrow(dat) == 0) {
        next
      }
      
      
     
      
      dat <- dat %>%
        mutate(
          
          count = as.numeric(count),
          
          pos = as.numeric(start),
          
          gene = gene_name,
          
          type = TSS,
          
          sample = sample_name,
          
          timepoint = timepoint
        ) %>%
        
        filter(
          !is.na(count),
          !is.na(pos),
          count >= MIN_READS
        )
      
      
      if (nrow(dat) == 0) {
        next
      }
      
      
      all_data[[length(all_data) + 1]] <- dat
    }
  }
  
  
  result <- bind_rows(
    all_data
  )
  
  
  if (nrow(result) == 0) {
    
    stop(
      "\nNo data remaining for ",
      sample_name,
      "\n"
    )
  }
  
  
  cat(
    "Records:",
    nrow(result),
    "\n"
  )
  
  cat(
    "Genes:",
    n_distinct(result$gene),
    "\n"
  )
  
  cat(
    "Gene × TSS-type combinations:",
    n_distinct(
      paste(
        result$gene,
        result$type,
        sep = "___"
      )
    ),
    "\n"
  )
  
  cat(
    "Total reads:",
    sum(result$count),
    "\n"
  )
  
  
  return(result)
}


# ============================================================
# 2. LOAD 5h AND 24h
# ============================================================

df_5h <- load_sample(
  sample_5h,
  "5h"
)

df_24h <- load_sample(
  sample_24h,
  "24h"
)


# ============================================================
# 3. BASIC DEPTH SUMMARY
# ============================================================

depth_summary <- bind_rows(
  
  df_5h %>%
    summarise(
      timepoint = "5h",
      records = n(),
      genes = n_distinct(gene),
      gene_TSS_pairs = n_distinct(
        paste(gene, type, sep = "___")
      ),
      total_reads = sum(count)
    ),
  
  df_24h %>%
    summarise(
      timepoint = "24h",
      records = n(),
      genes = n_distinct(gene),
      gene_TSS_pairs = n_distinct(
        paste(gene, type, sep = "___")
      ),
      total_reads = sum(count)
    )
)

write_tsv(
  depth_summary,
  file.path(
    plot_base,
    "TSS_Depth_Summary.tsv"
  )
)

cat("\nDepth summary:\n")
print(depth_summary)



classify_gene_type <- function(dat) {
  
  dat %>%
    
    group_by(
      gene,
      type
    ) %>%
    
    summarise(
      
      total_reads = sum(count),
      
      min_pos = min(pos),
      
      max_pos = max(pos),
      
      span_bp =
        max(pos) - min(pos),
      
      unique_positions =
        n_distinct(pos),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      class = case_when(
        
        span_bp <= 3 ~
          "Sharp",
        
        span_bp > 3 &
          span_bp <= 20 ~
          "Narrow",
        
        span_bp > 20 ~
          "Broad",
        
        TRUE ~
          NA_character_
      )
    )
}


gene_type_5h <-
  classify_gene_type(
    df_5h
  ) %>%
  
  rename(
    
    reads_5h = total_reads,
    span_5h = span_bp,
    class_5h = class
  )


gene_type_24h <-
  classify_gene_type(
    df_24h
  ) %>%
  
  rename(
    
    reads_24h = total_reads,
    
    span_24h = span_bp,
    
    class_24h_original = class
  )



original_comparison <-
  gene_type_5h %>%
  
  inner_join(
    
    gene_type_24h,
    
    by = c(
      "gene",
      "type"
    )
  )


write_tsv(
  original_comparison,
  file.path(
    plot_base,
    "TSS_Original_Gene_Type_Comparison.tsv"
  )
)



original_switches <-
  original_comparison %>%
  
  filter(
    !is.na(class_5h),
    !is.na(class_24h_original),
    class_5h != class_24h_original
  ) %>%
  
  mutate(
    
    n_classes = 2,
    
    classes_observed =
      paste(
        class_5h,
        class_24h_original,
        sep = " -> "
      )
  ) %>%
  
  select(
    
    gene,
    
    TSS_type = type,
    
    n_classes,
    
    classes_observed,
    
    class_5h,
    
    class_24h_original,
    
    span_5h,
    
    span_24h,
    
    reads_5h,
    
    reads_24h
  )


write_tsv(
  original_switches,
  file.path(
    plot_base,
    "TSS_Original_Class_Switches.tsv"
  )
)


# ============================================================
# 9. FIND COMMON GENE × TSS-TYPE PAIRS
# ============================================================

pairs_5h <-
  df_5h %>%
  distinct(
    gene,
    type
  )


pairs_24h <-
  df_24h %>%
  distinct(
    gene,
    type
  )


common_pairs <-
  inner_join(
    
    pairs_5h,
    
    pairs_24h,
    
    by = c(
      "gene",
      "type"
    )
  )


cat(
  "\nCommon gene × TSS-type pairs:",
  nrow(common_pairs),
  "\n"
)


write_tsv(
  common_pairs,
  file.path(
    plot_base,
    "TSS_Common_Gene_Type_Pairs.tsv"
  )
)



depth_match_gene_type <- function(
    data_24h,
    target_depth
) {
  
  total_24h <-
    sum(
      data_24h$count
    )
  
  

  if (
    total_24h <= target_depth
  ) {
    
    return(
      
      data_24h %>%
        mutate(
          sampled_count = count
        )
    )
  }
  
  
 
  
  probabilities <-
    data_24h$count /
    total_24h
  
  
  sampled_indices <-
    sample(
      
      seq_len(
        nrow(data_24h)
      ),
      
      size = target_depth,
      
      replace = TRUE,
      
      prob = probabilities
    )
  
  
  sampled_table <-
    table(
      sampled_indices
    )
  
  
  result <-
    data_24h %>%
    mutate(
      sampled_count = 0
    )
  
  
  idx <-
    as.integer(
      names(sampled_table)
    )
  
  
  result$sampled_count[idx] <-
    as.numeric(
      sampled_table
    )
  
  
  result
}




classify_matched <- function(dat) {
  
  dat %>%
    
    filter(
      sampled_count > 0
    ) %>%
    
    group_by(
      gene,
      type
    ) %>%
    
    summarise(
      
      reads_24h_matched =
        sum(sampled_count),
      
      min_pos_24h =
        min(pos),
      
      max_pos_24h =
        max(pos),
      
      span_24h_matched =
        max(pos) -
        min(pos),
      
      unique_positions_24h =
        n_distinct(pos),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      class_24h_matched =
        case_when(
          
          span_24h_matched <= 3 ~
            "Sharp",
          
          span_24h_matched > 3 &
            span_24h_matched <= 20 ~
            "Narrow",
          
          span_24h_matched > 20 ~
            "Broad",
          
          TRUE ~
            NA_character_
        )
    )
}


# ============================================================
# RUN 100 DEPTH-MATCHING ITERATIONS
# ============================================================

all_iterations <-
  vector(
    "list",
    N_ITERATIONS
  )


for (
  iteration in seq_len(N_ITERATIONS)
) {
  
  cat(
    "Depth-matching iteration",
    iteration,
    "/",
    N_ITERATIONS,
    "\n"
  )
  
  
  matched_list <-
    vector(
      "list",
      nrow(common_pairs)
    )
  
  
  for (
    i in seq_len(nrow(common_pairs))
  ) {
    
    gene_i <-
      common_pairs$gene[i]
    
    type_i <-
      common_pairs$type[i]
    
    
    # --------------------------------------------------------
    # 5h data for this EXACT gene × TSS type
    # --------------------------------------------------------
    
    data_5h_pair <-
      df_5h %>%
      
      filter(
        gene == gene_i,
        type == type_i
      )
    
    
    # --------------------------------------------------------
    # 24h data for this EXACT gene × TSS type
    # --------------------------------------------------------
    
    data_24h_pair <-
      df_24h %>%
      
      filter(
        gene == gene_i,
        type == type_i
      )
    
    
    if (
      nrow(data_5h_pair) == 0 ||
      nrow(data_24h_pair) == 0
    ) {
      next
    }
    
    
    # --------------------------------------------------------
    # 5h depth is the target
    # --------------------------------------------------------
    
    target_depth <-
      sum(
        data_5h_pair$count
      )
    
    
    # --------------------------------------------------------
    # Depth-match 24h
    # --------------------------------------------------------
    
    matched <-
      depth_match_gene_type(
        
        data_24h_pair,
        
        target_depth
      )
    
    
    matched_list[[i]] <-
      matched
  }
  
  
  matched_24h <-
    bind_rows(
      matched_list
    )
  
  
  if (
    nrow(matched_24h) == 0
  ) {
    
    warning(
      "No matched data in iteration ",
      iteration
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # Classify matched 24h
  # ----------------------------------------------------------
  
  iteration_result <-
    classify_matched(
      matched_24h
    ) %>%
    
    mutate(
      iteration = iteration
    )
  
  
  all_iterations[[iteration]] <-
    iteration_result
}



matched_results <-
  bind_rows(
    all_iterations
  )


if (
  nrow(matched_results) == 0
) {
  
  stop(
    "\nNo depth-matched results were generated.\n"
  )
}



matched_comparison <-
  matched_results %>%
  
  inner_join(
    
    gene_type_5h %>%
      
      select(
        
        gene,
        
        type,
        
        class_5h,
        
        span_5h,
        
        reads_5h
      ),
    
    by = c(
      "gene",
      "type"
    )
  )


matched_comparison <-
  matched_comparison %>%
  
  mutate(
    
    switched =
      !is.na(class_5h) &
      !is.na(class_24h_matched) &
      class_5h != class_24h_matched,
    
    
    transition =
      case_when(
        
        class_5h == "Sharp" &
          class_24h_matched == "Narrow" ~
          "Sharp -> Narrow",
        
        class_5h == "Sharp" &
          class_24h_matched == "Broad" ~
          "Sharp -> Broad",
        
        class_5h == "Narrow" &
          class_24h_matched == "Sharp" ~
          "Narrow -> Sharp",
        
        class_5h == "Narrow" &
          class_24h_matched == "Broad" ~
          "Narrow -> Broad",
        
        class_5h == "Broad" &
          class_24h_matched == "Sharp" ~
          "Broad -> Sharp",
        
        class_5h == "Broad" &
          class_24h_matched == "Narrow" ~
          "Broad -> Narrow",
        
        TRUE ~
          "No change"
      )
  )



write_tsv(
  matched_comparison,
  file.path(
    plot_base,
    "TSS_Depth_Matched_All_Iterations.tsv"
  )
)



gene_type_transition_summary <-
  matched_comparison %>%
  
  count(
    gene,
    type,
    class_5h,
    transition,
    name = "iterations"
  ) %>%
  
  group_by(
    gene,
    type
  ) %>%
  
  mutate(
    
    fraction_iterations =
      iterations /
      N_ITERATIONS
  ) %>%
  
  ungroup()


write_tsv(
  gene_type_transition_summary,
  file.path(
    plot_base,
    "TSS_Depth_Matched_Transition_Frequencies.tsv"
  )
)


# ============================================================
# 18. FIND MOST COMMON TRANSITION
# ============================================================

best_transition <-
  gene_type_transition_summary %>%
  
  group_by(
    gene,
    type
  ) %>%
  
  slice_max(
    
    order_by =
      fraction_iterations,
    
    n = 1,
    
    with_ties = FALSE
  ) %>%
  
  ungroup()



median_matched_values <-
  matched_results %>%
  
  group_by(
    gene,
    type
  ) %>%
  
  summarise(
    
    median_span_24h_matched =
      median(
        span_24h_matched,
        na.rm = TRUE
      ),
    
    median_reads_24h_matched =
      median(
        reads_24h_matched,
        na.rm = TRUE
      ),
    
    median_unique_positions_24h =
      median(
        unique_positions_24h,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )



robust_switches <-
  best_transition %>%
  
  filter(
    
    transition != "No change",
    
    fraction_iterations >=
      SWITCH_THRESHOLD
  ) %>%
  
  mutate(
    
    n_classes = 2,
    
    classes_observed =
      transition
  ) %>%
  
  
  left_join(
    
    gene_type_5h %>%
      
      select(
        
        gene,
        
        type,
        
        class_5h,
        
        span_5h,
        
        reads_5h
      ),
    
    by = c(
      "gene",
      "type"
    )
  ) %>%
  
  
  left_join(
    
    median_matched_values,
    
    by = c(
      "gene",
      "type"
    )
  ) %>%
  
  
  select(
    
    gene,
    
    TSS_type = type,
    
    n_classes,
    
    classes_observed,
    
    class_5h,
    
    span_5h,
    
    median_span_24h_matched,
    
    reads_5h,
    
    median_reads_24h_matched,
    
    median_unique_positions_24h,
    
    iterations,
    
    fraction_iterations
  ) %>%
  
  
  arrange(
    desc(fraction_iterations),
    gene,
    TSS_type
  )




write_tsv(
  
  robust_switches,
  
  file.path(
    plot_base,
    "TSS_Depth_Matched_Class_Switches.tsv"
  )
)




all_gene_type_results <-
  best_transition %>%
  
  left_join(
    
    gene_type_5h %>%
      
      select(
        
        gene,
        
        type,
        
        class_5h,
        
        span_5h,
        
        reads_5h
      ),
    
    by = c(
      "gene",
      "type"
    )
  ) %>%
  
  left_join(
    
    median_matched_values,
    
    by = c(
      "gene",
      "type"
    )
  ) %>%
  
  mutate(
    
    switched =
      transition != "No change",
    
    robust =
      switched &
      fraction_iterations >=
      SWITCH_THRESHOLD,
    
    n_classes =
      if_else(
        switched,
        2L,
        1L
      ),
    
    classes_observed =
      if_else(
        switched,
        transition,
        class_5h
      )
  ) %>%
  
  select(
    
    gene,
    
    TSS_type = type,
    
    n_classes,
    
    classes_observed,
    
    class_5h,
    
    span_5h,
    
    median_span_24h_matched,
    
    reads_5h,
    
    median_reads_24h_matched,
    
    median_unique_positions_24h,
    
    iterations,
    
    fraction_iterations,
    
    switched,
    
    robust
  ) %>%
  
  arrange(
    gene,
    TSS_type
  )

s
write_tsv(
  
  all_gene_type_results,
  
  file.path(
    plot_base,
    "TSS_Depth_Matched_All_Gene_Type_Results.tsv"
  )
)



direction_summary <-
  robust_switches %>%
  
  count(
    
    TSS_type,
    
    classes_observed,
    
    name = "gene_type_pair_count"
  ) %>%
  
  arrange(
    
    TSS_type,
    
    desc(gene_type_pair_count)
  )


write_tsv(
  
  direction_summary,
  
  file.path(
    plot_base,
    "TSS_Depth_Matched_Switch_Direction_Summary.tsv"
  )
)




switch_summary_by_type <-
  all_gene_type_results %>%
  
  group_by(
    TSS_type
  ) %>%
  
  summarise(
    
    gene_type_pairs =
      n(),
    
    switched_pairs =
      sum(switched),
    
    robust_switched_pairs =
      sum(robust),
    
    fraction_robust_switches =
      mean(
        robust
      ),
    
    .groups = "drop"
  )


write_tsv(
  
  switch_summary_by_type,
  
  file.path(
    plot_base,
    "TSS_Depth_Matched_Summary_by_TSS_Type.tsv"
  )
)


# ============================================================
# RANSITION MATRIX
# ============================================================

transition_matrix <-
  matched_comparison %>%
  
  count(
    
    type,
    
    class_5h,
    
    class_24h_matched,
    
    name = "iterations"
  ) %>%
  
  arrange(
    
    type,
    
    class_5h,
    
    class_24h_matched
  )


write_tsv(
  
  transition_matrix,
  
  file.path(
    plot_base,
    "TSS_Depth_Matched_Transition_Matrix.tsv"
  )
)




class_distribution_5h <-
  gene_type_5h %>%
  
  count(
    
    type,
    
    class_5h,
    
    name = "gene_type_pair_count"
  ) %>%
  
  arrange(
    
    type,
    
    class_5h
  )


write_tsv(
  
  class_distribution_5h,
  
  file.path(
    plot_base,
    "TSS_5h_Class_Distribution_by_Type.tsv"
  )
)


# ============================================================
# PRINT FINAL RESULTS
# ============================================================

cat(
  "\n5h total reads:",
  sum(df_5h$count),
  "\n"
)


cat(
  "24h total reads:",
  sum(df_24h$count),
  "\n"
)


cat(
  "\nCommon gene × TSS-type pairs:",
  nrow(common_pairs),
  "\n"
)


cat(
  "\nRobust switching gene × TSS-type pairs:",
  nrow(robust_switches),
  "\n"
)


cat(
  "\nRobust switch directions:\n"
)


if (
  nrow(direction_summary) > 0
) {
  
  print(
    direction_summary
  )
  
} else {
  
  cat(
    "No robust switches detected.\n"
  )
}


cat(
  "\n\nSummary by TSS type:\n"
)


print(
  switch_summary_by_type
)


cat(
  "\n\nRobust switching gene × TSS-type pairs:\n"
)


if (
  nrow(robust_switches) > 0
) {
  
  print(
    robust_switches
  )
  
} else {
  
  cat(
    "No robust switching pairs detected.\n"
  )
}


cat(
  "\n\nOutput directory:\n",
  plot_base,
  "\n"
)


cat(
  "\nMain result:\n",
  file.path(
    plot_base,
    "TSS_Depth_Matched_Class_Switches.tsv"
  ),
  "\n"
)


cat(
  "\n====================================================\n"
)
