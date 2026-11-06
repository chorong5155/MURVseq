library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(ggbeeswarm)

# --- 1. Define Paths ---
indir_15 <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output"
indir_30 <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output_2"
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/polyA_15vs30_comparison"
excel_path <- "Desktop/Desktop - Chorong’s MacBook Pro/MUSV/MUSV_3/sup/supplementary table2.xlsx"

dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load Gene List ---
if (!file.exists(excel_path)) stop("ERROR: Excel file not found.")

target_genes_df <- read_excel(excel_path, sheet = "Sheet1") %>%
  rename(Gene = 1) %>% 
  select(Gene) %>%
  distinct() %>%
  mutate(Gene_Num = as.numeric(str_extract(Gene, "\\d+"))) %>%
  arrange(Gene_Num)

ordered_genes <- target_genes_df$Gene

# --- 3. Data Collection Helper Function ---
load_sample_data <- function(base_dir, label) {
  # Find the folder containing "11-05"
  all_dirs <- list.dirs(base_dir, recursive = FALSE)
  target_path <- all_dirs[str_detect(basename(all_dirs), "11-05")]
  
  if(length(target_path) == 0) return(data.frame())
  
  polyA_folder <- file.path(target_path[1], "polyA_length_per_gene")
  if(!dir.exists(polyA_folder)) return(data.frame())
  
  sample_results <- data.frame()
  
  # List all available files once to speed up matching
  available_files <- list.files(polyA_folder, pattern = "_polyA_lengths.txt$")
  
  for (g in ordered_genes) {
    # Match gene name to start of filename (robust against different suffixes)
    match_file <- available_files[str_starts(available_files, paste0(g, "_"))]
    if(length(match_file) == 0) next
    
    file_path <- file.path(polyA_folder, match_file[1])
    lines <- readLines(file_path)
    idx <- grep("leader length counts:", lines, ignore.case = TRUE)
    if (length(idx) == 0) next
    
    df <- tryCatch({
      read.table(text = lines[(idx + 1):length(lines)], col.names = c("Count", "Length"))
    }, error = function(e) NULL)
    
    if(is.null(df) || nrow(df) == 0) next
    
    # Expand and build dataframe using base R subsetting to avoid select() errors
    df_expanded <- df[rep(seq_len(nrow(df)), df$Count), ]
    df_expanded$Gene <- g
    df_expanded$Method <- label
    
    # Explicitly keep only what we need
    df_expanded <- df_expanded[, c("Length", "Gene", "Method")]
    sample_results <- bind_rows(sample_results, df_expanded)
  }
  return(sample_results)
}

# --- 4. Process and Combine ---
data_15 <- load_sample_data(indir_15, "15mer")
data_30 <- load_sample_data(indir_30, "30mer")

if(nrow(data_15) == 0 && nrow(data_30) == 0) stop("No data found. Check folder names/paths.")

combined_data <- bind_rows(data_15, data_30) %>%
  mutate(
    Gene = factor(Gene, levels = ordered_genes),
    Method = factor(Method, levels = c("15mer", "30mer"))
  )

# Define pages (20 genes per page)
genes_per_page <- 20
plot_data <- combined_data %>%
  group_by(Gene, Method) %>%
  mutate(Total_Reads = n()) %>%
  ungroup() %>%
  mutate(Page = (as.numeric(Gene) - 1) %/% genes_per_page + 1)

# Add jitter for beeswarm
plot_data$Length_Jitter <- plot_data$Length + runif(nrow(plot_data), -0.4, 0.4)

# --- 5. Generate PDF ---
output_pdf <- file.path(plot_base, "Comparison_15mer_vs_30mer_Rep3.pdf")
pdf(output_pdf, width = 12, height = 12)

pages <- sort(unique(plot_data$Page))

for (p_idx in pages) {
  page_data <- plot_data %>% filter(Page == p_idx)
  
  p <- ggplot(page_data, aes(x = Method, y = Length)) +
    facet_wrap(~Gene, ncol = 5, nrow = 4, scales = "fixed") + 
    
    geom_quasirandom(aes(y = Length_Jitter, color = Method),
                     method = "smiley", width = 0.3, size = 0.15, alpha = 0.8) +
    
    geom_boxplot(fill = NA, color = "black", width = 0.1, 
                 size = 0.5, outlier.shape = NA) +
    
    # Display N-count
    geom_text(data = page_data %>% select(Gene, Method, Total_Reads) %>% distinct(),
              aes(y = 38, label = paste0("n=", Total_Reads)), 
              size = 2.8, color = "black", fontface = "bold", angle = 90, hjust = 1) +
    
    scale_color_manual(values = c("15mer" = "#404080", "30mer" = "#D35400")) +
    coord_cartesian(ylim = c(0, 40)) +
    theme_bw(base_size = 10) + 
    labs(title = "Poly-A Leader Length: 15mer vs 30mer Threshold",
         subtitle = paste("Sample 11-05 (Rep 3) | Page", p_idx),
         x = "15mer vs 30mer ", y = "Length (nt)") +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  
  print(p)
}

dev.off()
