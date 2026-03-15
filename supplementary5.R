
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(ggbeeswarm)

# 2. Define Paths
# ---------------
# Input folder containing the 'polyA_output' subdirectories
indir <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output"

# Output folder for the new plots
plot_base <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/polyA_beeswarm_final"

# Path to your supplementary gene list
excel_path <- "Desktop/Desktop - Chorong’s MacBook Pro/MUSV/MUSV_3/sup/supplementary table2.xlsx"


dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)

# Color scheme (Rep 1-4)
rep_colors <- c(
  "Rep 1" = "#404080", 
  "Rep 2" = "#8E44AD", 
  "Rep 3" = "#2980B9", 
  "Rep 4" = "#69b3a2"
)

# --- 2. Load Gene List ---
if (!file.exists(excel_path)) stop("ERROR: Excel file not found.")

target_genes_df <- read_excel(excel_path, sheet = "Sheet1") %>%
  select(Gene) %>%
  distinct() %>%
  mutate(Gene_Num = as.numeric(str_extract(Gene, "\\d+"))) %>%
  arrange(Gene_Num)

ordered_genes <- target_genes_df$Gene

# --- 3. Data Collection with Strict Timepoint Mapping ---
sample_dirs <- list.dirs(indir, recursive = FALSE, full.names = TRUE)
combined_data <- data.frame()

for (path in sample_dirs) {
  sample_id <- basename(path)
  
  # Strict mapping of folder strings to Timepoint and Replicate
  # We use temporary variables to store the results of the case_when
  mapping <- case_when(
    str_detect(sample_id, "10-04-01") ~ "24h_Rep 1",
    str_detect(sample_id, "11-04")    ~ "24h_Rep 2",
    str_detect(sample_id, "11-05")    ~ "24h_Rep 3",
    str_detect(sample_id, "12-04")    ~ "24h_Rep 4",
    str_detect(sample_id, "12-01")    ~ "5h_Rep 1",
    str_detect(sample_id, "12-02")    ~ "5h_Rep 2",
    str_detect(sample_id, "12-03")    ~ "5h_Rep 3",
    TRUE ~ "Other"
  )
  
  if(mapping == "Other") next
  
  # Split the mapping string back into TP and REP
  current_tp <- str_split(mapping, "_")[[1]][1]
  current_rep <- str_split(mapping, "_")[[1]][2]
  
  polyA_folder <- file.path(path, "polyA_length_per_gene")
  if(!dir.exists(polyA_folder)) next
  
  files <- list.files(polyA_folder, pattern = "_polyA_lengths.txt$", full.names = TRUE)
  
  for (f in files) {
    gene_name <- str_split(basename(f), "_")[[1]][1]
    if (!(gene_name %in% ordered_genes)) next
    
    lines <- readLines(f)
    idx <- grep("leader length counts:", lines, ignore.case = TRUE)
    if (length(idx) == 0) next
    
    df <- tryCatch({
      read.table(text = lines[(idx + 1):length(lines)], col.names = c("Count", "Length"))
    }, error = function(e) NULL)
    
    if(is.null(df)) next
    
    # Expand and assign labels
    df_expanded <- df[rep(seq_len(nrow(df)), df$Count), ] %>% select(Length)
    df_expanded$Gene <- gene_name
    df_expanded$Timepoint <- current_tp
    df_expanded$Replicate <- current_rep
    
    combined_data <- bind_rows(combined_data, df_expanded)
  }
}

# --- 4. Processing for Plotting ---
# --- 4. Processing for Plotting ---
genes_per_page <- 20 

plot_data_all <- combined_data %>%
  mutate(
    Gene = factor(Gene, levels = ordered_genes),
    Replicate = factor(Replicate, levels = names(rep_colors))
  ) %>%
  group_by(Timepoint, Gene, Replicate) %>%
  mutate(Total_Reads = n()) %>%
  ungroup() %>%
  # Create a Page variable that strictly groups 20 genes at a time
  group_by(Timepoint) %>%
  mutate(Page = (as.numeric(Gene) - 1) %/% genes_per_page + 1) %>%
  ungroup()

plot_data_all$Length_Jitter <- plot_data_all$Length + runif(nrow(plot_data_all), -0.4, 0.4)

# --- 5. Generate Two Separate PDFs (5h and 24h) ---
for (tp_to_plot in c("5h", "24h")) {
  
  tp_subset <- plot_data_all %>% filter(Timepoint == tp_to_plot)
  if(nrow(tp_subset) == 0) next
  
  output_pdf <- file.path(plot_base, paste0("PolyA_Comparison_", tp_to_plot, ".pdf"))
  pdf(output_pdf, width = 12, height = 12)
  
  # Get actual unique page numbers for this timepoint
  pages <- sort(unique(tp_subset$Page))
  
  for (p_idx in pages) {
    page_data <- tp_subset %>% filter(Page == p_idx)
    
    p <- ggplot(page_data, aes(x = Replicate, y = Length)) +
      # Use ncol AND nrow to force a consistent grid size on every page
      facet_wrap(~Gene, ncol = 5, nrow = 4, scales = "fixed") + 
      
      geom_quasirandom(aes(y = Length_Jitter, color = Replicate),
                       method = "smiley", width = 0.3, size = 0.15, alpha = 0.8) +
      
      geom_boxplot(aes(y = Length), fill = NA, color = "black", 
                   width = 0.1, size = 0.3, outlier.shape = NA) +
      
      geom_text(data = page_data %>% select(Gene, Replicate, Total_Reads) %>% distinct(),
                aes(y = 38, label = paste0("n=", Total_Reads)), 
                size = 3.0, color = "black", fontface = "bold", angle = 90, hjust = 1) +
      
      scale_color_manual(values = rep_colors) +
      coord_cartesian(ylim = c(1, 40)) +
      theme_bw(base_size = 9) + 
      labs(title = paste("Poly-A Leader Length:", tp_to_plot, "Replicate Comparison"),
           subtitle = paste("Page", p_idx, "| Numeric Order (VACWR)"),
           x = NULL, y = "Leader Length (nt)") +
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            strip.text = element_text(size = 8, face = "bold"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank())
    
    print(p)
  }
  
  dev.off()
}