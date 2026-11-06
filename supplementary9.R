
library(tidyverse)
library(ggbeeswarm)

# 1. Setup paths
target_paths <- c(
  "5ppp" = "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output/HAP1_VacV_5ppp_19215FL-11-01_S1_L001",
  "RppH_Rep2" = "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output/HAP1_VacV_RppH_24h_2_19215FL-11-04_S1_L001",
  "RppH_24h_1" = "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output/HAP1_VacV_RppH_24h_19215FL-10-04-01_S6_L008",
  "RppH_24h_3" = "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/polyA_output/HAP1_VacV_RppH_24h_3_19215FL-11-05_S2_L001",
  "RppH_24h_4" = "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-12_2/polyA_output/HAP1_VacV_RppH_24h_4_19215FL-12-04_S4_L003"
)

# New color scheme for the three groups
treatment_colors <- c(
  "5ppp" = "#69b3a2", 
  "RppH_Rep2" = "#8E44AD", 
  "RppH_Combined" = "#404080"
)

# 2. Collect and Group Data
combined_expanded_data <- data.frame()

for (i in seq_along(target_paths)) {
  path <- target_paths[i]
  group_name <- names(target_paths)[i]
  polyA_folder <- file.path(path, "polyA_length_per_gene")
  
  if(!dir.exists(polyA_folder)) next
  output_files <- list.files(polyA_folder, pattern = "_polyA_lengths.txt$", full.names = TRUE)
  
  for (f in output_files) {
    lines <- readLines(f)
    idx <- grep("leader length counts:", lines, ignore.case = TRUE)
    if (length(idx) == 0) next
    
    data_lines <- lines[(idx + 1):length(lines)]
    data_lines <- data_lines[nzchar(data_lines)]
    if (length(data_lines) == 0) next
    
    df <- read.table(text = data_lines, col.names = c("Count", "Length"))
    header <- lines[1]
    gene_match <- regmatches(header, regexec("(VACWR[0-9]+)", header))[[1]]
    
    if (length(gene_match) >= 2) {
      gene_id <- gene_match[2]
      df_expanded <- df[rep(seq_len(nrow(df)), df$Count), ] %>% select(Length)
      df_expanded$Gene <- gene_id
      
      # Logic for the three-way grouping:
      # 1. Add to 5ppp if it's the 5ppp sample
      if (group_name == "5ppp") {
        df_expanded$Group <- "5ppp"
        combined_expanded_data <- bind_rows(combined_expanded_data, df_expanded)
      } 
      
      # 2. Add to RppH_Combined if it's any RppH sample
      if (grepl("RppH", group_name)) {
        df_combined <- df_expanded
        df_combined$Group <- "RppH_Combined"
        combined_expanded_data <- bind_rows(combined_expanded_data, df_combined)
        
        # 3. ALSO add to RppH_Rep2 if it is specifically Rep 2
        if (group_name == "RppH_Rep2") {
          df_rep2 <- df_expanded
          df_rep2$Group <- "RppH_Rep2"
          combined_expanded_data <- bind_rows(combined_expanded_data, df_rep2)
        }
      }
    }
  }
}

# --- 3. Sorting and Paging ---
all_unique_genes <- sort(unique(combined_expanded_data$Gene))
genes_per_page <- 20 # Reduced for better visibility of 3 bars

plot_df_final <- combined_expanded_data %>%
  mutate(Gene = factor(Gene, levels = all_unique_genes),
         Group = factor(Group, levels = c("5ppp", "RppH_Rep2", "RppH_Combined"))) %>%
  group_by(Gene, Group) %>%
  mutate(Total_Reads = n()) %>%
  ungroup() %>%
  mutate(Page = (as.numeric(Gene) - 1) %/% genes_per_page + 1)

plot_df_final$Length_Jitter <- plot_df_final$Length + runif(nrow(plot_df_final), -0.4, 0.4)

# --- 4. Plotting to Multi-page PDF ---
final_pdf_path <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/All_Genes_PolyA_3Way_Comparison.pdf"
pdf(final_pdf_path, width = 11, height = 8.5) # Landscape for 3-way bars

num_pages <- max(plot_df_final$Page)

for (i in 1:num_pages) {
  page_data <- plot_df_final %>% filter(Page == i)
  if(nrow(page_data) == 0) next
  
  p <- ggplot(page_data, aes(x = Group, y = Length, fill = Group)) +
    facet_wrap(~Gene, ncol = 5, scales = "fixed") + 
    geom_violin(alpha = 0.2, scale = "width", color = NA) +
    geom_quasirandom(aes(y = Length_Jitter, color = Group),
                     method = "smiley", width = 0.25, size = 0.1, alpha = 0.1) +
    geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.7, color = "black") +
    geom_text(data = page_data %>% select(Gene, Group, Total_Reads) %>% distinct(),
              aes(y = 34, label = paste0("n=", Total_Reads)), 
              size = 1.8, color = "black") +
    scale_fill_manual(values = treatment_colors) +
    scale_color_manual(values = treatment_colors) +
    coord_cartesian(ylim = c(1, 35)) +
    theme_bw(base_size = 9) + 
    labs(title = "Poly-A Leader Comparison: 5ppp vs. RppH Rep 2 vs. Combined RppH",
         subtitle = paste("Page", i, "| Alphabetical Order"),
         x = NULL, y = "Leader Length (nt)") +
    theme(legend.position = "top",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 7, face = "bold"))
  
  print(p)
}

dev.off()
message("Success! 3-way comparison PDF saved.")
