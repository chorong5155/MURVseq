library(tidyverse)

# --- 1. Load the 5' Matrix ---
# Note: Ensure the file name matches exactly
data <- read_tsv("/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/merged_5prime_matrix.txt")

# Define the specific sample columns you requested
col_5ppp <- "HAP1_VacV_5ppp_19215FL-11-01_S1_L001_Aligned_properlyPaired.sortedByCoord.out"
col_rpph <- "HAP1_VacV_RppH_2_19215FL-11-04_S1_L001_Aligned_properlyPaired.sortedByCoord.out"

# --- 2. Filter and Clean ---
vacv_matrix <- data %>%
  filter(chr == "chrVacV") %>%
  select(chr, pos, all_of(c(col_5ppp, col_rpph))) %>%
  rename(count_5ppp = !!sym(col_5ppp),
         count_rpph = !!sym(col_rpph)) %>%
  # Filter out positions with zero signal in both to speed up processing
  filter(count_5ppp > 0 | count_rpph > 0)

# --- 3. Identify Top 6 TSS Regions to Show "Messiness" ---
# We find the positions with the highest RppH signal (canonical TSS)
top_tss <- vacv_matrix %>%
  arrange(desc(count_rpph)) %>%
  head(6) 

# Create windows of +/- 30bp around these top peaks
windows <- top_tss %>%
  rowwise() %>%
  do({
    center <- .$pos
    vacv_matrix %>%
      filter(pos >= center - 30, pos <= center + 30) %>%
      mutate(Peak_Label = paste0("Region at ", center, " bp"))
  }) %>%
  ungroup()

# --- 4. Normalize for Shape Comparison ---
# Normalizing by total signal in the window shows the 'jaggedness' shape 
# regardless of sequencing depth differences.
plot_df <- windows %>%
  pivot_longer(cols = starts_with("count"), 
               names_to = "Treatment", 
               values_to = "Count") %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c("count_5ppp", "count_rpph"), 
                            labels = c("5ppp", "RppH"))) %>%
  group_by(Peak_Label, Treatment) %>%
  mutate(Frequency = Count / sum(Count + 0.0001))

# --- 5. Generate Jaggedness Plot ---
p_jagged <- ggplot(plot_df, aes(x = pos, y = Frequency, fill = Treatment)) +
  # Using position="identity" and alpha=0.5 lets you see the overlap
  geom_col(position = "identity", alpha = 0.6, width = 1) +
  facet_wrap(~Peak_Label, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("5ppp" = "#69b3a2", "RppH" = "#404080")) +
  theme_classic(base_size = 12) +
  labs(
    title = "Genomic 5' End Jaggedness: VacV TSS Hotspots",
    subtitle = "Purple (RppH) shows precise peaks; Green (5ppp) shows dispersed, messy starts",
    x = "Genomic Position on chrVacV",
    y = "Relative Frequency"
  ) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Save the plot
print(p_jagged)
ggsave("VacV_Genomic_Jaggedness_Comparison.pdf", plot = p_jagged, width = 12, height = 8)



