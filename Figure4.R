

library(tidyverse)
library(readxl)

# --- 1. LOAD AND PREP ---
# Ensure you are reading Sheet 2
df <- read_excel("Desktop/Desktop - Chorong’s MacBook Pro/MUSV/MUSV_3/sup/supplementary table2.xlsx", sheet = 2)

# --- 1. PREP DATA (Using the exact names from your snippet) ---
plot_ready_df <- supplementary_table2 %>%
  select(
    Gene, 
    Phase, 
    `5h` = TP_Median_Length_5h, 
    `24h` = TP_Median_Length_24h
  ) %>%
  pivot_longer(
    cols = c(`5h`, `24h`), 
    names_to = "Timepoint", 
    values_to = "Median_Length"
  ) %>%
  mutate(
    Phase = factor(tolower(Phase), levels = c("early", "intermediate", "late")),
    Timepoint = factor(Timepoint, levels = c("5h", "24h")),
    Median_Length = as.numeric(Median_Length)
  )

# --- 2. COLORS (Must be defined for your plot to work) ---
phase_colors <- c(
  "early"        = "#A8DADC", 
  "intermediate" = "#457B9D", 
  "late"         = "#1D3557"
)

# --- 3. YOUR BAR CHART CODE ---
p1_median_bars <- ggplot(plot_ready_df, aes(x = reorder(Gene, Median_Length), 
                                            y = Median_Length, 
                                            fill = Phase)) +
  geom_col(color = "black", size = 0.1, position = position_dodge(0.8)) +
  
  geom_text(aes(label = round(Median_Length, 1)), 
            hjust = -0.2, size = 2.5, fontface = "bold",
            position = position_dodge(0.8)) +
  
  # This section keeps your requested layout
  facet_grid(Phase ~ Timepoint, scales = "free_y", space = "free_y") +
  
  coord_flip() +
  scale_fill_manual(values = phase_colors) +
  theme_minimal(base_size = 11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  labs(
    title = "Combined Median PolyA Lengths per Gene",
    subtitle = "Ranked by length within each infection phase",
    x = "Viral Gene ID",
    y = "Combined Median Tail Length (nt)"
  ) +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

print(p1_median_bars)
out_path <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/polyA_Median_Leader_Length_Grid.pdf"
ggsave(out_path, plot = p1_median_bars, width = 10, height = 15) 



library(tidyverse)
library(readxl)

df <- read_excel("Desktop/Desktop - Chorong’s MacBook Pro/MUSV/MUSV_3/sup/supplementary table3.xlsx", sheet = 3)

# --- 2. DATA CLEANING ---
plot_ready_df <- df %>%
  # Split 'Gene_Timepoint' into 'Gene' and 'Phase'
  separate(Gene_Timepoint, into = c("Gene", "Phase"), sep = "_") %>%
  # Select the new exact column names from your image
  select(
    Gene, 
    Phase, 
    `5h` = TP_Combined_Median_5h, 
    `24h` = TP_Combined_Median_24h
  ) %>%
  # Convert to long format
  pivot_longer(
    cols = c(`5h`, `24h`), 
    names_to = "Timepoint", 
    values_to = "Median_Length"
  ) %>%
  # Ensure factors are set correctly for your facet_grid
  mutate(
    Phase = factor(tolower(Phase), levels = c("early", "intermediate", "late")), # 'interme' based on your image crop
    Timepoint = factor(Timepoint, levels = c("5h", "24h")),
    Median_Length = as.numeric(Median_Length)
  )

# --- 3. COLORS (Clean Nature Palette) ---
phase_colors <- c(
  "early"   = "#A8DADC", 
  "intermediate" = "#457B9D", 
  "late"    = "#1D3557"
)

# --- 4. YOUR BAR CHART CODE ---
p1_median_bars <- ggplot(plot_ready_df, aes(x = reorder(Gene, Median_Length), 
                                            y = Median_Length, 
                                            fill = Phase)) +
  geom_col(color = "black", size = 0.1, position = position_dodge(0.8)) +
  
  geom_text(aes(label = round(Median_Length, 1)), 
            hjust = -0.2, size = 2.5, fontface = "bold",
            position = position_dodge(0.8)) +
  
  # Keeps your requested layout
  facet_grid(Phase ~ Timepoint, scales = "free_y", space = "free_y") +
  
  coord_flip() +
  scale_fill_manual(values = phase_colors) +
  theme_minimal(base_size = 11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  labs(
    title = "Combined Median PolyA Lengths per Gene",
    subtitle = "Ranked by length within each infection phase (Table 3)",
    x = "Viral Gene ID",
    y = "Combined Median Tail Length (nt)"
  ) +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

# --- 5. SAVE AND PRINT ---
print(p1_median_bars)

out_path <- "/Volumes/fsmresfiles/Basic_Sciences/DMI/WalshLab/Cho/enzymeseq_data/19215FL-11_2/rplot/NonpurepolyA_Median_Leader_Length_Grid.pdf"
ggsave(out_path, plot = p1_median_bars, width = 10, height = 15)