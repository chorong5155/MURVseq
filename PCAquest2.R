# -----------------------------
# Setup Library Paths
# -----------------------------
user_lib <- "/home/cpp4747/R/4.5.1"
.libPaths(c(user_lib, .libPaths()))

library(data.table)
library(ggplot2)
library(matrixStats)

infile <- "/projects/b1042/WalshLab/cpp4747/19215FL-11/matrix/merged_5prime_logCPM_batch_corrected.txt"
outdir <- "/projects/b1042/WalshLab/cpp4747/19215FL-11/matrix"

dt <- fread(infile)

# ⭐ Select RppH columns directly
rpph_cols <- grep("RppH", colnames(dt), value = TRUE)
mat <- as.matrix(dt[, ..rpph_cols])
rownames(mat) <- paste(dt$chr, dt$pos, sep=":")

# -------------------------------------------------
# ⭐ NEW FILTERING STRATEGY
# -------------------------------------------------
# Step 1: Filter out very low expression (removes background noise)
# Keep loci with at least a log2(CPM) of 1 (approx 2 CPM) in at least 3 samples
keep_expr <- rowSums(mat > 1) >= 3
mat_filtered <- mat[keep_expr, ]

# Step 2: From those, take the most variable
vars <- rowVars(mat_filtered)
keep_var <- vars > quantile(vars, 0.8) # Keep top 20% of expressed genes
mat_final <- mat_filtered[keep_var, ]

message("Using ", nrow(mat_final), " high-quality loci for PCA")

# -------------------------------------------------
# PCA Calculation
# -------------------------------------------------
# We use scale. = FALSE because log-data is already relatively homoscedastic
pca <- prcomp(t(mat_final), center = TRUE, scale. = TRUE)

pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)


df$Condition <- ifelse(grepl("mock", df$sample, ignore.case=TRUE), "Mock", "VacV")
df$Timepoint <- ifelse(grepl("5h", df$sample), "5h", "24h")

# Extract Replicate number/ID from the sample name (e.g., RppH_1, RppH_2)
df$Replicate <- sub(".*RppH_?(\\d+)?_.*", "\\1", df$sample)
df$Replicate[df$Replicate == ""] <- "1" # Default to 1 if no number found

# Create a unique ID for the legend/labels
df$ID <- paste(df$Condition, df$Timepoint, "Rep", df$Replicate, sep="_")

# -------------------------------------------------
# Visualization
# -------------------------------------------------
# -------------------------------------------------
# Visualization - DEBUG VERSION
# -------------------------------------------------
# 1. Add a tiny bit of "jitter" so overlapping points peek out
# 2. Add sample labels so you can count all 7
p <- ggplot(df, aes(PC1, PC2, color=Condition, shape=Timepoint)) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray") +
  geom_vline(xintercept = 0, linetype="dashed", color="gray") +
  # position_jitter prevents points from hiding exactly behind each other
  # Increase jitter significantly to physically push identical points apart
  geom_point(size=6, alpha=0.6, 
             position = position_jitter(width = 2, height = 2)) + 
  
  # Use ggrepel so the 10 labels don't stack on top of each other
  # Note: You'll need library(ggrepel)
  ggrepel::geom_text_repel(aes(label=sample), size=2, max.overlaps = Inf) +
  # This will label every point with its actual sample name
  geom_text(aes(label=sample), vjust=-1.5, size=3, check_overlap = FALSE) +
  theme_bw(base_size=14) + 
  scale_color_manual(values=c("Mock"="#4477AA","VacV"="#EE6677")) +
  labs(
    title="PCA: 5′ End (Finding the 7 Samples)",
    subtitle=paste0("PC1: ", pct[1], "% | PC2: ", pct[2], "%"),
    x=paste0("PC1 (", pct[1], "%)"),
    y=paste0("PC2 (", pct[2], "%)")
  )

# Increase width to make room for labels
ggsave(file.path(outdir, "PCA_RppH_Debug_Labels.pdf"), p, width=10, height=8)

# -------------------------------------------------
# CRITICAL CHECK: Print the coordinates to the log
# -------------------------------------------------
message("Coordinates for all samples:")
print(df[, c("sample", "PC1", "PC2", "Condition", "Timepoint")])