
user_lib <- "~/R/4.5.1"
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

pkgs <- c("data.table", "limma", "matrixStats")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

library(data.table)
library(limma)
library(matrixStats)


infile <- "/projects/b1042/WalshLab/cpp4747/19215FL-11/matrix/merged_5prime_matrix.txt"
outdir <- "/projects/b1042/WalshLab/cpp4747/19215FL-11/matrix"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

log_outfile <- file.path(outdir, "merged_5prime_logCPM_batch_corrected.txt")


chunk_size <- 50000
meta_cols <- c("chr","pos")
min_total_counts <- 5   # ⭐ aggressive filtering for sparse 5′ data


header <- fread(infile, nrows=0)
colnames_all <- colnames(header)
count_cols <- setdiff(colnames_all, meta_cols)

samples <- count_cols
batch <- sub(".*_(19215FL-\\d+)_.*", "\\1", samples)


# PASS 1 — GLOBAL LIBRARY SIZES

message("Computing global library sizes...")

lib_sizes <- setNames(rep(0, length(count_cols)), count_cols)

con <- file(infile, "r")
readLines(con, n = 1)  # skip header

repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0) break
  
  dt <- fread(text = lines, header = FALSE)
  setnames(dt, colnames_all)
  
  lib_sizes <- lib_sizes + colSums(dt[, ..count_cols])
}
close(con)

print(lib_sizes)


# PASS 2 — FILTER + BUILD MATRIX

message("Filtering low-signal loci and building matrix...")

keep_rows <- list()
meta_list <- list()

con <- file(infile, "r")
readLines(con, n = 1)

repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0) break
  
  dt <- fread(text = lines, header = FALSE)
  setnames(dt, colnames_all)
  
  rs <- rowSums(dt[, ..count_cols])
  keep <- rs >= min_total_counts
  
  if (any(keep)) {
    keep_rows[[length(keep_rows)+1]] <- as.matrix(dt[keep, ..count_cols])
    meta_list[[length(meta_list)+1]] <- dt[keep, ..meta_cols]
  }
}
close(con)

counts <- do.call(rbind, keep_rows)
meta <- rbindlist(meta_list)

rm(keep_rows, meta_list)
gc()

message("Rows kept: ", nrow(counts))


message("Computing CPM...")

cpm <- t(t(counts) / lib_sizes) * 1e6

rm(counts)
gc()


logCPM <- log2(cpm + 1)

rm(cpm)
gc()

message("Running limma batch correction with biological protection...")

# 1. Create grouping factor
groups <- rep("Mock_24h", length(samples))
groups[grepl("VacV", samples, ignore.case = TRUE)] <- "VacV_24h"
groups[grepl("5h", samples, ignore.case = TRUE)]   <- "VacV_5h"
groups[grepl("mock", samples, ignore.case = TRUE) & grepl("5h", samples)] <- "Mock_5h"

group_factor <- factor(groups)
message("Groups identified: ", paste(levels(group_factor), collapse=", "))

# 2. Create the design matrix
design <- model.matrix(~group_factor)

# 3. Run correction WITH the design argument
# This uses 'logCPM' which still exists because we haven't 'rm()'ed it yet!
logCPM_corrected <- removeBatchEffect(logCPM, batch = batch, design = design)

# NOW it is safe to remove the original logCPM
rm(logCPM)
gc()

# -----------------------------
# SAVE OUTPUT
# -----------------------------
message("Writing output...")
out_dt <- cbind(meta, as.data.table(logCPM_corrected))
fwrite(out_dt, log_outfile, sep="\t", quote=FALSE)

message("DONE — batch-corrected logCPM written.")