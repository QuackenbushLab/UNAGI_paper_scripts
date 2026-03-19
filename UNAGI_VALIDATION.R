# Setup
options(stringsAsFactors = FALSE)

library(Matrix)
library(igraph)
library(biomaRt)
library(data.table)  # for faster I/O
library(netZooR)

# Parameters
FAST_MODE <- TRUE

if (FAST_MODE) {
  TAIL_FRAC <- 0.05        # Keep only extreme 5% of edges
  COL_SAMPLE_FRAC <- 0.02  # Sample only 2% of columns for quantiles
  UNAGI_CAP <- 5e4         # Cap number of edges written
  BLOCK_SZ <- 5000         # Larger block size for speed
} else {
  TAIL_FRAC <- 0.05
  COL_SAMPLE_FRAC <- 1.0
  UNAGI_CAP <- 5e5
  BLOCK_SZ <- 1000
}

# Directories
coexpressionDir <- NULL
filteredCoexpressionDir <- NULL
unagiDir <- NULL
dir.create(filteredCoexpressionDir)
dir.create(unagiDir)

# Load sparse matrices
Skin_matrix     <- readRDS(paste0(coexpressionDir, "/skin_S.RDS"))
Skeletal_matrix <- readRDS(paste0(coexpressionDir, "/Skeletal_S.RDS"))
Adipose_matrix  <- readRDS(paste0(coexpressionDir, "/Adipose_S.RDS"))
Lung_matrix     <- readRDS(paste0(coexpressionDir, "/Lung_S.RDS"))
Artery_matrix   <- readRDS(paste0(coexpressionDir, "/Artery_S.RDS"))

ensure_dgc_with_names <- function(mat) {
  mat <- as(mat, "dgCMatrix")
  if (is.null(rownames(mat))) rownames(mat) <- as.character(seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- as.character(seq_len(ncol(mat)))
  mat
}

Skin_matrix     <- ensure_dgc_with_names(Skin_matrix)
Skeletal_matrix <- ensure_dgc_with_names(Skeletal_matrix)
Adipose_matrix  <- ensure_dgc_with_names(Adipose_matrix)
Lung_matrix     <- ensure_dgc_with_names(Lung_matrix)
Artery_matrix   <- ensure_dgc_with_names(Artery_matrix)

# Genes of interest
genes_of_interest <- c(
  "MB", "MYH2", "MYL2", "DES", "TNNC1", "TNNC2", "ENO3", "MYL3",
  "TTN", "TPM1", "TCAP", "MYL1", "TPM3", "TPM4", "COX5B",
  "COX5A", "CKMT2", "TUBA1A", "TUBA1B", "TUBA4A", "TUBA1C",
  "TUBA3C", "TUBA8", "TUBA3D", "TUBA3E", "TUBA4B",
  "PPARG", "FASN", "SREBF1", "SCD", "CEBPA", "ADIPOQ", "FABP4"
)
genes_of_interest2 <- c(
  "MB", "MYH2", "MYL2", "DES", "TNNC1", "TNNC2", "ENO3", "MYL3",
  "TTN", "TPM1", "TCAP", "MYL1", "TPM3", "TPM4", "COX5B",
  "COX5A", "CKMT2", "TUBA1A", "TUBA1B", "TUBA4A", "TUBA1C",
  "TUBA3C", "TUBA8", "TUBA3D", "TUBA3E", "TUBA4B"
)
genes_of_interest3 <- c("PPARG", "FASN", "SREBF1", "SCD", "CEBPA", "ADIPOQ", "FABP4")

# Map gene symbols → Ensembl IDs
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

map_symbols_to_ensembl <- function(symbols) {
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "hgnc_symbol",
    values = symbols,
    mart = mart
  )
  mapping <- subset(mapping, hgnc_symbol != "" & ensembl_gene_id != "")
  list(ids = unique(mapping$ensembl_gene_id), map = mapping)
}

goi1_map <- map_symbols_to_ensembl(genes_of_interest)
goi2_map <- map_symbols_to_ensembl(genes_of_interest2)
goi3_map <- map_symbols_to_ensembl(genes_of_interest3)

# FAST quantile estimation
estimate_quantiles_from_sample_fast <- function(mat, tail_frac = 0.05, col_sample_frac = 1.0, sample_cap = 1e6) {
  if (col_sample_frac < 1) {
    sampled_cols <- sample.int(ncol(mat), size = max(1, round(ncol(mat) * col_sample_frac)), replace = FALSE)
    mat <- mat[, sampled_cols, drop = FALSE]
  }
  nz_vals <- mat@x
  if (length(nz_vals) == 0) {
    return(list(qlo = NA_real_, qhi = NA_real_, constant = TRUE))
  }
  if (length(nz_vals) > sample_cap) {
    nz_vals <- sample(nz_vals, sample_cap)
  }
  qlo <- quantile(nz_vals, probs = tail_frac, type = 7, names = FALSE, na.rm = TRUE)
  qhi <- quantile(nz_vals, probs = 1 - tail_frac, type = 7, names = FALSE, na.rm = TRUE)
  list(qlo = as.numeric(qlo), qhi = as.numeric(qhi), constant = isTRUE(all.equal(qlo, qhi)))
}

# FAST filtering
stream_filter_to_csv_fast <- function(mat, out_csv, qlo, qhi, constant = FALSE,
                                      block_size = 1000, drop_diag = TRUE, write_cap = Inf) {
  edges <- vector("list", ceiling(ncol(mat) / block_size))
  total_written <- 0L
  idx <- 1L
  
  for (c0 in seq.int(1L, ncol(mat), by = block_size)) {
    if (total_written >= write_cap) break
    cols <- c0:min(c0 + block_size - 1L, ncol(mat))
    blk  <- mat[, cols, drop = FALSE]
    if (length(blk@x) == 0L) next
    
    row_ids <- blk@i + 1L
    col_ids <- rep.int(cols, diff(blk@p))
    vals    <- blk@x
    
    if (drop_diag) {
      keep_nd <- row_ids != col_ids
      row_ids <- row_ids[keep_nd]
      col_ids <- col_ids[keep_nd]
      vals    <- vals[keep_nd]
    }
    
    keep <- if (constant || is.na(qlo) || is.na(qhi)) rep(TRUE, length(vals)) else (vals <= qlo) | (vals >= qhi)
    sel <- which(keep)
    remain <- write_cap - total_written
    if (length(sel) > remain) sel <- sel[seq_len(remain)]
    
    if (length(sel) > 0) {
      edges[[idx]] <- data.table(
        source = rownames(mat)[row_ids[sel]],
        target = colnames(mat)[col_ids[sel]],
        score  = vals[sel]
      )
      total_written <- total_written + nrow(edges[[idx]])
      idx <- idx + 1L
    }
  }
  
  edges <- rbindlist(edges)
  fwrite(edges, out_csv)
  total_written
}

read_top_n_from_csv <- function(path, n) {
  if (!file.exists(path)) return(data.frame(source=character(), target=character(), score=numeric()))
  fread(path, nrows = if (is.finite(n)) n else -1)
}

# Main runner
build_and_run_unagi_from_matrix <- function(mat, goi_map, out_filtered_csv, out_unagi_csv,
                                            tail_frac = 0.05, block_size = 1000, col_sample_frac = 1.0,
                                            unagi_max_edges = 5e5, filtered_tmp_csv = tempfile(fileext = ".csv")) {
  message("\n===== START: ", out_unagi_csv, " =====")
  
  #Quantile estimation
  message("Estimating quantiles")
  t1 <- Sys.time()
  qs <- estimate_quantiles_from_sample_fast(mat, tail_frac, col_sample_frac)
  message("  → Quantiles: qlo=", signif(qs$qlo, 4), " | qhi=", signif(qs$qhi, 4))
  message("  → Estimation time: ", signif(difftime(Sys.time(), t1, units="secs"), 3), " sec")
  
  #Filtering
  message("Filtering edges")
  t2 <- Sys.time()
  total_written <- stream_filter_to_csv_fast(mat, filtered_tmp_csv, qs$qlo, qs$qhi, qs$constant,
                                             block_size, TRUE, unagi_max_edges)
  message("  → Edges kept: ", total_written)
  message("  → Filtering time: ", signif(difftime(Sys.time(), t2, units="secs"), 3), " sec")
  
  #Run UNAGI
  message("Running UNAGI")
  net_df <- read_top_n_from_csv(filtered_tmp_csv, unagi_max_edges)
  goi_present <- intersect(goi_map$ids, unique(c(net_df$source, net_df$target)))
  
  if (length(goi_present) == 0L) {
    saveRDS(net_df, out_filtered_csv)
    return(net_df)
  }
  
  subnetwork <- tryCatch({
    RunUNAGI(nodeSet = goi_present, network = net_df, hopConstraint = 2, verbose = TRUE, topX = NULL)
  }, error = function(e) {
    message("UNAGI error: ", e$message)
    data.frame()
  })
  
  if (nrow(subnetwork) > 0) {
    fwrite(subnetwork, out_unagi_csv)
    return(subnetwork)
  } else {
    fwrite(net_df, out_unagi_csv)
    return(net_df)
  }
}

# Run all tissues
subnetwork1 <- build_and_run_unagi_from_matrix(Skin_matrix,     goi1_map, paste0(filteredCoexpressionDir, "Skin_sm_S.RDS"),    
                                               paste0(unagiDir, "skin.csv"),
                                               tail_frac = TAIL_FRAC, block_size = BLOCK_SZ, col_sample_frac = COL_SAMPLE_FRAC, unagi_max_edges = UNAGI_CAP)
subnetwork2 <- build_and_run_unagi_from_matrix(Skeletal_matrix, goi1_map, paste0(filteredCoexpressionDir, "Skeletal_Muscle_sm_S.RDS"),
                                               paste0(unagiDir, "skeletal.csv"),
                                               tail_frac = TAIL_FRAC, block_size = BLOCK_SZ, col_sample_frac = COL_SAMPLE_FRAC, unagi_max_edges = UNAGI_CAP)
subnetwork3 <- build_and_run_unagi_from_matrix(Adipose_matrix,  goi1_map, paste0(filteredCoexpressionDir, "Adipose_Subcutaneous_sm_S.RDS"),
                                               paste0(unagiDir, "subcutaneous.csv"),
                                               tail_frac = TAIL_FRAC, block_size = BLOCK_SZ, col_sample_frac = COL_SAMPLE_FRAC, unagi_max_edges = UNAGI_CAP)
subnetwork4 <- build_and_run_unagi_from_matrix(Lung_matrix,     goi1_map, paste0(filteredCoexpressionDir, "Lung_sm_S.RDS"),     
                                               paste0(unagiDir, "lung.csv"),
                                               tail_frac = TAIL_FRAC, block_size = BLOCK_SZ, col_sample_frac = COL_SAMPLE_FRAC, unagi_max_edges = UNAGI_CAP)
subnetwork5 <- build_and_run_unagi_from_matrix(Artery_matrix,   goi1_map, paste0(filteredCoexpressionDir, "Artery_Aorta_sm_S.RDS"),
                                               paste0(unagiDir, "aorta.csv"),
                                               tail_frac = TAIL_FRAC, block_size = BLOCK_SZ, col_sample_frac = COL_SAMPLE_FRAC, unagi_max_edges = UNAGI_CAP)

# Summary function
summarize_tissue <- function(net_df, goi_map, label) {
  cat("\n=== ", label, " ===\n", sep = "")
  if (nrow(net_df) == 0) {
    cat("No edges.\n")
    return()
  }
  found_ids <- intersect(goi_map$ids, unique(c(net_df$source, net_df$target)))
  found_symbols <- goi_map$map$hgnc_symbol[match(found_ids, goi_map$map$ensembl_gene_id)]
  cat("GOIs found: ", length(found_ids), "/", length(goi_map$ids), "\n")
  if (length(found_ids) > 0) {
    print(data.frame(Ensembl_ID = found_ids, Gene_Symbol = found_symbols))
  }
  cat("\nHead of subnetwork:\n")
  print(head(net_df))
}

summarize_tissue(subnetwork1, goi1_map, "Skin")
summarize_tissue(subnetwork2, goi2_map, "Skeletal")
summarize_tissue(subnetwork3, goi3_map, "Adipose")
summarize_tissue(subnetwork4, goi1_map, "Lung")
summarize_tissue(subnetwork5, goi1_map, "Artery")
