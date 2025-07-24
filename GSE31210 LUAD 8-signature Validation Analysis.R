################################################################################
# Title: Validate 8-Gene LUAD Signature in GSE31210 (microarray)
# Goal: Load, process and assess GSE31210 for assessing 8-gene signature in dataset
# Dataset: GSE31210 (Affymetrix Human Genome U133A Array, GPL96)
# Objective: Assess survival stratification using 8-gene signature
################################################################################

# Setup
rm(list = ls())
options(stringsAsFactors = FALSE)
setwd('/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Nature Scientific Reports/Analyses/')

# Load libraries
library(survival)
library(survminer)
library(tidyverse)
library(pROC)
library(ggplot2)
library(pheatmap)
library(reportROC)
library(verification)

# Define file paths
series_file <- "GSE31210_series_matrix.txt"
annot_file <- "GPL570-55999.txt"
################################################################################

# Part 1: 
# 1. Extract sample IDs and survival info from raw lines
raw_lines <- readLines(series_file)

sample_ids <- gsub('"', '', strsplit(raw_lines[29], "\t")[[1]][-1])

# Vital_status (Alive or Dead) → create binary event indicator (1 = Dead, 0 = Alive)
surv_status_raw50 <- strsplit(raw_lines[50], "\t")[[1]][-1]
surv_status_raw51 <- strsplit(raw_lines[51], "\t")[[1]][-1]
surv_status_raw52 <- strsplit(raw_lines[52], "\t")[[1]][-1]

status_matrix <- data.frame(
  Line50 = surv_status_raw50,
  Line51 = surv_status_raw51,
  Line52 = surv_status_raw52,
  stringsAsFactors = FALSE, row.names = sample_ids
)

# Create event column: 1 if any row for that sample contains "Dead", 0 if only "Alive", NA otherwise
status_flags <- apply(status_matrix, 1, function(x) {
  if (any(grepl("Dead", x, ignore.case = TRUE))) {
    return(1)
  } else if (any(grepl("Alive", x, ignore.case = TRUE))) {
    return(0)
  } else {
    return(NA)
  }
})

# Doing just to be safe but if you check the numbers its likely not needed (can just take the numbers from line 52)
extract_exact_days <- function(x) {
  # Remove quotes
  x_clean <- gsub('"', '', x)
  
  # Look for exact match to "days before death/censor: <number>"
  match <- grep("^days before death/censor: ", x_clean, value = TRUE)
  
  if (length(match) == 1) {
    return(as.numeric(gsub("^days before death/censor: ", "", match)))
  } else {
    return(NA)
  }
}

status_matrix$`days before death/censor` <- apply(status_matrix[, 1:3], 1, extract_exact_days)
# Combine with status_matrix for clarity (optional)
status_matrix <- as.data.frame(status_matrix)
status_matrix$Event <- status_flags

pheno <- data.frame(
  SampleID = rownames(status_matrix),
  Event = status_matrix$Event,
  Time = status_matrix$`days before death/censor`,
  stringsAsFactors = FALSE
) # edit above raw files if Age or other pheno data is needed

# 2. Load and process expression data
expr <- read.delim(series_file, skip = 75, header = TRUE, check.names = FALSE)
expr <- expr[!grepl("^!", expr$ID_REF), ] #just removes last row

# 3. Load annotation and map probe IDs to gene symbols
annot <- read.delim(annot_file, header = TRUE, sep = "\t", comment.char = "#", quote = "")
probe_map <- annot[, c("ID", "Gene.Symbol")]
colnames(probe_map) <- c("ID_REF", "Gene")
probe_map$Gene <- gsub("///", "|", probe_map$Gene, fixed = TRUE)

# Merge expression with gene symbols
expr_annot <- merge(probe_map, expr, by = "ID_REF")
expr_annot <- expr_annot[expr_annot$Gene != "", ]

# Collapse to gene-level (mean across probes mapping to the same gene)
dim(expr_annot) #21225, 464
length(unique(expr_annot$Gene)) #13515
expr_gene <- expr_annot %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
  ungroup()

# Format matrix
expr_mat <- as.data.frame(expr_gene)
rownames(expr_mat) <- expr_mat$Gene
expr_mat <- expr_mat[, -1]
################################################################################

# Dataset cleaning & cluster 3.0 normlization
expr_mat[expr_mat == 0] <- NA # Replace zero values with NA to prepare for log2 transformation

# Remove rows with more than 50% & 90% missing values
cleanDat.removeRows.50percent.idx <- which(apply(expr_mat, 1, function(x) length(which(is.na(x)))) > ncol(expr_mat) / 2)
cleanDat.removeRows.90percent.idx <- which(apply(expr_mat, 1, function(x) length(which(is.na(x)))) > ncol(expr_mat) * 0.90)
print(c(cleanDat.removeRows.50percent.idx,cleanDat.removeRows.90percent.idx))

# No values w/missingness as expected in array data
#cleanDat.upTo50percentNA <- expr_mat[-cleanDat.removeRows.50percent.idx, ]
#cleanDat.upTo50percentNA <- log2(cleanDat.upTo50percentNA)  # Log-transform
#cleanDat.upTo90percentNA <- expr_mat[-cleanDat.removeRows.90percent.idx, ]
#cleanDat.upTo90percentNA <- log2(cleanDat.upTo90percentNA)  # Log-transform

# Equivalent Cluster 3.0 method (log2 was done in missing code for CPTACT data)
expr_mat_log <- log2(expr_mat)  # Log-transform
# Normalization step 1: Median subtraction for each sample column
norm.expr_mat <- apply(expr_mat_log, 2, function(x) x - median(x, na.rm = TRUE))
# Normalization step 2: Scale each column so the sum of the squares is 1
norm.scaled.expr_mat <- apply(norm.expr_mat, 2, function(x) x / sqrt(sum(x^2, na.rm = TRUE)))

# Create PDF to visualize histograms (didn't create)
pdf("Histograms-Normalized-Scaled-exper_mat.pdf", width = 16, height = 7)
par(mfrow = c(1, 2))  # Side-by-side plots

# Histogram for up to 90% missing data
h<-hist(norm.scaled.expr_mat, breaks = 250, main = "Up to 90% missing log2(FPKM) values", xlab = "Expression Value", col = "lightblue")
max(h$counts[0:233]) #63491
which(h$counts==63491) #154
h$mids[154] #0.0039
h$mids[162] #0.0055
abline(v=0.0039,col="red")

nonNoisePeak.shift <- 0.0039
#norm.tumor.exp.scaled.recentered<-2^(norm.scaled.upTo90percentNA - nonNoisePeak.shift) 
norm.tumor.exp.scaled.recentered<-2^(norm.scaled.expr_mat) 

hist(norm.tumor.exp.scaled.recentered, breaks = 250, main = "Up to 90% missing unlogged shifted log2(FPKM) values", xlab = "Expression Value", col = "lightblue") #this is a histogram of an distribution of an unlogged non-noise peak center at 1

# Histogram for up to 50% missing data
# hist(norm.scaled.upTo50percentNA, breaks = 250, main = "Up to 50% missing log2(FPKM) values", xlab = "Expression Value", col = "lightgreen")

dev.off()
################################################################################

# Subset to 8-gene signature (SVBP = CCDC23)
dim(norm.tumor.exp.scaled.recentered) #23520, 246

signature_genes <- c("ATP6V0E1", "CCDC23", "HSDL1", "UBTD1", 
                     "GNPNAT1", "XRCC2", "TFAP2A", "PPP1R13L")

signature_expr <- norm.tumor.exp.scaled.recentered[rownames(norm.tumor.exp.scaled.recentered) %in% signature_genes, ]
cat("Genes detected in GSE68465:\n")
print(rownames(signature_expr))
rownames(signature_expr)[rownames(signature_expr) == "CCDC23"] <- "SVBP"
print(rownames(signature_expr))

# Transpose and merge with survival data
signature_df <- as.data.frame(t(signature_expr))
signature_df$SampleID <- rownames(signature_df)

merged_data <- merge(pheno, signature_df, by = "SampleID")
#Input: merged_data should already contain expression + survival data (SampleID, Event, Time, and all 8 genes)

# Removing samples with missing survival time
dim(merged_data)
dim(norm.tumor.exp.scaled.recentered)
na.idx <- which(is.na(merged_data$Time))
if (length(na.idx)>0) {
  merged_data <- merged_data[-na.idx,]
  norm.tumor.exp.scaled.recentered <- norm.tumor.exp.scaled.recentered[, match(rownames(merged_data), colnames(norm.tumor.exp.scaled.recentered))]
}
dim(merged_data)
################################################################################

# Define gene lists (already directionally specified)
loExpr.hiSurv <- list(thisStudy = c('GNPNAT1','XRCC2','TFAP2A','PPP1R13L'))  # higher expression = better survival
hiExpr.hiSurv <- list(thisStudy = c('ATP6V0E1','SVBP','HSDL1','UBTD1'))  

# Subset gene expression
numerator.genes <- hiExpr.hiSurv[["thisStudy"]]
denominator.genes <- loExpr.hiSurv[["thisStudy"]]

numerator.expr <- merged_data[, numerator.genes]
denominator.expr <- merged_data[, denominator.genes]

# Compute signature ratio
merged_data$ratio_thisStudy <- rowMeans(numerator.expr) / rowMeans(denominator.expr)

# Create binary 18-month OS endpoint (numericMeta from previous work was in months)
merged_data$OS_18month <- NA
merged_data$OS_18month[merged_data$Time >= (18*30.44)] <- "Survived"
merged_data$OS_18month[merged_data$Event == 1 & merged_data$Time < (18*30.44)] <- "Died"

# Drop NA
valid_idx <- which(!is.na(merged_data$OS_18month))
Survived <- ifelse(merged_data$OS_18month[valid_idx] == "Survived", 1, 0)
predictor <- merged_data$ratio_thisStudy[valid_idx]

# Fit logistic model
glm.fit <- glm(Survived ~ predictor, family = binomial)

# Plot ROC
roc_obj <- roc(Survived, glm.fit$fitted.values,
               plot=TRUE, legacy.axes=TRUE, percent=TRUE,
               xlab="False Positive Percentage", ylab="True Positive Percentage",
               col="magenta", lwd=2, main="ROC: 8-Gene Signature Predicting 18-mo OS")

# AUC and 95% CI
auc_val <- auc(roc_obj)
ci_val <- ci(roc_obj)

# Accuracy, sensitivity, specificity
report <- reportROC(gold = Survived, predictor = glm.fit$fitted.values, plot = FALSE)
report <- report[report$AUC == max(report$AUC), ]  # Select correct orientation (Died vs Survived)

# Output
cat("\n--- Performance Summary ---\n")
cat("AUC:", round(auc_val, 4), "\n")
cat("95% CI:", paste0(round(ci_val[1:2], 4), collapse = " - "), "\n")
cat("Accuracy:", report$ACC, "\n")
cat("Sensitivity:", report$SEN, "\n")
cat("Specificity:", report$SPE, "\n")

################################################################################

# # Calculate Prognostic/predictor ratios based on signature gene lists from 4 studies 
# 
# #Current LUAD Study ROC Optimization 8 gene list:  ATP6V0E1+SVBP+HSDL1+UBTD1 / GNPNAT1+XRCC2+TFAP2A+PPP1R13L
# loExpr.hiSurv[["thisStudy"]]<-c('GNPNAT1','XRCC2','TFAP2A','PPP1R13L')
# hiExpr.hiSurv[["thisStudy"]]<-c('ATP6V0E1','SVBP','HSDL1','UBTD1')
# #The above 8 genes are based on AUC outcomes of ratioed combinations of genes tested in Round 3 (17 denom; 15 numerator)
# #loExpr.hiSurv <- c('hsa-miR-143-3p','hsa-miR-30d-3p','hsa-miR-99a-3p','hsa-miR-497-5p','hsa-miR-1254','hsa-miR-193a-5p','hsa-miR-141-3p','hsa-miR-141-5p','hsa-miR-7704','hsa-miR-221-5p','hsa-miR-550a-5p','TFAP2A','PPP1R13L','GNPNAT1','GRK3','CITED2','CDK1')
# #hiExpr.hiSurv <- c('hsa-miR-143-3p','hsa-miR-30c-2-3p','hsa-miR-29b-2-5p','hsa-miR-4709-3p','hsa-miR-3065-3p','hsa-miR-135b-5p','hsa-miR-664a-3p','hsa-miR-22-5p','hsa-miR-200c-3p','BEND5','HSDL1','ATP6V0E1','ATAT1','SVBP','GADD45GIP1')
# #note: Xena TCGA portal does not have SVBP curated RNA-Seq norm_counts.
# #Paste the following signature into https://xenabrowser.net/heatmap/
# # =ATP6V0E1+SVBP+HSDL1+UBTD1-GNPNAT1-XRCC2-TFAP2A-PPP1R13L
# ################################################################################
# 
# # Confirm gene ordering (per your correction)
# hiExpr.hiSurv <- c('ATP6V0E1','SVBP','HSDL1','UBTD1')      # ↑ expression = ↓ survival
# loExpr.hiSurv <- c('GNPNAT1','XRCC2','TFAP2A','PPP1R13L')  # ↑ expression = ↑ survival
# 
# # Ensure all genes exist in data
# all_genes <- c(hiExpr.hiSurv, loExpr.hiSurv)
# stopifnot(all(all_genes %in% colnames(merged_data)))
# 
# # Compute Ratio Signature Score
# dat <- merged_data %>%
#   mutate(
#     hi_mean = rowMeans(select(., all_of(hiExpr.hiSurv)), na.rm = TRUE),
#     lo_mean = rowMeans(select(., all_of(loExpr.hiSurv)), na.rm = TRUE),
#     SignatureRatio = hi_mean / lo_mean
#   )
# 
# # Dichotomize signature (median cutoff)
# median_cutoff <- median(dat$SignatureRatio, na.rm = TRUE)
# dat <- dat %>%
#   mutate(SignatureGroup = ifelse(SignatureRatio > median_cutoff, "High", "Low"))
# 
# # Kaplan-Meier Survival Analysis
# km_fit <- survfit(Surv(Time, Event) ~ SignatureGroup, data = merged_data)