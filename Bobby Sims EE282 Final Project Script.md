## Load necessary libraries

library(dplyr)
library(tibble)
library(affy)
library(GEOquery)
library(writexl)
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(FactoMineR)

## Retrieve Dataset from GEO omnibus
GetGEOSuppFiles(“GSE164460”) 

## Unpack files (“untar”)

untar("GSE164460/GSE164460_RAW.tar", exdir = 'data'/)

## .tar file now unpacked to .CEL.gz files


##Further unpack into Affybatch object

raw.data <- ReadAffy(celfile.path = 'data/')

##Inspect raw.data stats

raw.data

## Data normalization
normalized.data <- rma(raw.data) 

## Obtain normalized values
normalized.expression <- as.data.frame( exprs(normalized.data))


## Map gene IDs to gene symbols
gse <- getGEO("GSE164460", GSEMatrix = TRUE)

## Obtain feature data 
feature.data <-gse$GSE164460_series_matrix.txtgz@featureData@data

## Subset
feature.data <- feature.data[,c(1,11)]

# Set sample names

# Create sample_names vector from column names 3 to 22
sample_names <- colnames(normalized.expression)[3:22]

colData <- data.frame(sample= sample_names, group=factor(Groups))
rownames(colData) <-sample_names

#Build design matrix
design <- model.matrix(~0 + group, data = colData)
colnames(design) <- levels(colData$group)

# Make Contrast Matrix
contrast_matrix <- makeContrasts(
  TregN_Sp_vs_TconvN_Sp      =  TregN_Sp - TconvN_Sp ,                  
  TregEAE_Sp_vs_TregN_Sp    =  TregEAE_Sp - TregN_Sp,                
  TregEAE_CNS_vs_TregN_Sp    = TregEAE_CNS - TregN_Sp,            
  TregEAE_CNS_vs_TregEAE_Sp  = TregEAE_CNS - TregN_Sp,         
  TconvEAE_Sp_vs_TconvN_Sp = TconvEAE_Sp - TconvN_Sp,   
  TconvEAE_CNS_vs_TconvN_Sp = TconvEAE_CNS - TconvN_Sp,
  TconvEAE_CNS_vs_TconvEAE_Sp = TconvEAE_CNS - TconvEAE_Sp,
  levels=colnames(design)
)   



# Apply contrasts and empirical Bayes
fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)



# Extract results to toptable
res_TregN_Sp_vs_TconvN_Sp <- topTable(fit2, coef = "TregN_Sp_vs_TconvN_Sp", number=Inf, adjust.method = "BH")

# Preview Top table
head(res_TregN_Sp_vs_TconvN_Sp)

# Write to .csv
write.csv(res_TregN_Sp_vs_TconvN_Sp, "limma_DEG_TregN_Sp_vs_TconvN_Sp.csv")


# Set res0
res0 <- res_TregN_Sp_vs_TconvN_Sp 


# 1) Check column names and classes
print(names(res0))
print(sapply(res0[, intersect(names(res0), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res0[, intersect(names(res0), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))

# 2) Auto-detect sensible x (fold change) and y (adj p) columns
xcol <- if("logFC" %in% names(res0)) "logFC" else if("log2FoldChange" %in% names(res0)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")
ycol <- if("adj.P.Val" %in% names(res0)) "adj.P.Val" else if("padj" %in% names(res0)) "padj" else if("P.Value" %in% names(res0)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")
cat("Using x =", xcol, " and y =", ycol, "\n")

# 3) Coerce x and y to numeric safely (remove commas/extra chars)
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res0[[xcol]] <- safe_num(res0[[xcol]])
res0[[ycol]] <- safe_num(res0[[ycol]])

# 4) Report how many NAs resulted from coercion
cat("NAs in", xcol, ":", sum(is.na(res0[[xcol]])), "\n")
cat("NAs in", ycol, ":", sum(is.na(res0[[ycol]])), "\n")
if (sum(is.na(res0[[xcol]]))>0) print(head(res0[is.na(res0[[xcol]]), ], 10))


# 5) Filter out rows with missing/Inf values for plotting
keep <- !is.na(res0[[xcol]]) & is.finite(res0[[xcol]]) & !is.na(res0[[ycol]]) & is.finite(res0[[ycol]])
res0_plot <- res0[keep, ]


# 6) Call EnhancedVolcano (use exact column names detected)

top20_grp0 <- rownames(res_TregN_Sp_vs_TconvN_Sp)[order(res_TregN_Sp_vs_TconvN_Sp$adj.P.Val)][1:20]



EnhancedVolcano(
  res_TregN_Sp_vs_TconvN_Sp,                      
  selectLab = top20_grp0, 
  lab = rownames(res_TregN_Sp_vs_TconvN_Sp),  
  x = 'logFC',                  
  y = 'adj.P.Val',          
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = "TregN_Sp_vs_TconvN_Sp",
  subtitle = "Significant outliers highlighted"
)

# Extract results to toptable
res_TregEAE_Sp_vs_TregN_Sp <- topTable(fit2, coef = "TregEAE_Sp_vs_TregN_Sp", number=Inf, adjust.method = "BH")


# Preview Top table
head(res_TregEAE_Sp_vs_TregN_Sp)

# Write to .csv
write.csv(res_TregEAE_Sp_vs_TregN_Sp, "limma_DEG_TregEAE_Sp_vs_TregN_Sp.csv")


# Set res1
res1 <- res_TregEAE_Sp_vs_TregN_Sp

# 1) Check column names and classes
print(names(res1))
print(sapply(res1[, intersect(names(res1), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res1[, intersect(names(res1), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))

# 2) Auto-detect sensible x (fold change) and y (adj p) columns
xcol <- if("logFC" %in% names(res1)) "logFC" else if("log2FoldChange" %in% names(res1)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")
ycol <- if("adj.P.Val" %in% names(res1)) "adj.P.Val" else if("padj" %in% names(res1)) "padj" else if("P.Value" %in% names(res1)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")
cat("Using x =", xcol, " and y =", ycol, "\n")

# 3) Coerce x and y to numeric safely (remove commas/extra chars)
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res1[[xcol]] <- safe_num(res1[[xcol]])
res1[[ycol]] <- safe_num(res1[[ycol]])

# 4) Report how many NAs resulted from coercion
cat("NAs in", xcol, ":", sum(is.na(res1[[xcol]])), "\n")
cat("NAs in", ycol, ":", sum(is.na(res1[[ycol]])), "\n")
if (sum(is.na(res1[[xcol]]))>0) print(head(res1[is.na(res1[[xcol]]), ], 10))


# 5) Filter out rows with missing/Inf values for plotting
keep <- !is.na(res1[[xcol]]) & is.finite(res1[[xcol]]) & !is.na(res1[[ycol]]) & is.finite(res1[[ycol]])
res1_plot <- res1[keep, ]


# 6) Call EnhancedVolcano (use exact column names detected

top20_grp1 <- rownames(res_TregEAE_Sp_vs_TregN_Sp)[order(res_TregEAE_Sp_vs_TregN_Sp$adj.P.Val)][1:20]

EnhancedVolcano(
  res_TregEAE_Sp_vs_TregN_Sp,                     
  selectLab = top20_grp1,   
  lab = rownames(res_TregEAE_Sp_vs_TregN_Sp),  
  x = 'logFC',                  
  y = 'adj.P.Val',          
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = "TregEAE_Sp_vs_TregN_Sp",
  subtitle = "Significant outliers highlighted"
)

# Preview Top table
head(res_TregEAE_CNS_vs_TregN_Sp)

# Write to .csv
write.csv(res_TregEAE_CNS_vs_TregN_Sp, "limma_DEG_TregEAE_CNS_vs_TregN_Sp ")


# Set res2
res2 <- res_TregEAE_CNS_vs_TregN_Sp



# 1) Check column names and classes
print(names(res2))
print(sapply(res2[, intersect(names(res2), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res2[, intersect(names(res2), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))

# 2) Auto-detect sensible x (fold change) and y (adj p) columns
xcol <- if("logFC" %in% names(res2)) "logFC" else if("log2FoldChange" %in% names(res2)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")
ycol <- if("adj.P.Val" %in% names(res2)) "adj.P.Val" else if("padj" %in% names(res2)) "padj" else if("P.Value" %in% names(res2)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")
cat("Using x =", xcol, " and y =", ycol, "\n")

# 3) Coerce x and y to numeric safely (remove commas/extra chars)
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res2[[xcol]] <- safe_num(res2[[xcol]])
res2[[ycol]] <- safe_num(res2[[ycol]])

# 4) Report how many NAs resulted from coercion
cat("NAs in", xcol, ":", sum(is.na(res2[[xcol]])), "\n")
cat("NAs in", ycol, ":", sum(is.na(res2[[ycol]])), "\n")
if (sum(is.na(res2[[xcol]]))>0) print(head(res2[is.na(res2[[xcol]]), ], 10))


# 5) Filter out rows with missing/Inf values for plotting
keep <- !is.na(res2[[xcol]]) & is.finite(res2[[xcol]]) & !is.na(res2[[ycol]]) & is.finite(res2[[ycol]])
res2_plot <- res2[keep, ]


# 6) Call EnhancedVolcano (use exact column names detected)

top20_grp2 <- rownames(res_TregEAE_CNS_vs_TregN_Sp)[order(res_TregEAE_CNS_vs_TregN_Sp$adj.P.Val)][1:20]

EnhancedVolcano(
  res_TregEAE_CNS_vs_TregN_Sp,                     
  selectLab = top20_grp2,  
  lab = rownames(res_TregEAE_CNS_vs_TregN_Sp),  
  x = 'logFC',                  
  y = 'adj.P.Val',          
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = " TregEAE_CNS_vs_TregN_Sp ",
  subtitle = "Significant outliers highlighted"
)

# Preview Top table
head(res_TregEAE_CNS_vs_TregEAE_Sp)

#Write to .csv
write.csv(res_TregEAE_CNS_vs_TregEAE_Sp, "limma_DEG_TregEAE_CNS_vs_TregEAE_Sp")


# Set res3
res3 <- res_TregEAE_CNS_vs_TregEAE_Sp



# 1) Check column names and classes
print(names(res3))
print(sapply(res3[, intersect(names(res3), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res3[, intersect(names(res3), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))

# 2) Auto-detect sensible x (fold change) and y (adj p) columns
xcol <- if("logFC" %in% names(res3)) "logFC" else if("log2FoldChange" %in% names(res3)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")
ycol <- if("adj.P.Val" %in% names(res3)) "adj.P.Val" else if("padj" %in% names(res3)) "padj" else if("P.Value" %in% names(res3)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")
cat("Using x =", xcol, " and y =", ycol, "\n")

# 3) Coerce x and y to numeric safely (remove commas/extra chars)
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res3[[xcol]] <- safe_num(res3[[xcol]])
res3[[ycol]] <- safe_num(res3[[ycol]])

# 4) Report how many NAs resulted from coercion
cat("NAs in", xcol, ":", sum(is.na(res3[[xcol]])), "\n")
cat("NAs in", ycol, ":", sum(is.na(res3[[ycol]])), "\n")
if (sum(is.na(res3[[xcol]]))>0) print(head(res3[is.na(res3[[xcol]]), ], 10))


# 5) Filter out rows with missing/Inf values for plotting
keep <- !is.na(res3[[xcol]]) & is.finite(res3[[xcol]]) & !is.na(res3[[ycol]]) & is.finite(res3[[ycol]])
res3_plot <- res3[keep, ]


# 6) Call EnhancedVolcano (use exact column names detected)

top20_grp3 <- rownames(res_TregEAE_CNS_vs_TregEAE_Sp)[order(res_TregEAE_CNS_vs_TregEAE_Sp$adj.P.Val)][1:20]

EnhancedVolcano(
  res_TregEAE_CNS_vs_TregEAE_Sp,                     
  selectLab = top20_grp3,  
  lab = rownames(res_TregEAE_CNS_vs_TregEAE_Sp),  
  x = 'logFC',                  
  y = 'adj.P.Val',          
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = " TregEAE_CNS_vs_TregEAE_Sp ",
  subtitle = "Significant outliers highlighted"
)



## Generate PCA Plot
# Expression matrix (genes in rows, samples in columns)
expr_matrix <- as.matrix(normalized.expression[, 3:22])   # assuming cols 3–22 are samples

pca <- prcomp(t(expr_matrix), scale. = TRUE)

# Extract PCA scores
pca_data <- as.data.frame(pca$x)

# Add sample info (from your colData)
pca_data$sample <- rownames(pca_data)
pca_data$group  <- colData$group   # assuming you built colData earlier

library(ggplot2)

ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCA of Treg | Tconv DEGs",
       x = paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% variance)"),
       y = paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% variance)")) +
  theme_minimal()
