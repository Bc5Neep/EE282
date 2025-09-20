---
title: "BS_EE282_Microarray Project Script F24"
author: "Bobby Sims"
date: "`r Sys.Date()`"
output: html_document
---
## Load necessary libraries

```library(dplyr)
library(tibble)
library(affy)
library(GEOquery)
library(mouse4302.db)
library(AnnotationDbi)
library(writexl)
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(FactoMineR)
library(ggrepel)
library(msigdbr)
library(GSEABase)
library(enrichplot)
```

## Retrieve Dataset from GEO omnibus

```GetGEOSuppFiles("GSE164460")```

### Unpack files ("untar")

```untar("GSE164460/GSE164460_RAW.tar", exdir = 'data'/)```

### .tar file now unpacked to .CEL.gz files

### Further unpack into Affybatch object

```raw.data <- ReadAffy(celfile.path = 'data/')```

### Inspect raw.data stats

```raw.data```

### Data normalization

```normalized.data <- rma(raw.data)```

### Obtain normalized values

```normalized.expression <- as.data.frame( expr(normalized.data))```

### Map gene IDs to gene symbols

```gse <- getGEO("GSE164460", GSEMatrix = TRUE)```

### Obtain feature data

```feature.data <-gse$GSE164460_series_matrix.txtgz@featureData@data```

### Subset

```feature.data <- feature.data[,c(1,11)]```

```normalized.expression <- normalized.expression %>%```

```rownames_to_column (var = 'ID') %>%```

```inner_join(feature.data, by = 'ID')```

### Move Gene Symbol to first column

```normalized.expression <- normalized.expression[, c(22, setdiff(1:ncol(normalized.expression), 22))]```

### Remove ID column #

```normalized.expression$ID <- NULL```

### Create sample_names vector from column names

```sample_names <- colnames(normalized.expression [2:21])```

## Create groups

```Groups <- c( rep("TregN_Sp", 3), rep("TconvN_Sp", 3), rep ("TregEAE_Sp", 4), rep("TconvEAE_Sp", 3), rep ("TregEAE_CNS", 3), rep ("TconvEAE_CNS", 4))```

```colData <- data.frame(sample= sample_names, group=factor(Groups))```
```rownames(colData) <-sample_names```

### Build design matrix

```design <- model.matrix(~0 + group, data = colData)```

```colnames(design) <- levels(colData$group)```

### Make Contrast Matrix

```contrast_matrix <- makeContrasts(
  TregN_Sp_vs_TconvN_Sp      =  TregN_Sp - TconvN_Sp ,                  
  TregEAE_Sp_vs_TregN_Sp    =  TregEAE_Sp - TregN_Sp,                
  TregEAE_CNS_vs_TregN_Sp    = TregEAE_CNS - TregN_Sp,            
  TregEAE_CNS_vs_TregEAE_Sp  = TregEAE_CNS - TregN_Sp,         
  TconvEAE_Sp_vs_TconvN_Sp = TconvEAE_Sp - TconvN_Sp,   
  TconvEAE_CNS_vs_TconvN_Sp = TconvEAE_CNS - TconvN_Sp,
  TconvEAE_CNS_vs_TconvEAE_Sp = TconvEAE_CNS - TconvEAE_Sp,
  levels=colnames(design)
)
```

### Build expression matrix

```expr_matrix <- as.matrix(normalized.expression[,2:21])```

### Apply contrasts and empirical Bayes

```fit <- lmFit(normalized.expression, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2) 
```

### Group 1

### Extract results to toptable

```res_TregN_Sp_vs_TconvN_Sp<- topTable(fit2, coef = "TregN_Sp_vs_TconvN_Sp", number = Inf)```

```rownames(res_TregN_Sp_vs_TconvN_Sp) <- res_TregN_Sp_vs_TconvN_Sp[["Gene Symbol"]]```

### Preview Top table

```head(res_TregN_Sp_vs_TconvN_Sp)```

### Write to .csv

```write.csv(res_TregN_Sp_vs_TconvN_Sp, "limma_DEG_TregN_Sp_vs_TconvN_Sp.csv")```

### Set res0

```res0 <- res_TregN_Sp_vs_TconvN_Sp```

## Prepare differential expression data output

### 1 Check column names and classes

```print(names(res0))
print(sapply(res0[, intersect(names(res0), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res0[, intersect(names(res0), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))
```

### 2 Auto-detect sensible x (fold change) and y (adj p) columns

```xcol <- if("logFC" %in% names(res0)) "logFC" else if("log2FoldChange" %in% names(res0)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")```

```ycol <- if("adj.P.Val" %in% names(res0)) "adj.P.Val" else if("padj" %in% names(res0)) "padj" else if("P.Value" %in% names(res0)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")```

```cat("Using x =", xcol, " and y =", ycol, "\n")```

### 3 Coerce x and y to numeric safely (remove commas/extra chars)

```
safe_num <- function(x){
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)     
 
x2 <- gsub("\\(|\\)", "", x2)
  x2 <- trimws(x2)
  as.numeric(x2)
}
res0[[xcol]] <- safe_num(res0[[xcol]])
res0[[ycol]] <- safe_num(res0[[ycol]])

```

### 4 Report how many NAs resulted from coercion```cat("NAs in", xcol, ":", sum(is.na(res0[[xcol]])), "\n")

```cat("NAs in", ycol, ":", sum(is.na(res0[[ycol]])), "\n")```

```if (sum(is.na(res0[[xcol]]))>0) print(head(res0[is.na(res0[[xcol]]), ], 10))```

```cat("NAs in", xcol, ":", sum(is.na(res0[[xcol]])), "\n")```

```cat("NAs in", ycol, ":", sum(is.na(res0[[ycol]])), "\n")```

```if (sum(is.na(res0[[xcol]]))>0) print(head(res0[is.na(res0[[xcol]]), ], 10))```

### 5 Filter out rows with missing/Inf values for plotting

```keep <- !is.na(res0[[xcol]]) & is.finite(res0[[xcol]]) &```
 ```!is.na(res0[[ycol]]) & is.finite(res0[[ycol]])```

```res0_plot <- res0[keep, ]```

## Prepare volcano plot output

```sig_0 <- res_TregN_Sp_vs_TconvN_Sp[res_TregN_Sp_vs_TconvN_Sp$adj.P.Val < 0.05, ]```

```ord <- order(res_TregN_Sp_vs_TconvN_Sp$logFC, decreasing = TRUE)```

### Top 5 upregulated (highest logFC)

```top5_up_0 <- sig_0[ord[1:5], c("Gene.Symbol")]```

### Top 5 downregulated (lowest logFC)

```top5_down_0 <- res_TregN_Sp_vs_TconvN_Sp$Gene.Symbol[ord[(length(ord)-4):length(ord)]]```

### Combine

```top10_genes_0 <- c(top5_up_0, top5_down_0)```

```
top_summary_0 <- rbind(
  cbind(Direction = "Up", top5_up_0),
  cbind(Direction = "Down", top5_down_0)
)
```

### 6 Call EnhancedVolcano Plot

```
EnhancedVolcano(
  res_TregN_Sp_vs_TconvN_Sp,
  selectLab = top10_genes_0,       
  lab = res_TregN_Sp_vs_TconvN_Sp$Gene.Symbol,         
  x = 'logFC',
  y = 'adj.P.Val',
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = "TregN_Sp_vs_TconvN_Sp",
  subtitle = "Top 5 up & down genes highlighted",
  drawConnectors = TRUE,            
  widthConnectors = 0.5,            
  colConnectors = 'grey30',     
  max.overlaps = Inf)
  ```

### Group 2

### Extract results to toptable

```
res_TregEAE_Sp_vs_TregN_Sp<- topTable(fit2, coef = "TregEAE_Sp_vs_TregN_Sp", number = Inf)
```

```
rownames(res_TregEAE_Sp_vs_TregN_Sp) <- res_TregEAE_Sp_vs_TregN_Sp[["Gene Symbol"]]
```

### Preview Top table

```head(res_TregEAE_Sp_vs_TregN_Sp)```

### Write to .csv

```write.csv(res_TregEAE_Sp_vs_TregN_Sp, "limma_DEG_TregEAE_Sp_vs_TregN_Sp.csv")```

### Set res1

```res1 <- res_TregEAE_Sp_vs_TregN_Sp```


### 1 Check column names and classes

```print(names(res1))
print(sapply(res1[, intersect(names(res1), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res1[, intersect(names(res1), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))
```

### 2 Auto-detect sensible x (fold change) and y (adj p) columns

```xcol <- if("logFC" %in% names(res1)) "logFC" else if("log2FoldChange" %in% names(res1)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")```

```ycol <- if("adj.P.Val" %in% names(res1)) "adj.P.Val" else if("padj" %in% names(res1)) "padj" else if("P.Value" %in% names(res1)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")```

```cat("Using x =", xcol, " and y =", ycol, "\n")```

### 3 Coerce x and y to numeric safely (remove commas/extra chars)

```
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res1[[xcol]] <- safe_num(res1[[xcol]])
res1[[ycol]] <- safe_num(res1[[ycol]])
```

### 4 Report how many NAs resulted from coercion

```
cat("NAs in", xcol, ":", sum(is.na(res1[[xcol]])), "\n")
cat("NAs in", ycol, ":", sum(is.na(res1[[ycol]])), "\n")
if (sum(is.na(res1[[xcol]]))>0) print(head(res1[is.na(res1[[xcol]]), ], 10))
```

### 5 Filter out rows with missing/Inf values for plotting

```keep <- !is.na(res1[[xcol]]) & is.finite(res1[[xcol]]) & !is.na(res1[[ycol]]) & is.finite(res1[[ycol]])```

```res1_plot <- res1[keep, ]```


### Prepare volcano plot output

```sig_1 <- res_TregEAE_Sp_vs_TregN_Sp[res_TregEAE_Sp_vs_TregN_Sp$adj.P.Val < 0.05, ]```

```ord_1 <- order(res_TregEAE_Sp_vs_TregN_Sp$logFC, decreasing = TRUE)```

### Top 5 upregulated (highest logFC)

```top5_up_1 <- sig_1[ord_1[1:5], c("Gene.Symbol")]```

### Top 5 downregulated (lowest logFC)

```top5_down_1 <- res_TregN_Sp_vs_TconvN_Sp$Gene.Symbol[ord_1[(length(ord_1)-4):length(ord_1)]]```

```top10_genes_1 <- c(top5_up_1, top5_down_1)```

```
top_summary_1 <- rbind(
  cbind(Direction = "Up", top5_up_1),
  cbind(Direction = "Down", top5_down_1)
)
```

### 6 Call EnhancedVolcano Plot

```
EnhancedVolcano(
  res_TregEAE_Sp_vs_TregN_Sp,
  selectLab = top10_genes_1,        
  lab = res_TregEAE_Sp_vs_TregN_Sp$Gene.Symbol,  
  x = 'logFC',
  y = 'adj.P.Val',
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = "TregEAE_Sp_vs_TregN_Sp",
  subtitle = "Top 5 up & down genes highlighted",
  
  drawConnectors = TRUE,            
  widthConnectors = 0.5,            
  colConnectors = 'grey30',     
  max.overlaps = Inf                
)
```

### Group 3

### Extract results to toptable

```res_TregEAE_CNS_vs_TregN_Sp<- topTable(fit2, coef = "TregEAE_CNS_vs_TregN_Sp", number = Inf)```

```rownames(res_TregEAE_CNS_vs_TregN_Sp) <- res_TregEAE_CNS_vs_TregN_Sp[["Gene Symbol"]]```

### Preview Top table

```head(res_TregEAE_CNS_vs_TregN_Sp)```

### Write to .csv

```write.csv(res_TregEAE_CNS_vs_TregN_Sp, "limma_DEG_TregEAE_CNS_vs_TregN_Sp.csv")```

### Set res2

```res2 <- res_TregEAE_CNS_vs_TregN_Sp```

### 1 Check column names and classes

```
print(names(res2))
print(sapply(res2[, intersect(names(res2), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res2[, intersect(names(res2), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))
```

### 2 Auto-detect sensible x (fold change) and y (adj p) columns

```xcol <- if("logFC" %in% names(res2)) "logFC" else if("log2FoldChange" %in% names(res2)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")```

```ycol <- if("adj.P.Val" %in% names(res2)) "adj.P.Val" else if("padj" %in% names(res2)) "padj" else if("P.Value" %in% names(res2)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")```

```cat("Using x =", xcol, " and y =", ycol, "\n")```

### 3 Coerce x and y to numeric safely (remove commas/extra chars)

```
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res2[[xcol]] <- safe_num(res2[[xcol]])
res2[[ycol]] <- safe_num(res2[[ycol]])
```

### 4 Report how many NAs resulted from coercion

```cat("NAs in", xcol, ":", sum(is.na(res2[[xcol]])), "\n")
cat("NAs in", ycol, ":", sum(is.na(res2[[ycol]])), "\n")
if (sum(is.na(res2[[xcol]]))>0) print(head(res2[is.na(res2[[xcol]]), ], 10))
```

### 5 Filter out rows with missing/Inf values for plotting

```keep <- !is.na(res2[[xcol]]) & is.finite(res2[[xcol]]) & !is.na(res2[[ycol]]) & is.finite(res2[[ycol]])```

```res2_plot <- res2[keep, ]```


### Prepare volcano plot output

```sig_2 <- res_TregEAE_CNS_vs_TregN_Sp[res_TregEAE_CNS_vs_TregN_Sp$adj.P.Val < 0.05, ]```


```ord_2 <- order(res_TregEAE_CNS_vs_TregN_Sp$logFC, decreasing = TRUE)```

### Top 5 upregulated (highest logFC)

```top5_up_2 <- sig_2[ord_2[1:5], c("Gene.Symbol")]```


### Top 5 downregulated (lowest logFC)

```top5_down_2 <- res_TregEAE_CNS_vs_TregN_Sp$Gene.Symbol[ord_2[(length(ord_2)-4):length(ord_2)]]```

### Combine

```top10_genes_2 <- c(top5_up_2, top5_down_2)```

```
top_summary_2 <- rbind(
  cbind(Direction = "Up", top5_up_2),
  cbind(Direction = "Down", top5_down_2)
)
```

### 6 Call EnhancedVolcano Plot

```
EnhancedVolcano(
  res_TregEAE_CNS_vs_TregN_Sp ,
  selectLab = top10_genes_2,        # Top 5 up + Top 5 down
  lab = res_TregEAE_CNS_vs_TregN_Sp$Gene.Symbol,    
  x = 'logFC',
  y = 'adj.P.Val',
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = "TregEAE_Sp_vs_TregN_Sp",
  subtitle = "Top 5 up & down genes highlighted",
  drawConnectors = TRUE,            
  widthConnectors = 0.5,            
  colConnectors = 'grey30',     
  max.overlaps = Inf                
)
```

### Group 4

### Extract results to toptable

```res_TregEAE_CNS_vs_TregEAE_Sp<- topTable(fit2, coef = "TregEAE_CNS_vs_TregEAE_Sp", number = Inf)```

```rownames(res_TregEAE_CNS_vs_TregEAE_Sp) <- res_TregEAE_CNS_vs_TregEAE_Sp[["Gene Symbol"]]```

### Preview Top table

```head(res_TregEAE_CNS_vs_TregEAE_Sp)```


### Write to .csv

```write.csv(res_TregEAE_CNS_vs_TregEAE_Sp, "limma_DEG_TregEAE_CNS_vs_TregEAE_Sp.csv")```

### Set res3

```res3 <- res_TregEAE_CNS_vs_TregEAE_Sp```


### 1 Check column names and classes

```
print(names(res3))
print(sapply(res3[, intersect(names(res3), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))

print(sapply(res3[, intersect(names(res3), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))
```


### 2 Auto-detect sensible x (fold change) and y (adj p) columns

```xcol <- if("logFC" %in% names(res3)) "logFC" else if("log2FoldChange" %in% names(res3)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")```

```ycol <- if("adj.P.Val" %in% names(res3)) "adj.P.Val" else if("padj" %in% names(res3)) "padj" else if("P.Value" %in% names(res3)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")```

```cat("Using x =", xcol, " and y =", ycol, "\n")```


### 3 Coerce x and y to numeric safely (remove commas/extra chars)

```
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res3[[xcol]] <- safe_num(res3[[xcol]])
res3[[ycol]] <- safe_num(res3[[ycol]])
```

### 4 Report how many NAs resulted from coercion

```cat("NAs in", xcol, ":", sum(is.na(res3[[xcol]])), "\n")```

```cat("NAs in", ycol, ":", sum(is.na(res3[[ycol]])), "\n")```

```if (sum(is.na(res3[[xcol]]))>0) print(head(res3[is.na(res3[[xcol]]), ], 10))```

### 5 Filter out rows with missing/Inf values for plotting


```keep <- !is.na(res3[[xcol]]) & is.finite(res3[[xcol]]) & !is.na(res3[[ycol]]) & is.finite(res3[[ycol]])```

```res3_plot <- res3[keep, ]```

### Prepare volcano plot output

```sig_3 <- res_TregEAE_CNS_vs_TregEAE_Sp[res_TregEAE_CNS_vs_TregEAE_Sp$adj.P.Val < 0.05, ]```

```ord_3 <- order(res_TregEAE_CNS_vs_TregEAE_Sp$logFC, decreasing = TRUE)```


### Top 5 upregulated (highest logFC)

```top5_up_3 <- sig_3[ord_3[1:5], c("Gene.Symbol")]```


### Top 5 downregulated (lowest logFC)

```top5_down_3 <- res_TregEAE_CNS_vs_TregEAE_Sp$Gene.Symbol[ord_3[(length(ord_3)-4):length(ord_3)]]```


### Combine

```top10_genes_3 <- c(top5_up_3, top5_down_3)```

```  
top_summary_3 <- rbind(
  cbind(Direction = "Up", top5_up_3),
  cbind(Direction = "Down", top5_down_3))
  ```

### 6 Call EnhancedVolcano Plot

```
EnhancedVolcano(
  res_TregEAE_CNS_vs_TregEAE_Sp ,
  selectLab = top10_genes_3,       
  lab = res_TregEAE_CNS_vs_TregEAE_Sp$Gene.Symbol),
  x = 'logFC',
  y = 'adj.P.Val',
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = "TregEAE_Sp_vs_TregEAE_Sp",
  subtitle = "Top 5 up & down genes highlighted",
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey30',
  max.overlaps = Inf
  )
  ```


### Group 5

### Extract results to toptable

```res_TconvEAE_CNS_vs_TconvN_Sp <- topTable(fit2, coef = "TconvEAE_CNS_vs_TconvN_Sp", number = Inf)```

```rownames(res_TconvEAE_CNS_vs_TconvN_Sp) <- res_TconvEAE_CNS_vs_TconvN_Sp[["Gene Symbol"]]```

### Preview Top table

```head(res_TconvEAE_CNS_vs_TconvN_Sp)```

### Write to .csv

```write.csv(res_TconvEAE_CNS_vs_TconvN_Sp, "limma_DEG_TconvEAE_CNS_vs_TconvN_Sp.csv")```


### Set res4

```res4 <- res_TconvEAE_CNS_vs_TconvN_Sp```

### 1 Check column names and classes

```print(names(res4))
print(sapply(res4[, intersect(names(res4), c("logFC","log2FoldChange","log2FoldChange","logFC")) , drop=FALSE], class))
print(sapply(res4[, intersect(names(res4), c("adj.P.Val","padj","P.Value")), drop=FALSE], class))
```

### 2 Auto-detect sensible x (fold change) and y (adj p) columns

```
xcol <- if("logFC" %in% names(res4)) "logFC" else if("log2FoldChange" %in% names(res4)) "log2FoldChange" else stop("No fold-change column found (logFC or log2FoldChange).")
```

```
ycol <- if("adj.P.Val" %in% names(res4)) "adj.P.Val" else if("padj" %in% names(res4)) "padj" else if("P.Value" %in% names(res4)) "P.Value" else stop("No p-value column found (adj.P.Val, padj or P.Value).")
```
```cat("Using x =", xcol, " and y =", ycol, "\n")```

### 3 Coerce x and y to numeric safely (remove commas/extra chars)

```
safe_num <- function(x) {
  x2 <- as.character(x)
  x2 <- gsub(",", "", x2)         # remove thousands separators
  x2 <- gsub("\\(|\\)", "", x2)   # remove parentheses if any
  x2 <- trimws(x2)
  as.numeric(x2)
}
res4[[xcol]] <- safe_num(res4[[xcol]])
res4[[ycol]] <- safe_num(res4[[ycol]])
```

### 4 Report how many NAs resulted from coercion

```
cat("NAs in", xcol, ":", sum(is.na(res4[[xcol]])), "\n")
cat("NAs in", ycol, ":", sum(is.na(res4[[ycol]])), "\n")
if (sum(is.na(res4[[xcol]]))>0) print(head(res4[is.na(res4[[xcol]]), ], 10))
```

### 5 Filter out rows with missing/Inf values for plotting

```keep <- !is.na(res4[[xcol]]) & is.finite(res4[[xcol]]) & !is.na(res4[[ycol]]) & is.finite(res4[[ycol]])```

```res4_plot <- res4[keep, ]```


### Prepare volcano plot output

```sig_4 <- res_TconvEAE_CNS_vs_TconvN_Sp[res_TconvEAE_CNS_vs_TconvN_Sp$adj.P.Val < 0.05, ]```

```ord_4 <- order(res_TconvEAE_CNS_vs_TconvN_Sp$logFC, decreasing = TRUE)```


### Top 5 upregulated (highest logFC)

```top5_up_4 <- sig_4[ord_3[1:5], c("Gene.Symbol")]```


### Top 5 downregulated (lowest logFC)

```
top5_down_4 <- res_TregEAE_CNS_vs_TregEAE_Sp$Gene.Symbol[ord_4[(length(ord_4)-4):length(ord_4)]]
```
### Combine

```top10_genes_4 <- c(top5_up_4, top5_down_4)```

```
top_summary_4 <- rbind(
  cbind(Direction = "Up", top5_up_4),
  cbind(Direction = "Down", top5_down_4)
)
```

### 6 Call EnhancedVolcano Plot

```
EnhancedVolcano(
  res_TconvEAE_CNS_vs_TconvN_Sp ,
  selectLab = top10_genes_4,        # Top 5 up + Top 5 down
  lab = res_TconvEAE_CNS_vs_TconvN_Sp$Gene.Symbol,    
  x = 'logFC',
  y = 'adj.P.Val',
  FCcutoff = 1.5,
  pCutoff = 0.05,
  col = c('black','green','blue','red'),
  title = "TconvEAE_CNS_vs_TconvN_Sp",
  subtitle = "Top 5 up & down genes highlighted",
  drawConnectors = TRUE,            
  widthConnectors = 0.5,            
  colConnectors = 'grey30',     
  max.overlaps = Inf                
)
```

# MA plot

### Group 1

### Filtered data for the top 10 genes

```highlight <- res_TregN_Sp_vs_TconvN_Sp[res_TregN_Sp_vs_TconvN_Sp$Gene.Symbol %in% top10_genes_0, ]```

### MA plot

```
ggplot(res_TregN_Sp_vs_TconvN_Sp, aes(x = AveExpr, y = logFC)) +
  geom_point(color = "grey", size = 1) +
  geom_point(data = highlight, aes(x = AveExpr, y = logFC), color = "red", size = 2) +
  geom_text_repel(data = highlight, aes(x = AveExpr, y = logFC, label = Gene.Symbol),
 size = 3, color = "blue") +
  labs(x = "Average Expression (A)",
 y = "Log2 Fold Change (M)",
title = "TregN_Sp_vs_TconvN_Sp",
subtitle = "Top 5 up & downregulated genes highlighted",
theme_minimal()
  )
  ```

### Group 2

### Filtered data for the top 10 genes

```highlight <- res_TregEAE_Sp_vs_TregN_Sp[res_TregEAE_Sp_vs_TregN_Sp$Gene.Symbol %in% top10_genes_1, ]```

### MA plot

```
ggplot(res_TregEAE_Sp_vs_TregN_Sp, aes(x = AveExpr, y = logFC)) +
  geom_point(color = "grey", size = 1) +
  geom_point(data = highlight, aes(x = AveExpr, y = logFC), color = "red", size = 2) +
  geom_text_repel(data = highlight, aes(x = AveExpr, y = logFC, label = Gene.Symbol),
  size = 3, color = "blue") +
  labs(x = "Average Expression (A)",
  y = "Log2 Fold Change (M)",
  title="TregEAE_Sp_vs_TregN_Sp",
  subtitle = "Top 5 up & downregulated genes highlighted",
  theme_minimal()
  )
  ```

### Group 3

### Filtered data for the top 10 genes

```highlight <- res_TregEAE_CNS_vs_TregN_Sp[res_TregEAE_CNS_vs_TregN_Sp$Gene.Symbol %in% top10_genes_2,]```

### MA plot

```ggplot(res_TregEAE_CNS_vs_TregN_Sp, aes(x = AveExpr, y = logFC)) +
  geom_point(color = "grey", size = 1) +
  geom_point(data = highlight, aes(x = AveExpr, y = logFC), color = "red", size = 2) +
  geom_text_repel(data = highlight, aes(x = AveExpr, y = logFC, label = Gene.Symbol),
 size = 3, color = "blue") +
  labs(x = "Average Expression (A)",
  y = "Log2 Fold Change (M)",
  title = "TregEAE_CNS_vs_TregN_Sp",
  subtitle = "Top 5 up & downregulated genes highlighted",
  theme_minimal()
  )
  ```

### Group 4

### Filtered data for the top 10 genes

```highlight <- res_TregEAE_CNS_vs_TregEAE_Sp[res_TregEAE_CNS_vs_TregEAE_Sp$Gene.Symbol %in% top10_genes_2, ]```

### MA plot

```ggplot(res_TregEAE_CNS_vs_TregN_Sp, aes(x = AveExpr, y = logFC)) +
  geom_point(color = "grey", size = 1) +
  geom_point(data = highlight, aes(x = AveExpr, y = logFC), color = "red", size = 2) +
  geom_text_repel(data = highlight, aes(x = AveExpr, y = logFC, label = Gene.Symbol),
  size = 3, color = "blue") +
  labs(x = "Average Expression (A)",
  y = "Log2 Fold Change (M)",
  title = "TregEAE_CNS_vs_TregN_Sp",
  subtitle = "Top 5 up & downregulated genes highlighted",
  theme_minimal()
  )
  ```

### Group 5

### Filtered data for the top 10 genes

```highlight <- res_TconvEAE_CNS_vs_TconvN_Sp[res_TconvEAE_CNS_vs_TconvN_Sp$Gene.Symbol %in% top10_genes_4, ]```

### MA plot

```ggplot(res_TconvEAE_CNS_vs_TconvN_Sp, aes(x = AveExpr, y = logFC)) +
  geom_point(color = "grey", size = 1) +
  geom_point(data = highlight, aes(x = AveExpr, y = logFC), color = "red", size = 2) +
  geom_text_repel(data = highlight, aes(x = AveExpr, y = logFC, label = Gene.Symbol),
  size = 3, color = "blue") +
  labs(x = "Average Expression (A)",
  y = "Log2 Fold Change (M)",
 title = "TconvEAE_CNS_vs_TconvN_Sp",
 subtitle = "Top 5 up & downregulated genes highlighted",
theme_minimal()
  )
  ```
  
# Generate PCA Plot

### Expression matrix (genes in rows, samples in columns)

```Expr_Matrix <- as.matrix(normalized.expression[,-1])```

```pca <- prcomp(t(Expr_Matrix), scale. = TRUE)```

### Extract PCA scores

```pca_data <- as.data.frame(pca$x)```


### Add sample info (from your colData)

```pca_data$sample <- rownames(pca_data)```

```pca_data$group <- colData$group```


## PCA Plot

```ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = group)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCA of CD4 DEGs: Control vs EAE",
  x = paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% variance)"),
  y = paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% variance)")) +
  theme_minimal()
  ```

# GSEA

## GSEA group#1

### Get all contrast names

```contrast_names <- colnames(fit2)```


### Create a list to store results

```topTables <- list()```

### Loop through contrasts and extract DEG tables

```
for (contrast in contrast_names) {
  topTables[[contrast]] <- topTable(
    fit2,
    coef = contrast,
    number = Inf,
    sort.by = "P"   # sort by adjusted p-values
  )
}
```

### Source contrast

```res_TregEAE_CNS_vs_TregN_Sp<- topTable(fit2, coef = "TregEAE_CNS_vs_TregN_Sp", number = Inf)```

### Clear of any NA values

```res_clean1 <- res_TregEAE_CNS_vs_TregN_Sp [!is.na(res_TregEAE_CNS_vs_TregN_Sp $Gene.Symbol) & res_TregEAE_CNS_vs_TregN_Sp $Gene.Symbol != "", ]```


### Fix for duplicates

```res_dedup <- res_clean[!is.na(res_clean$Gene.Symbol) & res_clean$Gene.Symbol != "", ]```


### For genes with multiple probes, keep the one with the highest absolute logFC

```res_dedup <- res_dedup[order(abs(res_dedup$logFC), decreasing = TRUE), ]```

```res_dedup <- res_dedup[!duplicated(res_dedup$Gene.Symbol), ]```


### Build geneList

```geneList <- res_dedup$logFC```

```names(geneList) <- res_dedup$Gene.Symbol```

### Sort

```geneList <- sort(geneList, decreasing = TRUE)```


### Check for dupes

```any(duplicated(names(geneList)))```

## Run GSEA

```
gsea_ir <- GSEA(
  geneList    = geneList,
  TERM2GENE   = ir_term2gene,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
```

### View results
```head(gsea_ir@result)```

### Plot

```gseaplot2(gsea_ir, geneSetID = 1)```

## GSEA group#2

### Source contrast

```res_TconvEAE_CNS_vs_TconvN_Sp <- topTable(fit2, coef = "TconvEAE_CNS_vs_TconvN_Sp", number = Inf)```

### Clear of any NA values

```res_clean1 <- res_TconvEAE_CNS_vs_TconvN_Sp [!is.na(res_TconvEAE_CNS_vs_TconvN_Sp $Gene.Symbol) & res_TconvEAE_CNS_vs_TconvN_Sp $Gene.Symbol != "", ]```


### Fix for duplicates

```res_dedup1 <- res_clean1[!is.na(res_clean1$Gene.Symbol) & res_clean1$Gene.Symbol != "", ]```

### For genes with multiple probes, keep the one with the highest absolute logFC

```res_dedup1 <- res_dedup1[order(abs(res_dedup1$logFC), decreasing = TRUE), ]```

```res_dedup1 <- res_dedup1[!duplicated(res_dedup1$Gene.Symbol), ]```


### Build geneList

```geneList1 <- res_dedup1$logFC```

```names(geneList1) <- res_dedup1$Gene.Symbol```

### Sort

```geneList1 <- sort(geneList1, decreasing = TRUE)```

### Check for dupes

```any(duplicated(names(geneList1)))```

## Run GSEA

```
gsea_ir1 <- GSEA(
  geneList    = geneList1,
  TERM2GENE   = ir_term2gene,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
```

### View results

```head(gsea_ir1@result)```

### Plot

```gseaplot2(gsea_ir1, geneSetID = 1)```
