install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "edgeR",
  "limma",
  "pheatmap",
  "EnhancedVolcano"
))
library(GEOquery)
library(edgeR)
library(limma)
library(pheatmap)
library(EnhancedVolcano)

gset <- getGEO("GSEXXXX", GSEMatrix = TRUE)

gset <- getGEO("GSE96058", GSEMatrix = TRUE)
length(gset)
gset1 <- gset[[1]]
expr <- exprs(gset1)
pheno <- pData(gset1)
head(pheno)

expr <- exprs(gset[[1]])
pheno <- pData(gset[[1]])

table(pheno$`er status`)
table(pheno$`pam50 subtype`)

group <- factor(pheno$`pam50 subtype`)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr, design)
contrast <- makeContrasts(Basal - LumA, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTable(fit2)

EnhancedVolcano(fit2, lab=rownames(expr))

pheatmap(expr[1:50, ])
table(pheno$`pam50 subtype`)

group <- pheno$`pam50 subtype`

keep <- group %in% c("Basal", "LumA")
expr2 <- expr[, keep]
group2 <- factor(group[keep])

class(expr)
str(expr[1:5,1:5])

expr <- as.matrix(exprs(gset[[1]]))
pheno <- pData(gset[[1]])

group <- pheno$`pam50 subtype`

keep <- group %in% c("Basal", "LumA")

expr2 <- expr[, keep]
group2 <- factor(group[keep])

class(expr)
dim(expr)
expr <- exprs(gset[[1]])
expr <- as.matrix(expr)
dim(expr)
expr <- t(expr)
pheatmap(expr[1:50, 1:10])
expr <- exprs(gset[[1]])

dim(expr)
class(expr)
expr <- as.matrix(expr)

dim(expr)
expr[1:5, 1:5]
if (nrow(expr) < ncol(expr)) {
  expr <- t(expr)
}
n_genes <- min(50, nrow(expr))
n_samples <- min(10, ncol(expr))

expr_small <- expr[1:n_genes, 1:n_samples]
pheatmap(expr[seq_len(min(50, nrow(expr))), seq_len(min(10, ncol(expr)))])

pheno <- pData(gset[[1]])
group <- pheno$`pam50 subtype`

table(group)
rm(list = ls())
library(GEOquery)
library(limma)
library(pheatmap)
library(EnhancedVolcano)

gset <- getGEO("GSE42568", GSEMatrix = TRUE)
eset <- gset[[1]]

expr <- exprs(eset)
pheno <- pData(eset)
dim(expr)
colnames(pheno)
table(pheno$`tissue:ch1`)
group <- factor(pheno$`tissue:ch1`)

table(group)
design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

design
fit <- lmFit(expr, design)
group <- factor(pheno$`tissue:ch1`)

levels(group)
levels(group) <- c("breast_cancer", "normal_breast")
design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

colnames(design)
contrast.matrix <- makeContrasts(
  breast_cancer - normal_breast,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf)

head(results)
EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1
)
topgenes <- rownames(results)[1:50]

pheatmap(
  expr[topgenes, ],
  scale = "row"
)
sig <- results[results$adj.P.Val < 0.05, ]

dim(sig)
head(sig)
top_up <- head(sig[order(sig$logFC, decreasing = TRUE), ], 10)
top_down <- head(sig[order(sig$logFC), ], 10)

top_up
top_down
write.csv(results, "DEG_results_breast_cancer.csv")
write.csv(sig, "significant_genes.csv")
pca <- prcomp(t(expr), scale. = TRUE)

plot(pca$x[,1], pca$x[,2],
     col = as.factor(group),
     pch = 19,
     xlab = "PC1",
     ylab = "PC2")
legend("topright", legend = levels(group),
       col = 1:length(levels(group)), pch = 19)
expr_var <- apply(expr, 1, var)

expr_filtered <- expr[expr_var > 0, ]
pca <- prcomp(t(expr_filtered), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     col = as.factor(group),
     pch = 19,
     xlab = "PC1",
     ylab = "PC2")

legend("topright",
       legend = levels(group),
       col = 1:length(levels(group)),
       pch = 19)
sig <- results[results$adj.P.Val < 0.05, ]

dim(sig)
head(sig)
up_genes <- sig[sig$logFC > 1, ]
down_genes <- sig[sig$logFC < -1, ]

nrow(up_genes)
nrow(down_genes)
top_up <- head(up_genes[order(up_genes$logFC, decreasing = TRUE), ], 10)
top_down <- head(down_genes[order(down_genes$logFC), ], 10)

top_up
top_down
genes_of_interest <- c("ESR1", "ERBB2", "PGR", "MKI67")

results[genes_of_interest, ]
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
deg <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
head(deg)
top_genes <- rownames(deg)[1:min(50, nrow(deg))]
pheatmap(expr[top_genes, , drop = FALSE])
expr_clean <- expr[apply(expr, 1, var) > 0, ]

pca <- prcomp(t(expr_clean), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], col = as.factor(group))
group <- pheno$`pam50 subtype`
group <- as.factor(group)
table(group)
head(deg[order(deg$logFC, decreasing = TRUE), ])
head(deg[order(deg$logFC), ])

#How do Basal vs LumA tumors differ at the transcriptomic level?

table(group)
group <- factor(group)
library(matrixStats)

subtype_means <- sapply(levels(group), function(g) {
  rowMeans(expr[, group == g, drop = FALSE])
})
basal_luma <- results[rownames(results) %in% rownames(deg), ]
up_in_basal <- rownames(deg[deg$logFC > 1, ])
up_in_luma  <- rownames(deg[deg$logFC < -1, ])
sig_genes <- unique(c(up_in_basal[1:20], up_in_luma[1:20]))

pheatmap(expr[sig_genes, , drop = FALSE],
         annotation_col = data.frame(Subtype = group))
length(sig_genes)
sum(sig_genes %in% rownames(expr))
sig_genes <- unique(c(up_in_basal[1:20], up_in_luma[1:20]))

sig_genes <- sig_genes[sig_genes %in% rownames(expr)]
length(sig_genes)
mat <- expr[sig_genes, , drop = FALSE]
dim(mat)
pheatmap(mat,
         annotation_col = data.frame(Subtype = group))
sig_genes
length(sig_genes)
sum(sig_genes %in% rownames(expr))
dim(expr)
sig_genes <- rownames(results)[1:40]
mat <- expr[sig_genes, , drop = FALSE]
dim(mat)
pheatmap(mat, annotation_col = data.frame(Subtype = group))
dim(mat)
head(mat)
sum(is.na(mat))
# remove low-variance genes
vars <- apply(mat, 1, var)
mat_filt <- mat[vars > 0, ]
dim(mat_filt)
summary(vars)
pca <- prcomp(t(mat_filt), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     col = as.factor(group),
     pch = 19)
library(pheatmap)

top_genes <- head(order(apply(mat_filt, 1, var), decreasing = TRUE), 30)

pheatmap(mat_filt[top_genes, ],
         scale = "row")
dim(expr)
library(EnhancedVolcano)

EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Basal vs LumA Differential Expression",
  subtitle = "Breast Cancer (GEO dataset)",
  xlab = "Log2 Fold Change",
  ylab = "-log10 Adjusted P-value"
)
library(pheatmap)

# select top 30 significant genes
top_genes <- rownames(results)[1:30]

# extract expression safely
heat_mat <- expr[top_genes, ]

# handle possible duplicates or missing genes
heat_mat <- heat_mat[complete.cases(heat_mat), ]

pheatmap(
  heat_mat,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Top DE Genes Heatmap (Basal vs LumA)"
)
# remove zero-variance genes
vars <- apply(expr, 1, var)
expr_filt <- expr[vars > 0, ]

# PCA
pca <- prcomp(t(expr_filt), scale. = TRUE)

# plot
group <- pheno$`pam50 subtype`

plot(pca$x[,1], pca$x[,2],
     col = as.factor(group),
     pch = 19,
     main = "PCA: Breast Cancer Subtypes")

legend("topright",
       legend = unique(group),
       col = 1:length(unique(group)),
       pch = 19)
colnames(pheno)
table(pheno$`pam50 subtype`)
group <- pheno[, grep("pam50|subtype", colnames(pheno), ignore.case = TRUE)]

head(group)
pheno$characteristics_ch1[1:10]
pheno$description[1:10]
group <- pheno$characteristics_ch1

group <- sub(".*pam50 subtype: ", "", group)
group <- sub(" .*", "", group)   # remove anything after subtype

table(group)
pca <- prcomp(t(expr_filt), scale. = TRUE)

plot(pca$x[,1], pca$x[,2],
     col = as.factor(group),
     pch = 19,
     main = "PCA: Breast Cancer Subtypes")

legend("topright",
       legend = unique(group),
       col = 1:length(unique(group)),
       pch = 19)
pca <- prcomp(t(expr_filt), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     col = as.factor(group),
     pch = 19,
     main = "PCA: Breast Cancer Subtypes",
     xlab = "PC1",
     ylab = "PC2")
legend("topright",
       legend = levels(as.factor(group)),
       col = 1:length(unique(group)),
       pch = 19)

#Transcriptomic profiling of breast cancer (GSE dataset) identified distinct molecular clustering between PAM50 subtypes and revealed differentially expressed genes between Basal and Luminal A tumors. 

library(EnhancedVolcano)

EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Basal vs LumA Differential Expression",
  subtitle = "PAM50 subtype comparison",
  caption = "GSE42568 analysis"
)
top_genes <- results[order(results$adj.P.Val), ][1:10, ]
top_genes
library(pheatmap)

top50 <- rownames(results)[1:50]

mat_top <- mat[top50, ]

pheatmap(
  mat_top,
  scale = "row",
  show_rownames = FALSE,
  annotation_col = data.frame(Group = group)
)
common_genes <- intersect(rownames(results), rownames(mat))
length(common_genes)
top50 <- head(common_genes, 50)

mat_top <- mat[top50, , drop = FALSE]
library(pheatmap)

pheatmap(
  mat_top,
  scale = "row",
  show_rownames = FALSE
)

save.image("breast_cancer_analysis.RData")
write.csv(results, "DE_results.csv")
write.csv(mat, "expression_matrix.csv")
png("volcano_plot.png", width = 800, height = 600)
EnhancedVolcano(results,
                lab = rownames(results),
                x = "logFC",
                y = "adj.P.Val")
dev.off()
png("heatmap.png", width = 800, height = 600)
pheatmap(mat_top, scale = "row")
dev.off()