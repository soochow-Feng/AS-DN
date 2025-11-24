# AS-DN
1. Data Download and Preprocessing
GSE20129 Processing (Blood Samples)
r
# Load required packages
library(data.table)
library(tibble)
library(dplyr)
library(GEOquery)
library(limma)

# Download and process GSE20129
gset <- getGEO('GSE20129', destdir = '.', getGPL = FALSE)
exprSet <- exprs(gset[[2]])
pdata <- pData(gset[[2]])

# Normalization
exprSet <- normalizeBetweenArrays(exprSet)
exprSet <- as.data.frame(exprSet)

# Probe annotation using GPL6104 platform
anno <- fread("GPL6104-6.txt", data.table = FALSE, sep = '\t')
anno <- anno[, c(1, 6)]  # Select ID and gene symbol columns
anno <- anno[!anno$`ILMN_Gene` == "", ]

# Filter and annotate probes
tmp <- rownames(exprSet) %in% anno[, 1]
exprSet <- exprSet[tmp, ]
anno <- anno[match(rownames(exprSet), anno$ID), ]

# Keep probe with maximum expression for each gene
MAX <- by(exprSet, anno[, 2], 
          function(x) rownames(x)[which.max(rowMeans(x))])
MAX <- as.character(MAX)
exprSet <- exprSet[rownames(exprSet) %in% MAX, ]
rownames(exprSet) <- anno[match(rownames(exprSet), anno[, 1]), 2]

write.table(exprSet, file = 'Peripheral_Blood-GSE20129.txt', 
            sep = '\t', col.names = NA, quote = FALSE)
Additional Datasets Processing
Similar processing pipelines were applied to GSE43292 and GSE28829 datasets using their respective platform annotation files (GPL6244 and GPL570).

2. Differential Expression Analysis
r
# Differential analysis function
perform_DE_analysis <- function(expFile, conFile, treatFile, output_dir) {
  setwd(output_dir)
  
  # Read and process expression data
  rt <- read.table(expFile, sep = "\t", header = TRUE, check.names = FALSE)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp))
  rownames(data) <- rownames(exp)
  colnames(data) <- colnames(exp)
  data <- avereps(data)
  
  # Sample grouping
  sample1 <- read.table(conFile, sep = "\t", header = FALSE, check.names = FALSE)
  sample2 <- read.table(treatFile, sep = "\t", header = FALSE, check.names = FALSE)
  conData <- data[, as.vector(sample1[, 1])]
  treatData <- data[, as.vector(sample2[, 1])]
  rt <- cbind(conData, treatData)
  
  conNum <- ncol(conData)
  treatNum <- ncol(treatData)
  
  # Limma differential analysis
  Type <- c(rep("Control", conNum), rep("Treatment", treatNum))
  design <- model.matrix(~0 + factor(Type))
  colnames(design) <- c("Control", "Treatment")
  
  fit <- lmFit(rt, design)
  cont.matrix <- makeContrasts(Treatment - Control, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # Output results
  allDiff <- topTable(fit2, adjust = 'fdr', number = 200000)
  return(allDiff)
}

# Application example
de_results <- perform_DE_analysis("GSE96804_intersection.txt", 
                                 "Normal_samples.txt", 
                                 "DN_samples.txt", 
                                 "./DE_Analysis/")
3. Weighted Gene Co-expression Network Analysis (WGCNA)
r
library(WGCNA)

perform_WGCNA <- function(expFile, cliFile, output_dir) {
  setwd(output_dir)
  enableWGCNAThreads(30)  # Enable multi-threading
  
  # Data preparation
  rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp))
  rownames(data) <- rownames(exp)
  colnames(data) <- colnames(exp)
  data <- avereps(data)
  
  datExpr0 <- t(data)
  
  # Soft threshold selection
  powers <- c(1:30)
  sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
  softPower <- sft$powerEstimate
  
  # Network construction
  adjacency <- adjacency(datExpr0, power = softPower)
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM
  
  # Module identification
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              minClusterSize = 50, deepSplit = 2)
  dynamicColors <- labels2colors(dynamicMods)
  
  # Module-trait relationships
  MEs <- moduleEigengenes(datExpr0, colors = dynamicColors)$eigengenes
  cli <- read.table(cliFile, header = TRUE, sep = "\t", 
                   check.names = TRUE, row.names = 1)
  
  sameSample <- intersect(row.names(MEs), row.names(cli))
  MEs <- MEs[sameSample, ]
  datTraits <- cli[sameSample, ]
  
  moduleTraitCor <- cor(MEs, datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr0))
  
  return(list(moduleColors = dynamicColors, 
              moduleTraitCor = moduleTraitCor,
              moduleTraitPvalue = moduleTraitPvalue))
}

# WGCNA execution
wgcna_results <- perform_WGCNA("GSE96804_intersection.txt", 
                              "clinical_traits.txt", 
                              "./WGCNA_Analysis/")
4. Machine Learning Feature Selection (LASSO)
r
library(glmnet)

set.seed(123)  # For reproducibility

perform_LASSO <- function(inputFile, output_dir) {
  setwd(output_dir)
  
  # Data preparation
  rt <- read.table(inputFile, header = TRUE, sep = "\t", 
                  check.names = FALSE, row.names = 1)
  rt <- t(rt)
  
  x <- as.matrix(rt)
  y <- gsub("(.*)\\_(.*)", "\\2", row.names(rt))
  
  # LASSO with cross-validation
  cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1,
                    type.measure = 'deviance', nfolds = 10)
  
  # Feature selection
  coef <- coef(cvfit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  lassoGene <- row.names(coef)[index]
  lassoGene <- lassoGene[-1]  # Remove intercept
  
  return(lassoGene)
}

# LASSO application
selected_genes <- perform_LASSO("GeneExpression_data.txt", "./Machine_Learning/")
5. Biomarker Validation
r
library(ggplot2)
library(ggpubr)
library(reshape2)

validate_biomarkers <- function(expFile, typeFile, geneFile, output_dir) {
  setwd(output_dir)
  
  # Data preparation
  rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp))
  rownames(data) <- rownames(exp)
  colnames(data) <- colnames(exp)
  data <- avereps(data)
  data <- t(data)
  
  # Target genes expression
  key_genes <- read.table(geneFile, header = FALSE, 
                         sep = "\t", check.names = FALSE)[, 1]
  data <- data[, intersect(colnames(data), key_genes)]
  
  # Sample information
  type <- read.table(typeFile, sep = "\t", header = FALSE, 
                    check.names = FALSE, row.names = 1)
  colnames(type) <- "Group"
  
  # Prepare data for plotting
  sameSample <- intersect(row.names(data), row.names(type))
  plot_data <- cbind(data[sameSample, ], type[sameSample, ])
  colnames(plot_data)[ncol(plot_data)] <- "Group"
  plot_data <- as.data.frame(plot_data)
  plot_data[, 1:(ncol(plot_data) - 1)] <- lapply(plot_data[, 1:(ncol(plot_data) - 1)], as.numeric)
  plot_data <- melt(plot_data, id.vars = c("Group"))
  colnames(plot_data) <- c("Group", "Gene", "Expression")
  
  # Generate boxplots
  boxplot <- ggboxplot(plot_data, x = "Gene", y = "Expression", fill = "Group",
                      palette = c("#0073C2FF", "#EFC000FF"),
                      xlab = "", ylab = "Expression Level",
                      legend.title = "Condition") +
    stat_compare_means(aes(group = Group), method = "wilcox.test",
                      label = "p.signif") +
    theme_bw() +
    rotate_x_text(45)
  
  return(boxplot)
}

# Biomarker validation
biomarker_plot <- validate_biomarkers("GSE96804_intersection.txt",
                                     "sample_groups.txt",
                                     "candidate_genes.txt",
                                     "./Biomarker_Validation/")
6. ROC Curve Analysis
r
library(pROC)

perform_ROC_analysis <- function(expFile, geneFile, output_dir) {
  setwd(output_dir)
  
  # Data preparation
  rt <- read.table(expFile, header = TRUE, sep = "\t", 
                  check.names = FALSE, row.names = 1)
  y <- gsub("(.*)\\_(.*)", "\\2", colnames(rt))
  y <- ifelse(y == "Control", 0, 1)
  
  gene_list <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
  
  # ROC analysis for each gene
  for(gene in as.vector(gene_list[, 1])) {
    roc_result <- roc(y, as.numeric(rt[gene, ]))
    ci_result <- ci.auc(roc_result, method = "bootstrap")
    
    pdf(file = paste0("ROC_", gene, ".pdf"), width = 6, height = 5)
    plot(roc_result, print.auc = TRUE, col = "#2A9D8E", 
         legacy.axes = TRUE, main = gene)
    text(0.4, 0.4, paste0("95% CI: ", sprintf("%.3f", ci_result[1]), 
                          " - ", sprintf("%.3f", ci_result[3])), 
         col = "#2A9D8E")
    dev.off()
  }
}

# ROC analysis execution
perform_ROC_analysis("normalized_expression.txt", 
                    "biomarker_genes.txt", 
                    "./ROC_Analysis/")
Reproducibility Information
Software Versions
R version: 4.2.0 or higher

data.table: 1.14.2

limma: 3.52.0

WGCNA: 1.72

glmnet: 4.1

pROC: 1.18.0

ggplot2: 3.4.0

Seed Setting
All random processes (particularly in machine learning algorithms) used set.seed(123) to ensure reproducibility.

Data Availability
All datasets used are publicly available from GEO:

GSE20129, GSE43292, GSE28829, GSE96804, GSE142153
