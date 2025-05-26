## DEA 
library(DESeq2)
library(stats)
library(ggplot2)
library(genefilter)
library(pheatmap)
library(dplyr)
library(EnhancedVolcano)

#### 1. Load counts & metadata ####
counts <-read.delim("/users/k20049686/new_matrix_totalcounts.tsv", header = TRUE, sep = '\t', row.names = 1)
colData <- read.delim("/scratch/prj/ccc_dcis_ncrna/data/samples_batches.txt", header = TRUE, sep = '\t', row.names = 1)
gene_symbols <- read.delim("/users/k20049686/gene_symbols.tsv", header = T, sep = '\t', row.names = 1) 
# removing the decimal after ENS geneIDs 
rownames(counts) <- sub("\\..*", "", rownames(counts))
#remove total sums row from my_matrix
counts <- counts[-86403,]

#### 2. Filtering out samples for DEA ####
# primary tissue samples
colData <- colData[colData$tissue == "Primary",]

# ipsilateral DCIS needs to become cases
colData <- colData %>%
  mutate(case_control = ifelse(first_subseq_event=="ipsilateral DCIS" & is.na(case_control), "case", case_control))
  
colData %>%
  group_by(first_subseq_event) %>%
  summarise(total = n())

#subselecting just the 165 cases & controls samples (experimental approach)
keep_samples <- colData$case_control %in% c('case', 'control')
colData <- colData[keep_samples,]
keep_cols <- rownames(colData)
counts <- counts[,keep_cols]

all(colnames(counts) == rownames(colData))

# saving list of DEA samples
write.table(colData, "DEA_samples.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
DEA_samples <-read.delim("/users/k20049686/DEA_samples.tsv", header = TRUE, sep = '\t', row.names = 1)


writeLines(rownames(colData), "DEA_sampleIDs.txt")





#### Cleaning - dup have been removed and gene_symbols updated ####
# 3 duplicate gene names in 'gene_symbols'
gene_symbols[duplicated(rownames(gene_symbols)),]

gene_symbols <- gene_symbols[-c(30165,59108,63353,67410),]
"ENSG00000290758" %in% rownames(my_matrix) #TRUE
my_matrix[rownames(my_matrix)=="ENSG00000290758",]

"ENSG00000272655" %in% rownames(my_matrix)
my_matrix[rownames(my_matrix)=="ENSG00000272655",] #0 counts so POLR2J4 removed 

trans_allmetrics <- read.delim("/users/k20049686/allmetrics_t.tsv", header = TRUE, sep = '\t')
allmetrics$unique_mapped_rate <- round(allmetrics$unique_mapped_rate, 1)
write.table(allmetrics, "allmetrics_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

length(rownames(colData))
DEA_sampleIDs <- rownames(colData)
allmetrics <- allmetrics[rownames(allmetrics) %in% DEA_sampleIDs,]

all(rownames(allmetrics) %in% DEA_sampleIDs)

# saving the QC metrics of thr 154 samples 
write.table(allmetrics, "DEA_samples_QCmetrics.tsv", sep = "\t", quote = FALSE, row.names = TRUE)


######################################################################
# removing the 35 outliers

my_excluded_samples <- sort(c("S110859", "S110883", "S110868", "S110872", "S110892", "S110874", "S110761", 
                              "S110864", "S110866", "S110885", "S110895", "S110760", "S110894", "S110854", 
                              "S110940", "S110754", "S110865", "S110855", "S110766", "S110762", 
                              "S110751", "S110797","S110879", "S110891", "S110880", "S110852", "S110763", "S110851", 
                              "S110884", "S110853", "S110856", "S110893", "S110764", "S110850", "S110882"))

counts <- subset(counts, select = !(names(counts) %in% my_excluded_samples))
colData <- colData[!(rownames(colData) %in% my_excluded_samples),]

all(my_excluded_samples %in% rownames(colData))
all(my_excluded_samples %in% names(counts))

######################################################################

#### 3. Refactor w/ correct levels #### 
colData$batch <- factor(colData$batch)
colData$case_control <- factor(colData$case_control, levels = c("control","case")) #control as ref 



#### 4. Creating deseq2 dataset #### 

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ batch + case_control)

levels(dds$case_control)
print(dds)



#### 5. Feature selection - removing genes w/ low counts #### 
  # at least 5 reads in at least 10 samples 
keep <- rowSums(counts(dds) >= 5) >= 10
dds <- dds[keep,]
dim(dds) 

dds <- DESeq(dds)

#### 6. Contrast for the condition & order by increasing p-value ####
  # table of DE status of case samples compared to controls (ref)
DEresults <- results(dds, contrast = c("case_control", "case", "control")) #ensures comparison is case to control (ref)
DEresults <- DEresults[order(DEresults$pvalue),]
DEresults_df <- as.data.frame(DEresults)
print(DEresults)
mcols(DEresults)$description
# baseMean = norm GE values relative to all samples
# +ve lfc = upreg in cases 
# lfcSE = lfc standard error # stat = used to cal pvalue # padj = pvalue adj for multiple testing 

#### 7. Queries ####

summary(DEresults)
sum(DEresults_df$padj < 0.1, na.rm = TRUE)
sum(DEresults_df$padj < 0.05, na.rm = TRUE)


## filter genes by adj p-value & lFC
dim(DEresults) # 40,896 genes
dim(filtered) # 128 genes 
head(filtered)

# filtering for DE genes by padj & FC 
filtered <- DEresults_df %>%
  filter(DEresults_df$padj < 0.05)

filtered <- filtered %>%
  filter(abs(filtered$log2FoldChange) > 1 )

# adding gene symbols to filtered DE genes 
all(rownames(gene_symbols) %in% filtered_genes) 
filtered_genes <- rownames(filtered)

gene_symbols <- gene_symbols[,c(1,2)]
gene_symbols <- gene_symbols[rownames(gene_symbols) %in% filtered_genes,]

filtered <- merge(filtered,gene_symbols, by = "row.names")
colnames(filtered)[10] <- "gene_symbol"
summary(duplicated(filtered$gene_symbol)) # 47 NAs

# filled blanks w/ NAs 
filtered$gene_symbol[filtered$gene_symbol==""] <- NA

filtered %>%
  group_by(gene_symbol) %>%
  summarise(total = n())


# saving filtered DE genes 
write.table(filtered, "filtered_DEgenes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
filtered_DEgenes <- read.delim("/users/k20049686/filtered_DEgenes.tsv", header = T, sep = '\t', row.names = 1)

# doing same for DEresults
all(rownames(gene_symbols) %in% rownames(DEresults_df))

rownames <- rownames(DEresults_df)
gene_symbols <- gene_symbols[,c(1,2)]
gene_symbols <- gene_symbols[rownames(gene_symbols) %in% rownames,]

DEresults_df <- merge(DEresults_df, gene_symbols, by = "row.names")
colnames(DEresults_df)[10] <- "gene_symbol"
summary(duplicated(DEresults_df$gene_symbol)) # 12,532 NAs

rownames(DEresults_df) <- DEresults_df[,1]
DEresults_df$Row.names <- NULL

DEresults_df %>%
  group_by(diffexpressed) %>%
  summarise(total = n())

#### 8. VISUALISATIONS ####

# check normalisation - MA plot - rel between lfc (M values) x mean norm counts (A values) 
DESeq2::plotMA(object = dds) # before shrinkage

# need to shrink lFC estimates to reduce noise & improve accuracy & visualisation
resultsNames(dds) # checking correct coeff used 
lfcshrunk_dds <- lfcShrink(dds, coef = "case_control_case_vs_control", type = "ashr") #ashr less aggressive than apeglm 

# improved MA plot
DESeq2::plotMA(object = lfcshrunk_dds)


### p-value distribution plot
ggplot(data = DEresults_df, aes(x = pvalue)) +
  geom_histogram(bins=100)

# dispersion plot - not normal
plotDispEsts(dds)

# volcano plot
dds <- as.data.frame(dds)

# label the genes
filtered$diffexpressed <- "NO" # new col - default is not DE
filtered$diffexpressed[filtered$log2FoldChange > 0.1 & filtered$padj < 0.05] <- "UP" # these genes labelled upreg
filtered$diffexpressed[filtered$log2FoldChange < 0.1 & filtered$padj < 0.05] <- "DOWN" # downreg genes

filtered$delabel <- NA

# plotting 
  # noise removed - lfcshrunk_dds
  # noise not removed - DEresults
  # just filtered genes - filtered df
ggplot(data = DEresults_df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_point() +
  xlim(-5,5) +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c('blue', 'black', 'red')) +
  theme(text = element_text(size = 20))

# enhanced volcano
  # default pCutoff = 0.00001 and FCcutoff = |1| (i think)
EnhancedVolcano(DEresults_df,
                lab = DEresults_df$gene_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5,5),
                pCutoff = 0.05,
                drawConnectors = FALSE)


######################################################


### extracts normalised counts from dds object
norm_dds <- DESeq2::counts(dds, normalized = TRUE)


#select top 500 var gene names 
selected_genes <- names(sort(apply(norm_dds,1,var),
                             decreasing = TRUE)[1:500])

### variance stabilisation - vst faster than rlog 
vsd <- vst(dds, blind = FALSE)

#plot PCA of top 500 most variable genes 
plotPCA(vsd, intgroup="case_control", ntop = 500) +
  ylim(-50,50) +
  theme_bw() # or other metadata column

####### another pca method
pca <- prcomp(t(assay(vsd)))
summary(pca$x)

# scree plot - decide no of PCs
plot(pca, type = "l")

#### plotting more PCs 
# pairiwise plots

pc_df <- as.data.frame(pca$x)
pc_df$batch <- colData(vsd)$batch
pc_df$case_control <- colData(vsd)$case_control

ggplot(pc_df, aes(x = PC3, y = PC4, colour = case_control)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA: PC3 vs PC4", x = "PC3", y = "PC4")

### heatmap

# this method identifies the top 20 most vairable genes - not necessarily most DE
  # high variance genes may not be bio relevant or stat significant
topVarGenes <- head(order(-rowVars(assay(vsd))), 20)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
annotation <- as.data.frame(colData(vsd)[, "case_control", drop = FALSE])
pheatmap(mat, annotation_col = annotation)


# selects top 20 most stat sign DE genes for case/control comparisons
  # can filter more by lFC
top_hits <- DEresults_df[order(DEresults_df$padj), ][1:20,] #ordered by adj pvalue & selected top 20
top_hits <- row.names(top_hits)
top_hits

pheatmap(assay(vsd)[top_hits,], 
        annotation_col = annotation, 
        cluster_rows = FALSE,
        cluster_cols = FALSE)


##################################################
# heatmap of sample-to-sample distance matrix 
sampleDists <- dist(t(assay(vsd)))








