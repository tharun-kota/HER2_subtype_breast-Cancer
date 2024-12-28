## Differential gene expression between HER2 sub-type breast cancer


setwd('C:/Users/tharun kota/OneDrive/Desktop/Bio_principles/project/brca_tcga_pan_can_atlas_2018') # set the path
data_patient_path = paste(getwd(),"data_clinical_patient.txt", sep = "/") # load patient data
data_patient = read.delim(data_patient_path) 

# we will skip 5 rows of column descriptions. 
data_patient = data_patient[5:dim(data_patient)[1],]


# loading RNA data
path_RNA = paste(getwd(),"data_mrna_seq_v2_rsem.txt", sep = "/")
data_Rnaseq = read.delim(path_RNA)

# loading CNA data
path_cna = paste(getwd(),'data_cna.txt', sep = "/")
data_CNA = read.delim(path_cna)

# Cols 1 and 2 are gene names.
assay = round(as.matrix(data_Rnaseq[,-c(1,2)])) 
rownames(assay) = data_Rnaseq[,1]



# Find the row where the first column has "ERBB2"
ERBB2_row = which(data_CNA[, 1] == "ERBB2")

# Extract the ERBB2 counts (assuming patient IDs start from the 3rd column)
ERBB2_data = matrix(data_CNA[ERBB2_row, 3:ncol(data_CNA)], 
                    ncol = 1, 
                    dimnames = list(colnames(data_CNA)[3:ncol(data_CNA)], "ERBB2_Count"))

# matching and subsetting pateint_Id with RNA_seq_data & CNA_data
clean_colnames <- gsub("\\.$", "", gsub("\\.", "-", sub("(\\d{2})$", "", colnames(data_CNA)[3:length(colnames(data_CNA))]))) # replacing charachter's like (.,-)
clean_colnames<- as.data.frame(clean_colnames)
clean_colnames$id <- substr(clean_colnames$clean_colnames, 1, nchar(clean_colnames$clean_colnames) - 1) # removing charachter (-)
clean_colnames$rna_seq_ids <- colnames(data_CNA)[3:length(colnames(data_CNA))]
pat_id<- as.data.frame(data_patient$X.Patient.Identifier)
colnames(pat_id)<- 'id'

# matching and omiting rows
patient_ids <- merge(clean_colnames,pat_id,by='id',all=TRUE)
patient_ids <- na.omit(patient_ids)

# Extract matching column names from assay (those that are in ERBB2_data's row names)
matching_columns <- patient_ids$rna_seq_ids

# Subset assay matrix to only include the matching columns
assay <- assay[, colnames(assay) %in% matching_columns]
matching_ids <- rownames(ERBB2_data) %in% colnames(assay)

# Subset ERBB2_data to keep only the matching IDs
ERBB2_filtered <- ERBB2_data[matching_ids, , drop = FALSE]

# Initialize metadata
metadata = matrix(0, nrow = nrow(ERBB2_filtered), ncol = 1)
colnames(metadata) = "HER2_Status"

# Determine HER2 status based on ERBB2 counts
metadata[, 1] = ifelse(as.numeric(ERBB2_filtered[, "ERBB2_Count"]) > 0, 1, 0)



# Data tranformation using DESEQ2
library(DESeq2)
assay[is.na(assay)] = 0  # Impute with zeros the NA
assay[assay<0] = 0

# filter out genes with too many missing values.
smallestGroupSize = 3
keep = rowSums(assay >= 10) >= smallestGroupSize
assay = assay[keep,]

###
dds =  DESeqDataSetFromMatrix(countData = assay,
                              colData = metadata,
                              design = ~ HER2_Status) # data formating for DESEQ2
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
vsd = vst(dds) # variance stabilised transformed

# PCA & scatter plot
par(mfrow = c(1, 2))
plotPCA(vsd, intgroup=c("HER2_Status"))

#results
res = results(dds)
res[order(res$padj)[1:10],]# print Top 10 most differentially expressed



# pathway enrichment analysis
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# get subset of differentially expressed genes.
res_sig = res[res$padj<0.05,]
head(res_sig)
# separate into over and under expressed using log2foldchange
DE_over = rownames(res_sig[res_sig$log2FoldChange>0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange<0,])

go_results_over = enrichGO(
  gene          = DE_over,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

dotplot(go_results_over, showCategory=10) + ggtitle("Gene Ontology Enrichment over Expressed")


go_results_under = enrichGO(
  gene          = DE_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

dotplot(go_results_under, showCategory=10) + ggtitle("Gene Ontology Enrichment Under Expressed")


BiocManager::install("ReactomePA",force=TRUE)

library(ReactomePA)
library(pathview)

# we need to map into entrez for Reactome and Keggs


gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_entrez_under <- bitr(
  DE_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)



dotplot(kegg_results_over, showCategory=10) + ggtitle("Kegg Pathway Enrichment Over Expressed")



dotplot(kegg_results_under, showCategory=10) + ggtitle("Kegg Pathway Enrichment Under Expressed")

```





```{r}

reactome_results_over =  enrichPathway(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)

reactome_results_under =  enrichPathway(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)

dotplot(reactome_results_over, showCategory=10) + ggtitle("Reactome Pathway Enrichment Over Expressed")



dotplot(reactome_results_under, showCategory=10) + ggtitle("Reactome Pathway Enrichment Under Expressed")






kegg_results_under_pw = pairwise_termsim(kegg_results_under)
treeplot(kegg_results_under_pw)+ ggtitle("KEGG Enrichment Under Expressed")


top_DE = order(res$padj)

vsd_DE = assay(vsd)[top_DE[1:20],]

library(pheatmap)


annotation_colors = list(HER2_status = c(HER2_amplified = "blue1", Other = "#33a02c"))



annotation_col = data.frame(HER2_status= as.matrix(metadata[,1]))
rownames(annotation_col) = colnames(vsd)


pheatmap(
  vsd_DE,
  cluster_rows = TRUE,      
  cluster_cols = TRUE,  
  scale = 'row',
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col)




DEG_df <- data.frame(res[res$padj<0.05,])
write.xlsx(DEG_df, "DEG.xlsx", rowNames = TRUE)


library(glmnet)
library(survival)

surv_data <- subset(patient_ids, subset= rna_seq_ids %in% colnames(assay))
surv <- subset(
  data_patient[, c('X.Patient.Identifier', 'Overall.Survival.Status', 'Overall.Survival..Months.')], 
  X.Patient.Identifier %in% surv_data$id
)
# Convert the column to numeric and subset the data
surv <- subset(surv, subset = as.numeric(Overall.Survival..Months.) > 0)

C_type <- ERBB2_filtered[rownames(ERBB2_filtered) %in% pat_id3$rna_seq_ids, , drop = FALSE]
C_type <- ifelse(as.numeric(C_type[, "ERBB2_Count"]) > 0, 2, 1)

pat_id3 <- subset(patient_ids, subset =  patient_ids$id %in% surv$X.Patient.Identifier)
x <- vsd_DE[, colnames(vsd_DE) %in% pat_id3$rna_seq_ids]
x <-t(x)
y <- subset(surv[,c('Overall.Survival..Months.','Overall.Survival.Status')])
colnames(y) <- c('time','Status')
y$Status <- ifelse(y$Status == "0:LIVING", 0, 1)
y$Status <- gsub(",", "", y$status)
y$Status <- as.numeric(y$Status)
y$time <- as.numeric(y$time)
time <- as.numeric(y$time)
Status <- as.numeric(y$Status)

y <- Surv(time, Status)         # Survival response
# Fit the Cox model with glmnet
fit <- glmnet(x, y, family = "cox")
plot(fit)
res.cox <- survival::survfit(fit, s = 0.05, x = x, y = y) ~ C_type
plot(survival::survfit(fit, s = 0.05, x = x, y = y) ~ C_type)


library("survival")
library("survminer")



plot(surv_group1, col = "blue", lty = 1, xlab = "Time", ylab = "Survival Probability",
     main = "Survival Curves for Groups")
lines(surv_group2, col = "red", lty = 2)






