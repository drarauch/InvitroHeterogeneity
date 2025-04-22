### Figure 2
# The  files are provided at the open science framework https://osf.io/s8nfb/ 

### Figure 2B
# UMAP plot
scobject <- qs::qread("SeuratObject.qs")
Result <- readRDS("MarkerGeneList.Rds")
library(Seurat)

DimPlot(scobject, raster=F)

# Stats on enriched and exclusive Markers
mat <- matrix(NA,ncol=3,nrow=length(levels(scobject$seurat_clusters)))
rownames(mat) <- paste("Cl_",levels(scobject$seurat_clusters),sep="")
colnames(mat) <- c('Exclusive','Enriched','Significant')
for (i in 1:length(Clusters)){
  tmp <- Result[[i]]
  mat[i,1] <- nrow(tmp[tmp$Marker =="Exclusive",])
  mat[i,2] <- nrow(tmp[tmp$Marker =="Enriched",])
  mat[i,3] <- nrow(tmp[tmp$Marker =="Significant",])
}
mat

# output
#         Exclusive Enriched Significant
#	0        11       53         469
#	1        30      109         489
#	2        26       51         498
#	3        63      116         248
#	4        18       79         415
rm(mat,scobject,Result)



### Figure 2C Quality metrices for the clusters
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)

# Violin Plot of quality measures
VlnPlot(scobject, "detected", pt.size = 0)
VlnPlot(scobject, "DoubletScore", pt.size = 0)
VlnPlot(scobject, "subsets_Mito_percent", pt.size = 0)

rm(scobject)



### Figure 2D Show cellular composition and order Donor based on number of cells in the data set
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)
library(scales)
library(gplots)
library(fields)

# Cluster abundance across donor samples
tbl <- table(scobject$seurat_clusters, scobject$Dataset)
for (i in 1:ncol(tbl)){
  tbl[,i] <- tbl[,i] / as.numeric(table(scobject$Dataset)[i])
}

# Number of cells per donor to order tbl
vec <- table(scobject$Dataset)[order(table(scobject$Dataset))]

# Bar plot with cluster abundancy
barplot(tbl[,names(vec)],legend=F, col=hue_pal()(nrow(tbl)), las=2, horiz=T)

# Heatmap for the number of cells per library to align with the barplot (done in illustrator)
mat_col <- designer.colors(n=1000, col=c('white','black'))
mat_col_breaks <- seq(0,max(vec),length=1001)
heatmap.2(cbind(vec,vec), breaks= mat_col_breaks, col=mat_col,main="",dendrogram="none", Rowv = F, Colv=F, scale='none', trace='none')

rm(mat_col, mat_col_breaks, vec, tbl, i, scobject)



### Figure 2E expression of MSC related or qualifying genes
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)

# Dot Plot of known MSC related genes
DP <- DotPlot(scobject, features=c("PDGFRA", "NT5E", "CD34", "ENG", "THY1", "CXCL12", "GREM1", "NGFR", "MCAM", "CXCR4", "LEPR", "SDC2"))
DP

# Generate dendogram
data <- DP$data
scaled_expression <- cbind(data[1:13,c(5)], data[14:26,c(5)], data[27:39,c(5)], data[40:52,c(5)], data[53:65,c(5)])
colnames(scaled_expression) <- 0:4
rownames(scaled_expression) <- rownames(data)[1:13]
HM <- heatmap.2(as.matrix(scaled_expression))
plot(HM$rowDendrogram)

# Combine dendrogram and dot plot in Illustrator
rm(DP,HM,scaled_expression,data,scobject)



### Figure 2F Marker gene expression across clusters, focusing on exclusive genes
scobject <- qs::qread("SeuratObject.qs")
Result <- readRDS("MarkerGeneList.Rds")
library(Seurat)
library(dplyr)

markers <- bind_rows(Result, .id = "column_label")
# Plot heatmap
DoHeatmap(scobject, features=markers[markers$Marker=="Exclusive",'Gene'],assay = "originalexp")

rm(scobject,markers,Result)



### Figure 2G Dot Plot of cluster-specific genes
scobject <- qs::qread("SeuratObject.qs")
Result <- readRDS("MarkerGeneList.Rds")
library(Seurat)
library(dplyr)

markers <- bind_rows(Result, .id = "column_label")
markers <- markers[markers$Marker=="Exclusive",]
markers <- markers[order(markers$column_label,markers$logFC),]

DotPlot(scobject, features=c("SPARC","FGF7","IGFBP7","LUM","IGFBP5","PI16","S100A4","MFAP5","H2AZ1","MT2A","CENPF","MAP2K5","CEMIP","NEK7","GLS","SMURF2","HFM1","ASIP","TPM2","SEC61G"))

rm(markers, Result, scobject)



### Figure 2H Gene ontology analysis based on enriched genes
All <- readRDS("AllGene.Rds")
Result <- readRDS("MarkerGeneList.Rds")
library(fields)
library(gplots)
library(goseq)
library(dplyr)

# Extract cluster specific genes and focus on the enriched ones
markers <- bind_rows(Result, .id = "column_label")
markers <- markers[markers$Marker=="Enriched",]
Clusters <- c("Cl_0","Cl_1","Cl_2","Cl_3","Cl_4")

Gene_groups <- list()
for (i in 1: length(Clusters)){
	Gene_groups[[i]] <- markers[markers$column_label == Clusters[i],'Gene']
}
names(Gene_groups) <- Clusters

## Go category of clusters
# Initial round with 1st cluster
Cl <- Gene_groups[[1]]
tmp_Cl <- make.names(Cl, unique = TRUE)

# Extracting all measured genes from scRNA-seq
tmp_All <- make.names(All, unique = TRUE)
# Converting DE genes to intergers, significant DE genes = 1, non-significant genes = 0
Cl_man <- as.integer(tmp_All %in% tmp_Cl)
names(Cl_man) <- tmp_All

# Fitting the probability weighting (PWF) function while correcting for gene count length
pwf_Cl <- nullp(Cl_man, "hg19", "geneSymbol")
# GO analysis
# Here, the test is limited to biological processes (GO:BP), change to GO:CC or GO:MF for cellular components or molecular function testing
GO.BP_Cl <- goseq(pwf_Cl, "hg19", "geneSymbol", test.cats = c("GO:BP"), use_genes_without_cat = TRUE)
# Applying multiple testing by Benjamini Hochberg to overenriched p values
GO.BP_Cl$over_represented_p.adjust <- p.adjust(GO.BP_Cl$over_represented_pvalue, method = "BH")
#Collecting GO terms
GO_cluster <- GO.BP_Cl[,c(1,6,8)]
names(GO_cluster)[3] <- names(Gene_groups)[1]
# Loop for the other clusters
for(i in 2:length(Gene_groups)){
  Cl <- Gene_groups[[i]]
  tmp_Cl <- make.names(Cl, unique = TRUE)
  if(length(tmp_Cl)>0){
    Cl_man <- as.integer(tmp_All %in% tmp_Cl)
    names(Cl_man) <- tmp_All
    pwf_Cl <- nullp(Cl_man, "hg19", "geneSymbol")
    GO.BP_Cl <- goseq(pwf_Cl, "hg19", "geneSymbol", test.cats = c("GO:BP"), use_genes_without_cat = TRUE)
    GO.BP_Cl$over_represented_p.adjust <- p.adjust(GO.BP_Cl$over_represented_pvalue, method = "BH")
    names(GO.BP_Cl)[8] <- names(Gene_groups)[i]
    GO_cluster <- merge(GO_cluster,GO.BP_Cl[,c(1,8)],by="category")
  } else {}
}
# Transform p-values to log scale and reduce to significant pathways
p <- cbind(GO_cluster[,1:2],-log10(GO_cluster[,c(3:ncol(GO_cluster))]))
p <- p[apply(p[,3:ncol(p)],1,max)>2,]

# Select pathways of interest, i.e. 
keep <- c("microtubule cytoskeleton organization","mitotic cell cycle","nuclear division","angiogenesis","cytoplasmic translation",
"rRNA processing","apoptotic process","microtubule-based process","chromosome segregation","cell adhesion","cell death",
"tissue development","gene expression","actin cytoskeleton organization","extracellular matrix organization","collagen fibril organization",
"ribosome biogenesis","extracellular structure organization","organelle fission","chromosome organization","supramolecular fiber organization")
p <- p[complete.cases(p),]
rownames(p) <- p$term
p <- p[keep,3:7]

# Apply a simple binary clustersing based on significance or not
q <- p
q[q < 0.05] <- 0
q[q > 0.05] <- 1
hr <- hclust(dist(q, method="binary"))

# Cutt off significances above 15 (for visualization)
p[p>15] <- 15

# Do heatmap
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(1,max(p),length=51))
heatmap.2(as.matrix(p),main="", Rowv = as.dendrogram(hr), Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', cexRow=0.5)

rm(Result,hg19.geneSymbol.LENGTH,q,p,keep,GO_cluster,GO.BP_Cl,pwf_Cl,Cl,Cl_man,tmp_Cl,tmp_All,All,i, mat_col, mat_col_breaks, Gene_groups, Clusters, markers)



### Figure 2I rWIKI-pathway analysis. NOTE: As WikiPathways is updated results may change, here we use "wikipathways-20240101-gmt-Homo_sapiens.gmt"
Result <- readRDS("MarkerGeneList.Rds")
library(dplyr)
library(fields)
library(gplots)

# Prepare for WikiPathways analysis
gmt <- "wikipathways-20240101-gmt-Homo_sapiens.gmt"
wp2gene <- clusterProfiler::read.gmt(gmt)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name", "version", "wpid", "org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME

library(org.Hs.eg.db)
Pathways <- list()

# Loop through all clusters
for (i in 1:length(Result)) {
  # Extract list
  Tmp <- Result[[i]]
  
  # Convert gene names
  Enriched_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Enriched",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  Exclusive_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Exclusive",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  Significant_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Significant",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Do the pathway analysis
  if (sum(wp2gene$gene %in% Enriched_Entrez[[2]]) > 0) { wiki_enriched <- clusterProfiler::enricher(Enriched_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_enriched <- data.frame() }
  if (sum(wp2gene$gene %in% Exclusive_Entrez[[2]]) > 0) { wiki_exclusive <- clusterProfiler::enricher(Exclusive_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_exclusive <- data.frame() }
  if (sum(wp2gene$gene %in% Significant_Entrez[[2]]) > 0) { wiki_significant <- clusterProfiler::enricher(Significant_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_significant <- data.frame() }
  
  # Adding gene symbols to the resulting pathway file
  if (nrow(wiki_enriched) > 0) { wiki_enriched <- as.data.frame(DOSE::setReadable(wiki_enriched, org.Hs.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive <- as.data.frame(DOSE::setReadable(wiki_exclusive, org.Hs.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_significant) > 0) { wiki_significant <- as.data.frame(DOSE::setReadable(wiki_significant, org.Hs.eg.db, keyType = "ENTREZID")) }
  
  # Calculate gene ratios (i.e. the fraction of genes in the pathway)
  if (nrow(wiki_enriched) > 0) { wiki_enriched$Enriched_GeneRatio <- as.numeric(substr(wiki_enriched$GeneRatio, 1, regexpr("/", wiki_enriched$GeneRatio)-1)) / as.numeric(substr(wiki_enriched$GeneRatio, regexpr("/", wiki_enriched$GeneRatio)+1, nchar(wiki_enriched$GeneRatio))) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive$Exclusive_GeneRatio <- as.numeric(substr(wiki_exclusive$GeneRatio, 1, regexpr("/", wiki_exclusive$GeneRatio)-1)) / as.numeric(substr(wiki_exclusive$GeneRatio, regexpr("/", wiki_exclusive$GeneRatio)+1, nchar(wiki_exclusive$GeneRatio))) }
  if (nrow(wiki_significant) > 0) { wiki_significant$Significant_GeneRatio <- as.numeric(substr(wiki_significant$GeneRatio, 1, regexpr("/", wiki_significant$GeneRatio)-1)) / as.numeric(substr(wiki_significant$GeneRatio, regexpr("/", wiki_significant$GeneRatio)+1, nchar(wiki_significant$GeneRatio))) }
  
  # Calculate pathway ratios (i.e. the fraction of the pathway in the genes)
  if (nrow(wiki_enriched) > 0) { wiki_enriched$Enriched_PathRatio <- as.numeric(substr(wiki_enriched$GeneRatio, 1, regexpr("/", wiki_enriched$GeneRatio)-1)) / as.numeric(substr(wiki_enriched$BgRatio, 1, regexpr("/", wiki_enriched$BgRatio)-1)) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive$Exclusive_PathRatio <- as.numeric(substr(wiki_exclusive$GeneRatio, 1, regexpr("/", wiki_exclusive$GeneRatio)-1)) / as.numeric(substr(wiki_exclusive$BgRatio, 1, regexpr("/", wiki_exclusive$BgRatio)-1)) }
  if (nrow(wiki_significant) > 0) { wiki_significant$Significant_PathRatio <- as.numeric(substr(wiki_significant$GeneRatio, 1, regexpr("/", wiki_significant$GeneRatio)-1)) / as.numeric(substr(wiki_significant$BgRatio, 1, regexpr("/", wiki_significant$BgRatio)-1)) }
  
  # Set column names
  if (nrow(wiki_enriched) > 0) { colnames(wiki_enriched)[c(8,9,11)] <- c("Enriched_Pvalue","Enriched_FDR","Enriched_Genes") }
  if (nrow(wiki_exclusive) > 0) { colnames(wiki_exclusive)[c(8,9,11)] <- c("Exclusive_Pvalue","Exclusive_FDR","Exclusive_Genes") }
  if (nrow(wiki_significant) > 0) { colnames(wiki_significant)[c(8,9,11)] <- c("Significant_Pvalue","Significant_FDR","Significant_Genes") }
  
  # Combine the results
  significant <- c(wiki_significant[ wiki_significant$Significant_FDR <= 0.1,1],wiki_enriched[ wiki_enriched$Enriched_FDR <= 0.1,1],  wiki_exclusive[ wiki_exclusive$Exclusive_FDR <= 0.1,1])
  if (length(significant) > 0) {
    # Subset from all pathways to only significant ones
    all <- wpid2name[ wpid2name[,1] %in% significant,]
    all <- all[ duplicated(all[,1])==F,]
    colnames(all) <- c("ID", "Description")
    
    # Merge statistics
    if (nrow(wiki_significant) > 0) { all <- merge(all, wiki_significant[,c(1,8,9,11)], all.x=T, by="ID") }
    if (nrow(wiki_enriched) > 0) { all <- merge(all, wiki_enriched[,c(1,8,9,11)], all.x=T, by="ID") }
    if (nrow(wiki_exclusive) > 0) { all <- merge(all, wiki_exclusive[,c(1,8,9,11)], all.x=T, by="ID") }
    
    # Handle NAs appropriately (ratios to 0, pvalues/FDRs to 1 and gene lists to "")
     
    for (m in grep("Pvalue", colnames(all))) { all[ is.na(all[,m]),m] <- 1 }
    for (m in grep("FDR", colnames(all))) { all[ is.na(all[,m]),m] <- 1 }
    for (m in grep("Genes", colnames(all))) { all[ is.na(all[,m]),m] <- "" }
    
    # Store the results
    Pathways[[length(Pathways)+1]] <- all
    names(Pathways)[i] <- names(Result)[i]
  }
}

# Collect information in a "plotable" format

for (i in 1:length(Pathways)){
	tmp <- Pathways[[i]]
	tmp <- unique(tmp[tmp$Enriched_FDR < 0.1 | tmp$Exclusive_FDR < 0.1,c(2,7,10)])
	colnames(tmp)[2:3] <- c(paste(names(Pathways[i]),c("_Enriched_FDR","_Exclusive_FDR"),sep=""))
	if (i ==1){
		p <- tmp
	} else {
		p <- merge(p, tmp, by="Description", all.x=T, all.y=T)
	}
}

# Handle NAs appropriately (FDRs to 1)
for (m in grep("FDR", colnames(p))) { p[ is.na(p[,m]),m] <- 1 }

# Focus on enriched genes, and patways with FDR < 0.05
rownames(p) <- p$Description
p <- p[,grep("Enriched",colnames(p))]
p <- p[apply(p,1,min) < 0.05,]

# Remove redundant terms and those completely related to other cell types (muscle, neurons, or cancer)
keep <- c("Aerobic glycolysis","Alpha 6 beta 4 signaling pathway","Angiotensin II receptor type 1 pathway","Burn wound healing","Cell migration and invasion through p75NTR",
"Cori cycle","Cytoplasmic ribosomal proteins","EGF EGFR signaling pathway","ErbB signaling pathway","Focal adhesion","Glycolysis and gluconeogenesis","Glycolysis in senescence","Hippo signaling regulation pathways",
"Inflammatory response pathway","Interactions between LOXL4 and oxidative stress pathway","Microtubule cytoskeleton regulation","miRNA targets in ECM and membrane receptors","Neural crest cell migration during development","NRF2 pathway",
"Photodynamic therapy induced HIF 1 survival signaling","PI3K Akt signaling pathway","Role of hypoxia angiogenesis and FGF pathway in OA chondrocyte hypertrophy",
"TGF beta signaling pathway","TGF Smad signaling pathway","VEGFA VEGFR2 signaling","Wnt signaling")

p <- p[keep,]

# Apply a simple binary clustersing based on significance or not
q <- p
q[q < 0.05] <- 0
q[q > 0.05] <- 1
hr <- hclust(dist(q, method="binary"))

# Transform into log scale and cutt off significances above 15 (for visualization)
p <- -log10(p)
p[p>15] <- 15

# Do heatmap
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(1,max(p),length=51))
heatmap.2(as.matrix(p),main="", Rowv = as.dendrogram(hr), Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', cexRow=0.5)

rm(Result,Pathways,mat_col, mat_col_breaks, p, q, hr, keep, m,i,wiki_significant,wiki_exclusive,wiki_enriched,all,Clusters,Enriched_Entrez,Exclusive_Entrez,gmt,significant,Significant_Entrez,tmp,Tmp,wp2gene,wpid2gene,wpid2name)



### Figure 2K checking the enriched markers genes for eBMD associated SNPs in the vicinity of the TSS
All <- readRDS("AllGene.Rds")
Result <- readRDS("MarkerGeneList.Rds")
# GWAS summary statistics from http://www.gefos.org/?q=content/data-release-2018 doi:10.1038/s41588-018-0302-x
GWAS_summary <- fread("Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",h=T)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(gplots)
library(fields)

# Define TSS of genes from RNA-seq on hg19 coordinates
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Get transcripts for hg19
GR <- data.frame(transcripts(txdb))
# Get symbols to UCSC IDs
GR2 <- select(Homo.sapiens,GR$tx_name, "SYMBOL","TXNAME")
GR <- merge(GR,GR2, by.x="tx_name", by.y="TXNAME")
rm(GR2)
GR <- GR[GR$seqnames %in% paste('chr',c(1:23,'X','Y','M'),sep=""),]
# Reduce to the gene expressed in the scRNAseq dataset
GR <- unique(GR[GR$SYMBOL %in% All,])
# Get lowest TSS for "+"-stranded
GR1 <- GR[GR$strand =="+",] 
GR1 <- GR1[order(GR1$SYMBOL,GR1$start),]
GR1$dup <- duplicated(GR1$SYMBOL)
GR1 <- GR1[GR1$dup =="FALSE",]
GR1$dup <- NULL
GR1 <- GR1[,c('SYMBOL','seqnames','start','strand')]
names(GR1) <- c("Symbol",'Chr',"TSS","Strand")
# Get highest TSS for "-"-stranded
GR2 <- GR[GR$strand =="-",] 
GR2 <- GR2[order(GR2$SYMBOL,-GR2$end),]
GR2$dup <- duplicated(GR2$SYMBOL)
GR2 <- GR2[GR2$dup =="FALSE",]
GR2$dup <- NULL
GR2 <- GR2[,c('SYMBOL','seqnames','end','strand')]
names(GR2) <- c("Symbol",'Chr',"TSS","Strand")
GR <- rbind(GR1,GR2)
rm(GR1,GR2, txdb)

# Reformat GWAS summary statistics
GWAS_summary <- GWAS_summary[,c('RSID','CHR','BP','P.NI','BETA')]
names(GWAS_summary)[1] <- 'SNP'
GWAS_summary$Chr <- paste("chr",GWAS_summary$CHR, sep="")
GWAS_summary$Start <- as.numeric(GWAS_summary$BP)
GWAS_summary <- GWAS_summary[complete.cases(GWAS_summary$Start),c('SNP','P.NI','Chr','Start','BETA')]
GWAS_summary <- data.frame(GWAS_summary)
GWAS_summary$Pval <- as.numeric(GWAS_summary$P.NI)

# Do overlap of SNPs in window of 5 Mb
grGenes <- with(unique(GR[,c("Symbol","Chr","TSS")]) , GRanges(Chr, IRanges(start=TSS - 5000000, end=TSS + 5000000, names=Symbol)))
grSummary <- with(unique(GWAS_summary[,c("SNP","Chr","Start")]) , GRanges(Chr, IRanges(start=Start, end=Start, names=SNP)))

hits = findOverlaps(grGenes,grSummary)
tmp2 <- cbind(data.frame(ranges(grGenes)[queryHits(hits)]),data.frame(ranges(grSummary)[subjectHits(hits)]))
colnames(tmp2) <- c('TSSminus','TSSplus','Window','Symbol','SNP_Start','SNP_End','SNP_Length','SNP')
head(tmp2)

# Calculate distance bewteen TSS and SNP
tmp2$Distance <- tmp2$TSSminus+5000000 - tmp2$SNP_Start
tmp2 <- tmp2[,c('Symbol','SNP','Distance')]
tmp2 <- merge(tmp2[,c('Symbol','SNP','Distance')], GWAS_summary[,c("SNP","Pval","BETA")], by="SNP")
tmp2 <- merge( GR,tmp2[,c('Symbol','SNP','Pval','Distance')], by="Symbol")

#Chisquare test for enrichment with different distances from the TSS of genes from the scRNA-seq clusters (focus on enriched genes extracted from the list "Result")
mat2 <- matrix(NA, ncol=5, nrow=13)
x <- 1
for (k in c(50000,100000,250000,seq(500000,5000000, length=10))){
  for (i in 1:length(c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4'))){
	tmp <- Result[[i]]
	tmp <- tmp[tmp$Marker=="Enriched",'Gene']
    a <- length(tmp2[tmp2$Symbol %in% tmp & tmp2$Pval < 5E-8 & abs(tmp2$Distance) < k ,'SNP'])
    b <- length(tmp2[tmp2$Symbol %in% tmp & tmp2$Pval > 5E-8 & abs(tmp2$Distance) < k,'SNP'])
    c <- length(tmp2[tmp2$Pval < 5E-8 & abs(tmp2$Distance) < k ,'SNP'])
    d <- length(tmp2[tmp2$Pval > 5E-8 & abs(tmp2$Distance) < k ,'SNP'])
    M <- as.table(rbind(c(a, b), c(c, d)))
    dimnames(M) <- list(group = c("in_cluster", "in_gwascat"),
                        categ = c("HBMD_assoc","not_HBMD_assoc"))
    mat2[x,i] <- chisq.test(M)$p.value
  }
  x <- x+1
}
rownames(mat2) <- paste('Dist_',c(50000,100000,250000,seq(500000,5000000, length=10)),sep="")
colnames(mat2) <- c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4')

saveRDS(mat2,file="eBMD_Heterogeneity.rds")

mat2 <- -log10(mat2[1:7,])
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(5,max(mat2),length=51))
heatmap.2(as.matrix(mat2),Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=rownames(mat2),labCol=colnames(mat2) )

rm(mat2,a,b,c,d,tmp,k,i,GR,mat_col, mat_col_breaks, GWAS_summary,hits,grGenes, grSummary,Result, All)


### Figure 2L Relation of cluster specific genes to adipocyte and osteoblast differentiation
Result <- readRDS("MarkerGeneList.Rds")
All <- readRDS("AllGene.Rds")
# Alignment with expression data (https://doi.org/10.1038/s41588-019-0359-1), hypergeometric test to test overlap with osteogenic and adiopgenic induced genes 
TERT <- read.delim("GeneExpression_TERT_ATERT.txt",h=T)
library(dplyr)
library(fields)
library(gplots)

# Combine markers of single RNA-seq in one data frame and split enriched markers in a list
markers <- bind_rows(Result, .id = "column_label")

Gene_groups <- list()
Gene_groups[[1]] <- markers[markers$column_label== "Cl_0" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[2]] <- markers[markers$column_label== "Cl_1" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[3]] <- markers[markers$column_label== "Cl_2" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[4]] <- markers[markers$column_label== "Cl_3" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[5]] <- markers[markers$column_label== "Cl_4" & markers$Marker =="Enriched",'Gene' ]

names(Gene_groups) <- c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4')

# Group the osteoblast and adipocyte TERT genes in a list
Gene_groups_2 <- list()
Gene_groups_2[[1]] <- TERT[TERT$pVal_RNA_Ob7d_vs_Msc0h < 0.01 & TERT$log2FC_Ob7d < 0,'SYMBOL']
Gene_groups_2[[2]] <- TERT[TERT$pVal_RNA_Ob7d_vs_Msc0h < 0.01 & TERT$log2FC_Ob7d > 0,'SYMBOL']
Gene_groups_2[[3]] <- TERT[TERT$pVal_RNA_Ad7d_vs_Msc0h < 0.01 & TERT$log2FC_Ad7d < 0,'SYMBOL']
Gene_groups_2[[4]] <- TERT[TERT$pVal_RNA_Ad7d_vs_Msc0h < 0.01 & TERT$log2FC_Ad7d > 0,'SYMBOL']
names(Gene_groups_2) <- c('Ob_down','Ob_up','Ad_down','Ad_up')

# Test the overlap of both gene groups using a hypergeometric test
mat <- matrix(NA, ncol=length(Gene_groups),nrow=length(Gene_groups_2))
colnames(mat) <- names(Gene_groups)
rownames(mat) <- names(Gene_groups_2)

for (i in 1:length(Gene_groups_2)){
  for (k in 1:length(Gene_groups)){
    tmp_all <- All[All %in% TERT$SYMBOL]
    tmp_i <- Gene_groups_2[[i]][Gene_groups_2[[i]] %in% tmp_all]
    tmp_k <- Gene_groups[[k]][Gene_groups[[k]] %in% tmp_all]
    x <- length(tmp_k[tmp_k %in% tmp_i])
    m <- length(tmp_k)
    n <- length(tmp_all[!tmp_all %in% tmp_k])
    y <- length(tmp_i)
    mat[i,k] <- phyper(x,m,n,y,lower.tail = F)
  }
}
# transform to log scale
mat <- -log10(mat)

# Show enrichment in a heatmap
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(mat),length=51))
heatmap.2(mat,Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')
rm(x,m,n,y,i,k,mat_col,mat_col_breaks,tmp_i,tmp_k,tmp_all,mat,markers, Result,Gene_groups,Gene_groups_2, All )


### Figure 2M Relation of cluster specific genes to adipocyte and osteoblast differentiation
Result <- readRDS("MarkerGeneList.Rds")
# Alignment with expression data (https://doi.org/10.1038/s41588-019-0359-1), PCA plot to show that gene groups demarcate differentiation samples
TERT <- read.delim("GeneExpression_TERT_ATERT.txt",h=T)
library(dplyr)
library(scales)

# Combine markers of single RNA-seq in one data frame and split enriched markers in a list
markers <- bind_rows(Result, .id = "column_label")

Gene_groups <- list()
Gene_groups[[1]] <- markers[markers$column_label== "Cl_0" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[2]] <- markers[markers$column_label== "Cl_1" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[3]] <- markers[markers$column_label== "Cl_2" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[4]] <- markers[markers$column_label== "Cl_3" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[5]] <- markers[markers$column_label== "Cl_4" & markers$Marker =="Enriched",'Gene' ]

names(Gene_groups) <- c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4')

# extract variance stablized expression data from TERT cells, day 7 of differentiation and undifferentiated
object <- TERT[,c(65,69,73,76,80,84,87,91,95)]
rownames(object) = TERT$SYMBOL

# Make a grouping object to define Timepoint, Replicate, and color for the data points
intgroup <- c('Timepoint', 'Replicate')
intgroup.df <- data.frame("Sample"=colnames(object),
                          "Timepoint"=rep(c('7dOb','Msc','7dAd'),3),
                          "Replicate"=c(rep('a',3),rep('b',3),rep('c',3)),
                          "color" = rep(c('blue','grey','red'),3))

# Run PCA analysis in a loop over the cluster specific genes and plot results
par(mfrow=c(3,2),pty="s")
for (i in 1:5){
  # Do the PCA
  pca <- prcomp(t(object[Gene_groups[[i]][Gene_groups[[i]] %in% rownames(object)],]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  # Check that the grouping objects matches
  if (!all(intgroup %in% names(intgroup.df))) {
    stop("the argument 'intgroup' should specify columns of intgroup.df")
  }
  rownames(intgroup.df) <- intgroup.df$Sample
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    intgroup.df[[intgroup]]
  }
  # Combine grouping variables and PCA results
  dat <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4], group=group, intgroup.df, name=colnames(object))
  attr(dat, "percentVar") <- percentVar[1:4]
  
  # plot the data
  plot(dat$PC1,dat$PC2,col=dat$color, pch=16, main=names(Gene_groups)[i])
}
rm(Result, TERT, markers, Gene_groups, object, intgroup, intgroup.df,pca,percentVar,group,dat,i)


### Figure 2N Relation of cluster specific genes to adipocyte and osteoblast differentiation
Result <- readRDS("MarkerGeneList.Rds")
# Alignment with expression data (https://doi.org/10.1038/s41588-019-0359-1), quantify changes at day 7 of Ob and Ad differentiation 
TERT <- read.delim("GeneExpression_TERT_ATERT.txt",h=T)
library(dplyr)
library(scales)

# Combine markers of single RNA-seq in one data frame and split enriched markers in a list
markers <- bind_rows(Result, .id = "column_label")

Gene_groups <- list()
Gene_groups[[1]] <- markers[markers$column_label== "Cl_0" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[2]] <- markers[markers$column_label== "Cl_1" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[3]] <- markers[markers$column_label== "Cl_2" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[4]] <- markers[markers$column_label== "Cl_3" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[5]] <- markers[markers$column_label== "Cl_4" & markers$Marker =="Enriched",'Gene' ]

names(Gene_groups) <- c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4')

# Sample 100 random genes for statistic comparison (roughly similar size to the cluster specific genes)
tmp <- sample(TERT$SYMBOL,100)

# Osteogenic differentiation day 7
boxplot(
  TERT[TERT$SYMBOL %in% Gene_groups[[1]],120],
  TERT[TERT$SYMBOL %in% Gene_groups[[2]],120],
  TERT[TERT$SYMBOL %in% Gene_groups[[3]],120],
  TERT[TERT$SYMBOL %in% Gene_groups[[4]],120],
  TERT[TERT$SYMBOL %in% Gene_groups[[5]],120], 
  TERT[TERT$SYMBOL %in% tmp,120],outline=F, col=c(hue_pal()(5),'grey')
)
abline(h=0)
# The p-values will be dependent on the random selection of the 100 genes and are therefore not precisely reproducible
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[1]],120],TERT[TERT$SYMBOL %in% tmp,120])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[2]],120],TERT[TERT$SYMBOL %in% tmp,120])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[3]],120],TERT[TERT$SYMBOL %in% tmp,120])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[4]],120],TERT[TERT$SYMBOL %in% tmp,120])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[5]],120],TERT[TERT$SYMBOL %in% tmp,120])

# Adipogenic differentiation day 7
boxplot(
  TERT[TERT$SYMBOL %in% Gene_groups[[1]],127],
  TERT[TERT$SYMBOL %in% Gene_groups[[2]],127],
  TERT[TERT$SYMBOL %in% Gene_groups[[3]],127],
  TERT[TERT$SYMBOL %in% Gene_groups[[4]],127],
  TERT[TERT$SYMBOL %in% Gene_groups[[5]],127], 
  TERT[TERT$SYMBOL %in% tmp,127],outline=F, col=c(hue_pal()(5),'grey')
)
abline(h=0)
# The p-values will be dependent on the random selection of the 100 genes and are therefore not precisely reproducible
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[1]],127],TERT[TERT$SYMBOL %in% tmp,127])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[2]],127],TERT[TERT$SYMBOL %in% tmp,127])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[3]],127],TERT[TERT$SYMBOL %in% tmp,127])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[4]],127],TERT[TERT$SYMBOL %in% tmp,127])
wilcox.test(TERT[TERT$SYMBOL %in% Gene_groups[[5]],127],TERT[TERT$SYMBOL %in% tmp,127])

rm(Result, TERT, markers, Gene_groups, tmp)


