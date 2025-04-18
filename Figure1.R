### Figure 1
# The  files are provided at the open science framework https://osf.io/s8nfb/ 

### Figure 1C showing differentiation potential of primary cells
library(ggplot2)
tbl_diff <- read.delim("PrimaryDifferentiation.txt")

#subset to samples were submitted for scRNA-seq
tbl_diff <- tbl_diff[tbl_diff$scRNAseq ==1,]

#Scatterplot with color gradient
ggplot() + 
  geom_point(data=tbl_diff, aes(x=ALP, y=ORO, colour=AZR)) + 
  #theme(legend.position="none") + 
  scale_colour_gradient(low="grey", high="black", na.value="red") 

#Get statistics
#ALP vs ORO
summary(lm(tbl_diff$ALP~tbl_diff$ORO), pch=16)
# adjusted R-squared is -0.04014, p-value is 0.7882
#AZR vs ORO
summary(lm(tbl_diff$AZR~tbl_diff$ORO), pch=16)
# adjusted R-squared is -0.04398, p-value is 0.8615
#ALP vs AZR
summary(lm(tbl_diff$ALP~tbl_diff$AZR), pch=16)
# adjusted R-squared is 0.662, p-value is 8.101e-07

rm(tbl_diff)


### Figure 1D, stats on the scRNA-seq data, cells per Donor
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)
library(dplyr)

# Extract stats on detected genes and transcripts
mat <- data.frame("Donor" = scobject$Dataset,"Genes" = scobject$detected,"Transcripts" = scobject$sum)

# Plot number of cells as barplot
barplot(table(mat$Donor)[order(table(mat$Donor))])

rm(mat, scobject)



### Figure 1E, stats on the scRNA-seq data, genes and transcripts per Donor
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)
library(dplyr)

# Extract stats on detected genes and transcripts
mat <- data.frame("Donor" = scobject$Dataset,"Genes" = scobject$detected,"Transcripts" = scobject$sum)

# Compute mean values by Donor
mat_summary <- mat %>%
  group_by(Donor) %>%
  summarise(
    Mean_Genes = mean(Genes, na.rm = TRUE),
    Mean_Transcripts = mean(Transcripts, na.rm = TRUE),
    SD_Genes = sd(Genes, na.rm = TRUE),
    SD_Transcripts = sd(Transcripts, na.rm = TRUE)
  )
mat_summary <- data.frame(mat_summary)

# Plot Genes against Transcripts with confidence intervalls using SD
plot(0,0,pch="",xlab="Detected Genes",ylab="Detected Transcripts", xlim=c(0,max(mat_summary$Mean_Genes + mat_summary$SD_Genes)), ylim=c(0,max(mat_summary$Mean_Transcripts + mat_summary$SD_Transcripts)))

arrows(
  x0 = mat_summary$Mean_Genes - mat_summary$SD_Genes,
  x1 = mat_summary$Mean_Genes + mat_summary$SD_Genes,
  y0 = mat_summary$Mean_Transcripts,
  y1 = mat_summary$Mean_Transcripts,
  angle = 90, code = 3, length = 0.025, col = "lightgrey"
)

arrows(
  x0 = mat_summary$Mean_Genes,
  x1 = mat_summary$Mean_Genes,
  y0 = mat_summary$Mean_Transcripts - mat_summary$SD_Transcripts,
  y1 = mat_summary$Mean_Transcripts + mat_summary$SD_Transcripts,
  angle = 90, code = 3, length = 0.025, col = "lightgrey"
)

points(mat_summary$Mean_Genes, mat_summary$Mean_Transcripts)

rm(scobject,mat_summary,mat)



### Figure 1E, stats on the scRNA-seq data, Detected genes, mitochondial gene content, and doublet score on UMAP plot
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)
library(ggplot2)
library(patchwork)

# FeaturePlot of quality measures
FeaturePlot(scobject, "detected", raster=F,pt.size = 0.01) + theme(aspect.ratio = 1)
FeaturePlot(scobject, "DoubletScore", raster=F,pt.size = 0.01) + theme(aspect.ratio = 1)
FeaturePlot(scobject, "subsets_Mito_percent", raster=F,pt.size = 0.01) + theme(aspect.ratio = 1)

rm(scobject)
