### Figure 4
# The  files are provided at the open science framework https://osf.io/s8nfb/ 

### Figure 4B
# Cluster specific expression levels of putative surface proteins, based on a mass spectrometry based surface protein atalas https://wlab.ethz.ch/cspa/#news
scobject <- qs::qread("SeuratObject.qs")
Result <- readRDS("MarkerGeneList.Rds")
Surface <- read.delim("Human_Surface_Proteom_cspa.txt",h=T)
library(Seurat)
library(dplyr)

# extract marker genes into a data frame
markers <- bind_rows(Result, .id = "column_label")

# keep those genes that are putative surface markers
markers_surface <- markers[markers$Gene %in% Surface$Gene_symbol,]

# reduce to enriched markers
markers_surface <- markers_surface[markers_surface$Marker=="Enriched",]

# Show Heatmap
DoHeatmap(scobject, features=markers_surface$Gene,assay = "originalexp")
rm(scobject,markers,Result,markers_surface,Surface)



### Figure 4C
# Dot plot of chosen examples for follow up
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)

# Dot Plot of examples
DotPlot(scobject, features=c("CD151", "NT5E", "ITGA11"))
rm(scobject)



### Figure 4D
# Correlation of mRNA expression level of ITGA11, NT5E (encoding CD73), and CD151 with osteogenic differentiation potnetial
tbl <- read.delim("qPCR_data.txt",h=T)
# mRNA expression levels of ITGA11,CD151, and CD73 in MSCs that were part of the clinical and cell morphology characterisation from previous in vitro studies (high content imaging: 10.1002/sctm.19-0171, clinical signature 10.1186/s13287-021-02338-1 )
Features <- read.delim("Justyna_Data.txt",h=T)
rownames(Features) <- paste("P",Features$Patient, sep="")
rownames(tbl) <- tbl$Donor

# Select differentiation assays
Features <- Features[,c(44,85,87)]

# cbind with feature data, seperated into clinical and cellular data
common <-  rownames(tbl)[rownames(tbl) %in% rownames(Features)]
tbl <- tbl[common,]
Features <- Features[common,]

# Run Pearson's correlation and linear regression between gene expression and differentiation features
cor <- data.frame(matrix(NA,ncol=3,nrow=ncol(Features)))
colnames(cor) <- colnames(tbl)[2:4]
rownames(cor) <- colnames(Features)

sig <- data.frame(matrix(NA,ncol=3,nrow=ncol(Features)))
colnames(sig) <- colnames(tbl)[2:4]
rownames(sig) <- colnames(Features)

for (i in 1:3) {
  for(k in 1:ncol(Features)){
    cor[k,i] <- cor(tbl[,i+1],Features[,k], use="pairwise.complete.obs")
    sig[k,i] <- summary(lm(Features[,k]~tbl[,i+1]))$coefficients[2,4]
  }
}

# show scatter plots for osteogenic differentiation (Alizarin Red)

# Alizarin Red ITGA11
plot(tbl$ITGA11 ,Features_cellular$OsteoDiff.I)
abline(lm(Features_cellular$OsteoDiff.I ~ tbl$ITGA11), lty=2)
summary(lm(Features_cellular$OsteoDiff.I ~ tbl$ITGA11))
# adjusted R-squared is 0.1979, p-value is 0.004888

plot(tbl$CD73 ,Features_cellular$OsteoDiff.I)
abline(lm(Features_cellular$OsteoDiff.I ~ tbl$CD73), lty=2)
summary(lm(Features_cellular$OsteoDiff.I ~ tbl$CD73))
# adjusted R-squared is 0.1121, p-value is 0.02987

plot(tbl$CD151 ,Features_cellular$OsteoDiff.I)
abline(lm(Features_cellular$OsteoDiff.I ~ tbl$CD151), lty=2)
summary(lm(Features_cellular$OsteoDiff.I ~ tbl$CD151))
# adjusted R-squared is 0.01288, p-value is 0.2404
rm(tbl, Features,i,k,cor, sig,common)


### Figure 4G
# Osteogenic differentiation capacity based on ITGA11, CD73, and CD151 protein expression
tbl_diff <- read.delim("PrimaryDifferentiation.txt")

# Alizarin Red ITGA11
plot(tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11'],tbl_diff[tbl_diff$ITGA11 ==1,'AZR'])
abline(lm(tbl_diff[tbl_diff$ITGA11 ==1,'AZR'] ~tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11']), lty=2)
summary(lm(tbl_diff[tbl_diff$ITGA11 ==1,'AZR'] ~tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11']))
# adjusted R-squared is 0.6007, p-value is 0.00862

# Alizarin Red CD73
plot(tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73'],tbl_diff[tbl_diff$CD73 ==1,'AZR'])
abline(lm(tbl_diff[tbl_diff$CD73 ==1,'AZR'] ~tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD73 ==1,'AZR'] ~tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73']))
# adjusted R-squared is 0.516, p-value is 0.01764

# Alizarin Red CD151
plot(tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151'],tbl_diff[tbl_diff$CD151 ==1,'AZR'])
abline(lm(tbl_diff[tbl_diff$CD151 ==1,'AZR'] ~tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD151 ==1,'AZR'] ~tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151']))
# adjusted R-squared is 0.225, p-value is 0.1111



# ALP activity ITGA11
plot(tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11'],tbl_diff[tbl_diff$ITGA11 ==1,'ALP'])
abline(lm(tbl_diff[tbl_diff$ITGA11 ==1,'ALP'] ~tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11']), lty=2)
summary(lm(tbl_diff[tbl_diff$ITGA11 ==1,'ALP'] ~tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11']))
# adjusted R-squared is 0.7036, p-value is 0.002896

# ALP activity CD73
plot(tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73'],tbl_diff[tbl_diff$CD73 ==1,'ALP'])
abline(lm(tbl_diff[tbl_diff$CD73 ==1,'ALP'] ~tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD73 ==1,'ALP'] ~tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73']))
# adjusted R-squared is 0.5229, p-value is 0.01672

# ALP activity CD151
plot(tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151'],tbl_diff[tbl_diff$CD151 ==1,'ALP'])
abline(lm(tbl_diff[tbl_diff$CD151 ==1,'ALP'] ~tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD151 ==1,'ALP'] ~tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151']))
# adjusted R-squared is 0.1191, p-value is 0.1923
rm(tbl_diff)



### Figure 4I
# Adipogenic differentiation capacity based on ITGA11, CD73, and CD151 protein expression
tbl_diff <- read.delim("PrimaryDifferentiation.txt")

plot(tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11'],tbl_diff[tbl_diff$ITGA11 ==1,'ORO'])
abline(lm(tbl_diff[tbl_diff$ITGA11 ==1,'ORO'] ~tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11']), lty=2)
summary(lm(tbl_diff[tbl_diff$ITGA11 ==1,'ORO'] ~tbl_diff[tbl_diff$ITGA11 ==1,'MFI_ITGA11']))
# adjusted R-squared is -0.1266, p-value is 0.6604

plot(tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73'],tbl_diff[tbl_diff$CD73 ==1,'ORO'])
abline(lm(tbl_diff[tbl_diff$CD73 ==1,'ORO'] ~tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD73 ==1,'ORO'] ~tbl_diff[tbl_diff$CD73 ==1,'MFI_CD73']))
# adjusted R-squared is 0.3857, p-value is 0.05924

plot(tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151'],tbl_diff[tbl_diff$CD151 ==1,'ORO'])
abline(lm(tbl_diff[tbl_diff$CD151 ==1,'ORO'] ~tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD151 ==1,'ORO'] ~tbl_diff[tbl_diff$CD151 ==1,'MFI_CD151']))
# adjusted R-squared is -0.07715, p-value is 0.5343
rm(tbl_diff)



### Figure 4K
# Ectopic bone formation capacity based on ITGA11, CD73, and CD151 protein expression
tbl_diff <- read.delim("PrimaryDifferentiation.txt")

plot(tbl_diff[tbl_diff$ITGA11 ==1 & tbl_diff$Implant ==1 ,'MFI_ITGA11'],tbl_diff[tbl_diff$ITGA11 ==1& tbl_diff$Implant ==1 ,'Quant_Implant'])
abline(lm(tbl_diff[tbl_diff$ITGA11 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'] ~tbl_diff[tbl_diff$ITGA11 ==1 & tbl_diff$Implant ==1 ,'MFI_ITGA11']), lty=2)
summary(lm(tbl_diff[tbl_diff$ITGA11 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'] ~tbl_diff[tbl_diff$ITGA11 ==1 & tbl_diff$Implant ==1 ,'MFI_ITGA11']))
# adjusted R-squared is 0.8294, p-value is 0.02022

plot(tbl_diff[tbl_diff$CD73 ==1 & tbl_diff$Implant ==1 ,'MFI_CD73'],tbl_diff[tbl_diff$CD73 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'])
abline(lm(tbl_diff[tbl_diff$CD73 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'] ~tbl_diff[tbl_diff$CD73 ==1 & tbl_diff$Implant ==1 ,'MFI_CD73']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD73 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'] ~tbl_diff[tbl_diff$CD73 ==1 & tbl_diff$Implant ==1 ,'MFI_CD73']))
# adjusted R-squared is 0.521, p-value is 0.1037

plot(tbl_diff[tbl_diff$CD151 ==1 & tbl_diff$Implant ==1 ,'MFI_CD151'],tbl_diff[tbl_diff$CD151 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'])
abline(lm(tbl_diff[tbl_diff$CD151 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'] ~tbl_diff[tbl_diff$CD151 ==1 & tbl_diff$Implant ==1 ,'MFI_CD151']), lty=2)
summary(lm(tbl_diff[tbl_diff$CD151 ==1 & tbl_diff$Implant ==1 ,'Quant_Implant'] ~tbl_diff[tbl_diff$CD151 ==1 & tbl_diff$Implant ==1 ,'MFI_CD151']))
# adjusted R-squared is 0.2212, p-value is 0.3066
rm(tbl_diff)



