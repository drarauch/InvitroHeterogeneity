### Figure 3
# The  files are provided at the open science framework https://osf.io/s8nfb/ 

### Figure 3B show single Donor in the data set
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)
library(ggplot2)

# Dimplot grouped by Donor
DimPlot(scobject, raster=F, group.by="Dataset",pt.size = 0.01) + theme(aspect.ratio = 1)

rm(scobject)



### Figure 3C Highlight each single Donor on a UMAP plot
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)
library(ggplot2)
library(patchwork)

# Make vectors for each Donor and highlight them in a Dimplot with the other cells being grey
Idents(object = scobject) <- "Dataset"
Donor <- names(table(scobject$Dataset))
Donor_list <- list()
x <- 1
for (k in Donor){
	Donor_list[[x]] <- WhichCells(scobject, idents = k)
	names(Donor_list)[[x]] <- k
	x <- x+1
}

# Set up for multiplotting
p1 <- DimPlot(scobject, cells.highlight= Donor_list[[1]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p2 <- DimPlot(scobject, cells.highlight= Donor_list[[2]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p3 <- DimPlot(scobject, cells.highlight= Donor_list[[3]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p4 <- DimPlot(scobject, cells.highlight= Donor_list[[4]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p5 <- DimPlot(scobject, cells.highlight= Donor_list[[5]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p6 <- DimPlot(scobject, cells.highlight= Donor_list[[6]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p7 <- DimPlot(scobject, cells.highlight= Donor_list[[7]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p8 <- DimPlot(scobject, cells.highlight= Donor_list[[8]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p9 <- DimPlot(scobject, cells.highlight= Donor_list[[9]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p10 <- DimPlot(scobject, cells.highlight= Donor_list[[10]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p11 <- DimPlot(scobject, cells.highlight= Donor_list[[11]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p12 <- DimPlot(scobject, cells.highlight= Donor_list[[12]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p13 <- DimPlot(scobject, cells.highlight= Donor_list[[13]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p14 <- DimPlot(scobject, cells.highlight= Donor_list[[14]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p15 <- DimPlot(scobject, cells.highlight= Donor_list[[15]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p16 <- DimPlot(scobject, cells.highlight= Donor_list[[16]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p17 <- DimPlot(scobject, cells.highlight= Donor_list[[17]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p18 <- DimPlot(scobject, cells.highlight= Donor_list[[18]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p19 <- DimPlot(scobject, cells.highlight= Donor_list[[19]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p20 <- DimPlot(scobject, cells.highlight= Donor_list[[20]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p21 <- DimPlot(scobject, cells.highlight= Donor_list[[21]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p22 <- DimPlot(scobject, cells.highlight= Donor_list[[22]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p23 <- DimPlot(scobject, cells.highlight= Donor_list[[23]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p24 <- DimPlot(scobject, cells.highlight= Donor_list[[24]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p25 <- DimPlot(scobject, cells.highlight= Donor_list[[25]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)
p26 <- DimPlot(scobject, cells.highlight= Donor_list[[26]], cols.highlight = c("darkblue"), cols= "grey95") + theme(aspect.ratio = 1)

combined_plot1 <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)
combined_plot2 <- p7 + p8 + p9 + p10 + p11 + p12 + plot_layout(ncol = 2)
combined_plot3 <- p13 + p14 + p15 + p16 + p17 + p18 + plot_layout(ncol = 2)
combined_plot4 <- p19 + p20 + p21 + p22 + p23 + p24 + plot_layout(ncol = 2)
combined_plot5 <- p25 + p26 + p1 + p1 + p1 + p1 + plot_layout(ncol = 2)

#Export as pdf
pdf("Donor_in_UMAP.pdf",paper="A4")
print(combined_plot1)
print(combined_plot2)
print(combined_plot3)
print(combined_plot4)
print(combined_plot5)
dev.off()

rm(Donor,Donor_list,x,k,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26)



### Figure 3D grouping of donors based on cluster abundancy
scobject <- qs::qread("SeuratObject.qs")
library(Seurat)
library(scales)

#Cluster abundance across donor samples
tbl <- table(scobject$seurat_clusters, scobject$Dataset)
for (i in 1:ncol(tbl)){
  tbl[,i] <- tbl[,i] / as.numeric(table(scobject$Dataset)[i])
}

# cluster the donor samples
hr <- hclust(dist(t(tbl)))

# plot the dendrogram
plot(hr)

# Use clustering order to show cluster abundancy across donors
barplot(tbl[,hr$order],legend=F, col=hue_pal()(nrow(tbl)), las=2, horiz=T)

rm(tbl,hr,scobject)



### Figure 3E Differentiation ability of donors grouped due to their subpopulation composition 
scobject <- qs::qread("SeuratObject.qs")
tbl_diff <- read.delim("PrimaryDifferentiation.txt")
library(Seurat)

# Extract cluster abundancy per donor
tbl <- table(scobject$seurat_clusters, scobject$Dataset)
for (i in 1:ncol(tbl)){
  tbl[,i] <- tbl[,i] / as.numeric(table(scobject$Dataset)[i])
}
tbl <- t(tbl)
tbl <- as.data.frame.matrix(tbl)

# Cluster the donor samples and order according to dendrogram
hr <- hclust(dist(tbl))
tbl <- tbl[hr$order,]
colnames(tbl) <- paste("Cl_",colnames(tbl),sep="")

# cbind with differentiation data 
rownames(tbl_diff) <- tbl_diff$Donor
tbl_diff <- tbl_diff[rev(rownames(tbl)),]
tbl <- cbind(tbl,tbl_diff[,2:4])

# Histogram plots
barplot(tbl$AZR, names=rownames(tbl), ylab="AlizarinRed", las=2, col='grey')
barplot(tbl$ALP, names=rownames(tbl), ylab="ALP", las=2, col='grey')
barplot(tbl$ORO, names=rownames(tbl), ylab="OilRedO", las=2, col='grey')

rm(tbl,hr,tbl_diff,scobject,i)

# Group donors based on hirarchical clustering (abitrary, here we used 5 clusters, called a,b,c,d,e)

a <- names(cutree(hr,5)[cutree(hr,5)==4])
b <- names(cutree(hr,5)[cutree(hr,5)==2])
c <- names(cutree(hr,5)[cutree(hr,5)==3])
d <- names(cutree(hr,5)[cutree(hr,5)==5])
e <- names(cutree(hr,5)[cutree(hr,5)==1])

# Plot the individual donors (keeping the order of the dendrogram) and show means and std for the subgroups a, b, c, d, and e
plot(1:5,c(mean(tbl[a,'AZR'],na.rm=T),mean(tbl[b,'AZR'],na.rm=T),mean(tbl[c,'AZR'],na.rm=T),mean(tbl[d,'AZR'],na.rm=T),mean(tbl[e,'AZR'],na.rm=T)), ylim=c(0,max(tbl[,'AZR'],na.rm=T)))
lines(c(1,1),c(mean(tbl[a,'AZR'],na.rm=T)-sd(tbl[a,'AZR'],na.rm=T),mean(tbl[a,'AZR'],na.rm=T)+sd(tbl[a,'AZR'],na.rm=T)))
lines(c(2,2),c(mean(tbl[b,'AZR'],na.rm=T)-sd(tbl[b,'AZR'],na.rm=T),mean(tbl[b,'AZR'],na.rm=T)+sd(tbl[b,'AZR'],na.rm=T)))
lines(c(3,3),c(mean(tbl[c,'AZR'],na.rm=T)-sd(tbl[c,'AZR'],na.rm=T),mean(tbl[c,'AZR'],na.rm=T)+sd(tbl[c,'AZR'],na.rm=T)))
lines(c(4,4),c(mean(tbl[d,'AZR'],na.rm=T)-sd(tbl[d,'AZR'],na.rm=T),mean(tbl[d,'AZR'],na.rm=T)+sd(tbl[d,'AZR'],na.rm=T)))
lines(c(5,5),c(mean(tbl[e,'AZR'],na.rm=T)-sd(tbl[e,'AZR'],na.rm=T),mean(tbl[e,'AZR'],na.rm=T)+sd(tbl[e,'AZR'],na.rm=T)))

plot(1:5,c(mean(tbl[a,'ALP'],na.rm=T),mean(tbl[b,'ALP'],na.rm=T),mean(tbl[c,'ALP'],na.rm=T),mean(tbl[d,'ALP'],na.rm=T),mean(tbl[e,'ALP'],na.rm=T)), ylim=c(0,max(tbl[,'ALP'],na.rm=T)))
lines(c(1,1),c(mean(tbl[a,'ALP'],na.rm=T)-sd(tbl[a,'ALP'],na.rm=T),mean(tbl[a,'ALP'],na.rm=T)+sd(tbl[a,'ALP'],na.rm=T)))
lines(c(2,2),c(mean(tbl[b,'ALP'],na.rm=T)-sd(tbl[b,'ALP'],na.rm=T),mean(tbl[b,'ALP'],na.rm=T)+sd(tbl[b,'ALP'],na.rm=T)))
lines(c(3,3),c(mean(tbl[c,'ALP'],na.rm=T)-sd(tbl[c,'ALP'],na.rm=T),mean(tbl[c,'ALP'],na.rm=T)+sd(tbl[c,'ALP'],na.rm=T)))
lines(c(4,4),c(mean(tbl[d,'ALP'],na.rm=T)-sd(tbl[d,'ALP'],na.rm=T),mean(tbl[d,'ALP'],na.rm=T)+sd(tbl[d,'ALP'],na.rm=T)))
lines(c(5,5),c(mean(tbl[e,'ALP'],na.rm=T)-sd(tbl[e,'ALP'],na.rm=T),mean(tbl[e,'ALP'],na.rm=T)+sd(tbl[e,'ALP'],na.rm=T)))

plot(1:5,c(mean(tbl[a,'ORO'],na.rm=T),mean(tbl[b,'ORO'],na.rm=T),mean(tbl[c,'ORO'],na.rm=T),mean(tbl[d,'ORO'],na.rm=T),mean(tbl[e,'ORO'],na.rm=T)), ylim=c(0,max(tbl[,'ORO'],na.rm=T)))
lines(c(1,1),c(mean(tbl[a,'ORO'],na.rm=T)-sd(tbl[a,'ORO'],na.rm=T),mean(tbl[a,'ORO'],na.rm=T)+sd(tbl[a,'ORO'],na.rm=T)))
lines(c(2,2),c(mean(tbl[b,'ORO'],na.rm=T)-sd(tbl[b,'ORO'],na.rm=T),mean(tbl[b,'ORO'],na.rm=T)+sd(tbl[b,'ORO'],na.rm=T)))
lines(c(3,3),c(mean(tbl[c,'ORO'],na.rm=T)-sd(tbl[c,'ORO'],na.rm=T),mean(tbl[c,'ORO'],na.rm=T)+sd(tbl[c,'ORO'],na.rm=T)))
lines(c(4,4),c(mean(tbl[d,'ORO'],na.rm=T)-sd(tbl[d,'ORO'],na.rm=T),mean(tbl[d,'ORO'],na.rm=T)+sd(tbl[d,'ORO'],na.rm=T)))
lines(c(5,5),c(mean(tbl[e,'ORO'],na.rm=T)-sd(tbl[e,'ORO'],na.rm=T),mean(tbl[e,'ORO'],na.rm=T)+sd(tbl[e,'ORO'],na.rm=T)))



### Figure 2F Check if differentiation potential is linked to cell abundancy across the clusters, aligning donors according to clustering in 2I or direct Pearson's correlation
scobject <- qs::qread("SeuratObject.qs")
tbl_diff <- read.delim("PrimaryDifferentiation.txt")
library(Seurat)
library(gplots)
library(fields)

# Extract cluster abundancy per donor
tbl <- table(scobject$seurat_clusters, scobject$Dataset)
for (i in 1:ncol(tbl)){
  tbl[,i] <- tbl[,i] / as.numeric(table(scobject$Dataset)[i])
}
tbl <- t(tbl)
tbl <- as.data.frame.matrix(tbl)
colnames(tbl) <- paste("Cl_",colnames(tbl),sep="")

# cbind with differentiation data 
rownames(tbl_diff) <- tbl_diff$Donor
tbl_diff <- tbl_diff[rownames(tbl),]
tbl <- cbind(tbl,tbl_diff[,2:4])

# Pearson's correlation between cluster abundance and differentiation potential
cor <- data.frame(matrix(NA,ncol=5,nrow=3))
colnames(cor) <- colnames(tbl)[1:5]
rownames(cor) <- colnames(tbl)[6:8]

for (i in 1:5) {
	for(k in 6:8){
		cor[k-5,i] <- cor(tbl[,i],tbl[,k], use="pairwise.complete.obs")
	}
}

# Heatmap
mat_col <- designer.colors(n=50, col=c('blue','white','red'))
mat_col_breaks <- c(seq(-1,1, length=51))
heatmap.2(as.matrix(cor), Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')

# Testing significance of correlation using linear regression models
sig <- data.frame(matrix(NA,ncol=5,nrow=3))
colnames(sig) <- colnames(tbl)[1:5]
rownames(sig) <- colnames(tbl)[6:8]

for (i in 1:5) {
	for(k in 6:8){
		sig[k-5,i] <- summary(lm(tbl[,k]~tbl[,i]))$coefficients[2,4]
	}
}

# Cluster 0 abundancy correlates significantly with ALP and Alizarin Red, while Cl_1 and Cl_2 correlate with OilRed O
rm(i,k,sig,mat_col, mat_col_breaks,cor,scobject,tbl, tbl_diff)



### Figure 3G Scatter plots of ALP and AlizarinRed against Cl0 abundancy as this reached significance in 2J
scobject <- qs::qread("SeuratObject.qs")
tbl_diff <- read.delim("PrimaryDifferentiation.txt")
library(Seurat)
library(gplots)
library(fields)

# Extract cluster abundancy per donor
tbl <- table(scobject$seurat_clusters, scobject$Dataset)
for (i in 1:ncol(tbl)){
  tbl[,i] <- tbl[,i] / as.numeric(table(scobject$Dataset)[i])
}
tbl <- t(tbl)
tbl <- as.data.frame.matrix(tbl)
colnames(tbl) <- paste("Cl_",colnames(tbl),sep="")

# cbind with differentiation data 
rownames(tbl_diff) <- tbl_diff$Donor
tbl_diff <- tbl_diff[rownames(tbl),]
tbl <- cbind(tbl,tbl_diff[,2:4])

# Plot Alizarin Red against Cl0
plot(tbl$Cl_0,tbl$AZR)
abline(lm(tbl$AZR ~tbl$Cl_0), lty=2)
summary(lm(tbl$AZR ~tbl$Cl_0))
# adjusted R-squared is 0.4062, p-value is 0.0003653
# Plot ALP against Cl0
plot(tbl$Cl_0,tbl$ALP)
abline(lm(tbl$ALP ~tbl$Cl_0), lty=2)
summary(lm(tbl$ALP ~tbl$Cl_0))
# adjusted R-squared is 0.4264, p-value is 0.000179
# Plot ORO against Cl1
plot(tbl$Cl_1,tbl$ORO)
abline(lm(tbl$ORO ~tbl$Cl_1), lty=2)
summary(lm(tbl$ORO ~tbl$Cl_1))
# adjusted R-squared is 0.1719, p-value is 0.02018
# Plot ORO against Cl2
plot(tbl$Cl_2,tbl$ORO)
abline(lm(tbl$ORO ~tbl$Cl_2), lty=2)
summary(lm(tbl$ORO ~tbl$Cl_2))
# adjusted R-squared is 0.2058, p-value is 0.01155

rm(tbl, sig, cor, i,k,mat_col, mat_col_breaks,tbl_diff,scobject)



### Figure 3H Pearson's correlation between cluster abundance and clinical and cellular features from previous in vitro studies (high content imaging: 10.1002/sctm.19-0171, clinical signature 10.1186/s13287-021-02338-1 )
scobject <- qs::qread("SeuratObject.qs")
Features <- read.delim("Justyna_Data.txt")
library(Seurat)
library(gplots)
library(fields)

rownames(Features) <- paste("P",Features$Patient, sep="")

# Extract cluster abundancy per donor
tbl <- table(scobject$seurat_clusters, scobject$Dataset)
for (i in 1:ncol(tbl)){
  tbl[,i] <- tbl[,i] / as.numeric(table(scobject$Dataset)[i])
}
tbl <- t(tbl)
tbl <- as.data.frame.matrix(tbl)
colnames(tbl) <- paste("Cl_",colnames(tbl),sep="")

# cbind with feature data, seperated into clinical and cellular data
common <-  rownames(tbl)[rownames(tbl) %in% rownames(Features)]
tbl <- tbl[common,]
Features <- Features[common,]
Features_clinical <- Features[,c(2:27)]
Features_cellular <- Features[,c(44,85,87,28:36,45:83)]

# Run Pearson's correlation and linear regression between cluster abundance and cellular features
cor_cellular <- data.frame(matrix(NA,ncol=5,nrow=ncol(Features_cellular)))
colnames(cor_cellular) <- colnames(tbl)[1:5]
rownames(cor_cellular) <- colnames(Features_cellular)

sig_cellular <- data.frame(matrix(NA,ncol=5,nrow=ncol(Features_cellular)))
colnames(sig_cellular) <- colnames(tbl)[1:5]
rownames(sig_cellular) <- colnames(Features_cellular)

for (i in 1:5) {
	for(k in 1:ncol(Features_cellular)){
		cor_cellular[k,i] <- cor(tbl[,i],Features_cellular[,k], use="pairwise.complete.obs")
		sig_cellular[k,i] <- summary(lm(Features_cellular[,k]~tbl[,i]))$coefficients[2,4]
	}
}

# Only keep pearson's correaltion for the significant assocaitions
cor_cellular[sig_cellular > 0.05] <- 0

# Heatmap (move differentition read out to the top)
mat_col <- designer.colors(n=50, col=c('blue','white','red'))
mat_col_breaks <- c(seq(-1,1, length=51))
heatmap.2(as.matrix(cor_cellular), Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')


# Run Pearson's correlation and linear regression between cluster abundance and clinical features
cor_clinical <- data.frame(matrix(NA,ncol=5,nrow=ncol(Features_clinical)))
colnames(cor_clinical) <- colnames(tbl)[1:5]
rownames(cor_clinical) <- colnames(Features_clinical)

sig_clinical <- data.frame(matrix(NA,ncol=5,nrow=ncol(Features_clinical)))
colnames(sig_clinical) <- colnames(tbl)[1:5]
rownames(sig_clinical) <- colnames(Features_clinical)

for (i in 1:5) {
	for(k in 1:ncol(Features_clinical)){
		cor_clinical[k,i] <- cor(tbl[,i],Features_clinical[,k], use="pairwise.complete.obs")
		sig_clinical[k,i] <- summary(lm(Features_clinical[,k]~tbl[,i]))$coefficients[2,4]
	}
}

# Only keep pearson's correaltion for the significant assocaitions
cor_clinical[sig_clinical > 0.05] <- 0

# Heatmap
mat_col <- designer.colors(n=50, col=c('blue','white','red'))
mat_col_breaks <- c(seq(-1,1, length=51))
heatmap.2(as.matrix(cor_clinical), Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')

rm(tbl, sig_clinical,sig_cellular, Features, Features_cellular, Features_clinical,scobject,cor_clinical,cor_cellular, i,k,mat_col, mat_col_breaks,tbl_diff,scobject)

