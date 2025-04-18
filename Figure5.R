### Figure 5
# The  files are provided at the open science framework https://osf.io/s8nfb/ 

### Figure 5B
# Flow cytometry based Surface marker expression on MSC (Standford) cultured in HPL 
Data <- read.delim("SurfaceMarker_MSC.txt",h=T)

# boxplot
vec <- c("CD73_percentage","ITGA11_percentage","CD151_percentage","PDPN_percentage","TIE2_percentage","CD31_percentage","CD146_percentage","CD164_percentage")
boxplot(Data[,vec])
# add points
for (i in 1:length(vec)){
  points(rep(i, nrow(Data)),Data[,vec[i]])
}
rm(Data,vec)

### Figure 5C
# Correlation of percentage positive cells for SSC markers and CD151, CD73, and ITGA11 on MSC (Standford) cultured in HPL 
library(gplots)
library(fields)
Data <- read.delim("SurfaceMarker_MSC.txt",h=T)

# Select "Percentage positive cells"
vec <- c("CD73_percentage","ITGA11_percentage","CD151_percentage","PDPN_percentage","TIE2_percentage","CD146_percentage","CD164_percentage")

# Run correlation and test for significance using linear regression model
y_cor <- cor(Data[,vec])
y_sig <- data.frame(matrix(NA,ncol=length(vec),nrow=length(vec)))
colnames(y_sig) <- colnames(y_cor)
rownames(y_sig) <- rownames(y_cor)

for (i in 1:length(vec)) {
  for(k in 1:length(vec)){
    y_sig[i,k] <- summary(lm(Data[,vec[k]]~Data[,vec[i]]))$coefficients[2,4]
  }
}

# Heatmap
mat_col <- designer.colors(n=50, col=c('blue','white','red'))
mat_col_breaks <- c(seq(-1,1, length=51))
heatmap.2(as.matrix(y_cor), Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')

rm(Data,vec,mat_col.mat_col_breaks,y_sig, y_cor)



### Figure 5D scatter plot of significant examples in 4C
# Correlation of percentage positive cells for SSC markers and CD151, CD73, and ITGA11 on MSC (Standford) cultured in HPL 
Data <- read.delim("SurfaceMarker_MSC.txt",h=T)

par(mfrow=c(1,2), pty="s")
plot(Data[,c("TIE2_percentage")],Data[,c("ITGA11_percentage")])
abline(lm(Data[,c("ITGA11_percentage")]~Data[,c("TIE2_percentage")]))
summary(lm(Data[,c("ITGA11_percentage")]~Data[,c("TIE2_percentage")]))
# adjusted R-squared is 0.634, p-value is 0.0003986

plot(Data[,c("PDPN_percentage")],Data[,c("CD73_percentage")])
abline(lm(Data[,c("CD73_percentage")]~Data[,c("PDPN_percentage")]))
summary(lm(Data[,c("CD73_percentage")]~Data[,c("PDPN_percentage")]))
# adjusted R-squared is 0.7843, p-value is 1.544e-05

rm(Data)



### Figure 5E
# Correlation of MFI for SSC markers and CD151, CD73, and ITGA11 on MSC (Standford) cultured in HPL 
library(gplots)
library(fields)
Data <- read.delim("SurfaceMarker_MSC.txt",h=T)

# Select "Percentage positive cells"
vec <- c("CD73_MFI","ITGA11_MFI","CD151_MFI","PDPN_MFI","TIE2_MFI","CD146_MFI","CD164_MFI")

# Run correlation and test for significance using linear regression model
y_cor <- cor(Data[,vec])
y_sig <- data.frame(matrix(NA,ncol=length(vec),nrow=length(vec)))
colnames(y_sig) <- colnames(y_cor)
rownames(y_sig) <- rownames(y_cor)

for (i in 1:length(vec)) {
  for(k in 1:length(vec)){
    y_sig[i,k] <- summary(lm(Data[,vec[k]]~Data[,vec[i]]))$coefficients[2,4]
  }
}

# Heatmap
mat_col <- designer.colors(n=50, col=c('blue','white','red'))
mat_col_breaks <- c(seq(-1,1, length=51))
heatmap.2(as.matrix(y_cor), Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')

rm(Data,vec,mat_col.mat_col_breaks,y_sig, y_cor)



### Figure 5H
# Correlation of ITGA11 on MSC and ITGA11 on gated SSC (Standford) cultured in HPL 
Data <- read.delim("SurfaceMarker_MSC.txt",h=T)

par(mfrow=c(1,2), pty="s")
# Scatter plot of ITGA11 percentage positive MSC versus ITGA11 percentage positive SSC
plot(Data[,c("ITGA11_percentage")],Data[,c("SSC_ITGA11_percentage")])
abline(lm(Data[,c("SSC_ITGA11_percentage")]~Data[,c("ITGA11_percentage")]))
summary(lm(Data[,c("SSC_ITGA11_percentage")]~Data[,c("ITGA11_percentage")]))
# adjusted R-squared is 0.6829, p-value is 0.0001641

rm(Data)



### Figure 5I
# Correlation of ITGA11 on MSC and gated SSC (Standford) cultured in HPL 
Data <- read.delim("SurfaceMarker_MSC.txt",h=T)

par(mfrow=c(1,2), pty="s")
# Scatter plot of ITGA11 percentage positive MSC versus percentage of SSC
plot(Data[,c("ITGA11_percentage")],Data[,c("SSC_percentage")])
abline(lm(Data[,c("SSC_percentage")]~Data[,c("ITGA11_percentage")]))
summary(lm(Data[,c("SSC_percentage")]~Data[,c("ITGA11_percentage")]))
# adjusted R-squared is 0.3434, p-value is 0.02075

rm(Data)