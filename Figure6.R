### Figure 6
# The  files are provided at the open science framework https://osf.io/s8nfb/ 

### Figure 6D
# Osteogenic potential of cells sorted into high and low expression of ITGA11, CD151, CD73 
Data <- read.delim("Sorting_MSC.txt",h=T)
Data$condition <- factor(Data$condition, levels=c('low','high'))

# selecting ITGA11 and ALP
Data_tmp <- Data[Data$assay =="ALP" & complete.cases(Data$ITGA11),]

boxplot(ITGA11~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$ITGA11)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$ITGA11)
}

t.test(Data_tmp[Data_tmp$condition == "low",'ITGA11'],
       Data_tmp[Data_tmp$condition == "high",'ITGA11'], paired = T)
# paired T-Test p = 0.02759

# selecting ITGA11 and AZR
Data_tmp <- Data[Data$assay =="AZR" & complete.cases(Data$ITGA11),]
boxplot(ITGA11~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$ITGA11)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$ITGA11)
}

t.test(Data_tmp[Data_tmp$condition == "low",'ITGA11'],
       Data_tmp[Data_tmp$condition == "high",'ITGA11'], paired = T)
# paired T-Test p = 0.01167

# selecting CD73 and ALP
Data_tmp <- Data[Data$assay =="ALP" & complete.cases(Data$CD73),]

boxplot(CD73~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$CD73)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$CD73)
}

t.test(Data_tmp[Data_tmp$condition == "low",'CD73'],
       Data_tmp[Data_tmp$condition == "high",'CD73'], paired = T)
# paired T-Test p = 0.02064

# selecting CD73 and AZR
Data_tmp <- Data[Data$assay =="AZR" & complete.cases(Data$CD73),]
boxplot(CD73~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$CD73)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$CD73)
}

t.test(Data_tmp[Data_tmp$condition == "low",'CD73'],
       Data_tmp[Data_tmp$condition == "high",'CD73'], paired = T)
# paired T-Test p = 0.02536

# selecting CD151 and ALP
Data_tmp <- Data[Data$assay =="ALP" & complete.cases(Data$CD151),]

boxplot(CD151~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$CD151)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$CD151)
}

t.test(Data_tmp[Data_tmp$condition == "low",'CD151'],
       Data_tmp[Data_tmp$condition == "high",'CD151'], paired = T)
# paired T-Test p = 0.3495

# selecting CD151 and AZR
Data_tmp <- Data[Data$assay =="AZR" & complete.cases(Data$CD151),]
boxplot(CD151~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$CD151)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$CD151)
}

t.test(Data_tmp[Data_tmp$condition == "low",'CD151'],
       Data_tmp[Data_tmp$condition == "high",'CD151'], paired = T)
# paired T-Test p = 0.6016
rm(tmp, Data_tmp, Data, i)



### Figure 6F
# Adipogenic potential of cells sorted into high and low expression of ITGA11, CD151, CD73 
Data <- read.delim("Sorting_MSC.txt",h=T)
Data$condition <- factor(Data$condition, levels=c('low','high'))

# selecting ITGA11 and ORO
Data_tmp <- Data[Data$assay =="ORO" & complete.cases(Data$ITGA11),]

boxplot(ITGA11~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$ITGA11)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$ITGA11)
}

t.test(Data_tmp[Data_tmp$condition == "low",'ITGA11'],
       Data_tmp[Data_tmp$condition == "high",'ITGA11'], paired = T)
# paired T-Test p = 0.06287

# selecting CD73 and ORO
Data_tmp <- Data[Data$assay =="ORO" & complete.cases(Data$CD73),]
boxplot(CD73~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$CD73)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$CD73)
}

t.test(Data_tmp[Data_tmp$condition == "low",'CD73'],
       Data_tmp[Data_tmp$condition == "high",'CD73'], paired = T)
# paired T-Test p = 0.5251

# selecting CD151 and ORO
Data_tmp <- Data[Data$assay =="ORO" & complete.cases(Data$CD151),]

boxplot(CD151~condition,Data_tmp, col="white")
points(as.numeric(Data_tmp$condition),Data_tmp$CD151)     
for (i in unique(Data_tmp$Replicate)){
  tmp <- Data_tmp[Data_tmp$Replicate == i,]
  lines(tmp$condition,tmp$CD151)
}

t.test(Data_tmp[Data_tmp$condition == "low",'CD151'],
       Data_tmp[Data_tmp$condition == "high",'CD151'], paired = T)
# paired T-Test p = 0.6305
rm(tmp, Data_tmp, Data, i)



### Figure 6G
# Osteogenic potential in FBS vs. HPL supplemented stromal cultures (Stanford)
Data <- read.delim("FBS_HPL_MSC.txt",h=T)

# collect all paired samples for HPL, FBS and Alizarin Red quantification
k=1
for (i in unique(Data$ID)){
  tmp <- Data[Data$ID == i,]
  if(nrow(tmp[complete.cases(tmp$AZR_Quantification),]) == 2){
   if(k == 1){
     Data_tmp <- Data[Data$ID == i,c('ID','Media','AZR_Quantification')]
   } else{
     Data_tmp <- rbind(Data_tmp,Data[Data$ID == i,c('ID','Media','AZR_Quantification')])
   } 
  } else{
  }
  k <- k+1
}
# Define factor level for media supplement
Data_tmp$Media <- factor(Data_tmp$Media, levels=c('FBS','HPL'))

# scatter plot of ALZ over serum supplementation
boxplot(AZR_Quantification~ Media,Data_tmp, col="white")
points(as.numeric(Data_tmp$Media),Data_tmp$AZR_Quantification)     
axis(1,at=c(1,2), labels = c('FBS','HPL'))
for (i in unique(Data_tmp$ID)){
  tmp <- Data_tmp[Data_tmp$ID == i,]
  lines(tmp$Media,tmp$AZR_Quantification)
}

t.test(Data_tmp[Data_tmp$Media == "FBS",'AZR_Quantification'],
       Data_tmp[Data_tmp$Media == "HPL",'AZR_Quantification'], paired = T)
# paired T-Test p = 0.09776
rm(i,k,Data_tmp, tmp)



### Figure 6I
# Surface marker expression (percentage positive cells) in FBS vs. HPL supplemented stromal cultures (Stanford)
Data <- read.delim("FBS_HPL_MSC.txt",h=T)

# collect all paired samples for HPL, FBS and ITGA11
k=1
for (i in unique(Data$ID)){
  tmp <- Data[Data$ID == i,]
  if(nrow(tmp[complete.cases(tmp$ITGA11_percentage),]) == 2){
    if(k == 1){
      Data_tmp <- Data[Data$ID == i,c('ID','Media','ITGA11_percentage')]
    } else{
      Data_tmp <- rbind(Data_tmp,Data[Data$ID == i,c('ID','Media','ITGA11_percentage')])
    } 
  } else{
  }
  k <- k+1
}
# Define factor level for media supplement
Data_tmp$Media <- factor(Data_tmp$Media, levels=c('FBS','HPL'))

# scatter plot of ITGA11 over serum supplementation
boxplot(ITGA11_percentage~ Media,Data_tmp, col="white")
points(as.numeric(Data_tmp$Media),Data_tmp$ITGA11_percentage)     
for (i in unique(Data_tmp$ID)){
  tmp <- Data_tmp[Data_tmp$ID == i,]
  lines(tmp$Media,tmp$ITGA11_percentage)
}

t.test(Data_tmp[Data_tmp$Media == "FBS",'ITGA11_percentage'],
       Data_tmp[Data_tmp$Media == "HPL",'ITGA11_percentage'], paired = T)
# paired T-Test p = 0.03515
  

# collect all paired samples for HPL, FBS and CD73
k=1
for (i in unique(Data$ID)){
  tmp <- Data[Data$ID == i,]
  if(nrow(tmp[complete.cases(tmp$CD73_percentage),]) == 2){
    if(k == 1){
      Data_tmp <- Data[Data$ID == i,c('ID','Media','CD73_percentage')]
    } else{
      Data_tmp <- rbind(Data_tmp,Data[Data$ID == i,c('ID','Media','CD73_percentage')])
    } 
  } else{
  }
  k <- k+1
}
# Define factor level for media supplement
Data_tmp$Media <- factor(Data_tmp$Media, levels=c('FBS','HPL'))

# scatter plot of ITGA11 over serum supplementation
boxplot(CD73_percentage~ Media,Data_tmp, col="white")
points(as.numeric(Data_tmp$Media),Data_tmp$CD73_percentage)     
for (i in unique(Data_tmp$ID)){
  tmp <- Data_tmp[Data_tmp$ID == i,]
  lines(tmp$Media,tmp$CD73_percentage)
}

t.test(Data_tmp[Data_tmp$Media == "FBS",'CD73_percentage'],
       Data_tmp[Data_tmp$Media == "HPL",'CD73_percentage'], paired = T)
# paired T-Test p = 0.4303

# collect all paired samples for HPL, FBS and CD151
k=1
for (i in unique(Data$ID)){
  tmp <- Data[Data$ID == i,]
  if(nrow(tmp[complete.cases(tmp$CD151_percentage),]) == 2){
    if(k == 1){
      Data_tmp <- Data[Data$ID == i,c('ID','Media','CD151_percentage')]
    } else{
      Data_tmp <- rbind(Data_tmp,Data[Data$ID == i,c('ID','Media','CD151_percentage')])
    } 
  } else{
  }
  k <- k+1
}
# Define factor level for media supplement
Data_tmp$Media <- factor(Data_tmp$Media, levels=c('FBS','HPL'))

# scatter plot of ITGA11 over serum supplementation
boxplot(CD151_percentage~ Media,Data_tmp, col="white")
points(as.numeric(Data_tmp$Media),Data_tmp$CD151_percentage)     
for (i in unique(Data_tmp$ID)){
  tmp <- Data_tmp[Data_tmp$ID == i,]
  lines(tmp$Media,tmp$CD151_percentage)
}

t.test(Data_tmp[Data_tmp$Media == "FBS",'CD151_percentage'],
       Data_tmp[Data_tmp$Media == "HPL",'CD151_percentage'], paired = T)
# paired T-Test p = 0.3358
rm(i,k,Data_tmp, tmp)	 




### Figure 6J
# Surface marker expression (MFI) in FBS vs. HPL supplemented stromal cultures (Stanford)
Data <- read.delim("FBS_HPL_MSC.txt",h=T)

# collect all paired samples for HPL, FBS and ITGA11
k=1
for (i in unique(Data$ID)){
  tmp <- Data[Data$ID == i,]
  if(nrow(tmp[complete.cases(tmp$ITGA11_MFI),]) == 2){
    if(k == 1){
      Data_tmp <- Data[Data$ID == i,c('ID','Media','ITGA11_MFI')]
    } else{
      Data_tmp <- rbind(Data_tmp,Data[Data$ID == i,c('ID','Media','ITGA11_MFI')])
    } 
  } else{
  }
  k <- k+1
}
# Define factor level for media supplement
Data_tmp$Media <- factor(Data_tmp$Media, levels=c('FBS','HPL'))

# scatter plot of ITGA11 over serum supplementation
boxplot(ITGA11_MFI~ Media,Data_tmp, col="white")
points(as.numeric(Data_tmp$Media),Data_tmp$ITGA11_MFI)     
for (i in unique(Data_tmp$ID)){
  tmp <- Data_tmp[Data_tmp$ID == i,]
  lines(tmp$Media,tmp$ITGA11_MFI)
}

t.test(Data_tmp[Data_tmp$Media == "FBS",'ITGA11_MFI'],
       Data_tmp[Data_tmp$Media == "HPL",'ITGA11_MFI'], paired = T)
# paired T-Test p = 0.03158
  

# collect all paired samples for HPL, FBS and CD73
k=1
for (i in unique(Data$ID)){
  tmp <- Data[Data$ID == i,]
  if(nrow(tmp[complete.cases(tmp$CD73_MFI),]) == 2){
    if(k == 1){
      Data_tmp <- Data[Data$ID == i,c('ID','Media','CD73_MFI')]
    } else{
      Data_tmp <- rbind(Data_tmp,Data[Data$ID == i,c('ID','Media','CD73_MFI')])
    } 
  } else{
  }
  k <- k+1
}
# Define factor level for media supplement
Data_tmp$Media <- factor(Data_tmp$Media, levels=c('FBS','HPL'))

# scatter plot of ITGA11 over serum supplementation
boxplot(CD73_MFI~ Media,Data_tmp, col="white")
points(as.numeric(Data_tmp$Media),Data_tmp$CD73_MFI)     
for (i in unique(Data_tmp$ID)){
  tmp <- Data_tmp[Data_tmp$ID == i,]
  lines(tmp$Media,tmp$CD73_MFI)
}

t.test(Data_tmp[Data_tmp$Media == "FBS",'CD73_MFI'],
       Data_tmp[Data_tmp$Media == "HPL",'CD73_MFI'], paired = T)
# paired T-Test p = 0.007304

# collect all paired samples for HPL, FBS and CD151
k=1
for (i in unique(Data$ID)){
  tmp <- Data[Data$ID == i,]
  if(nrow(tmp[complete.cases(tmp$CD151_MFI),]) == 2){
    if(k == 1){
      Data_tmp <- Data[Data$ID == i,c('ID','Media','CD151_MFI')]
    } else{
      Data_tmp <- rbind(Data_tmp,Data[Data$ID == i,c('ID','Media','CD151_MFI')])
    } 
  } else{
  }
  k <- k+1
}
# Define factor level for media supplement
Data_tmp$Media <- factor(Data_tmp$Media, levels=c('FBS','HPL'))

# scatter plot of ITGA11 over serum supplementation
boxplot(CD151_MFI~ Media,Data_tmp, col="white")
points(as.numeric(Data_tmp$Media),Data_tmp$CD151_MFI)     
for (i in unique(Data_tmp$ID)){
  tmp <- Data_tmp[Data_tmp$ID == i,]
  lines(tmp$Media,tmp$CD151_MFI)
}

t.test(Data_tmp[Data_tmp$Media == "FBS",'CD151_MFI'],
       Data_tmp[Data_tmp$Media == "HPL",'CD151_MFI'], paired = T)
# paired T-Test p = 0.07663
rm(i,k,Data_tmp, tmp)	 