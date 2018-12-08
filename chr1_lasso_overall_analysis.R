#########################################################################################################
## Loading required libraries
library(glmnet)
library(dplyr)
library(data.table)

## Loading the p value matrix from fitHiC output
chr1_Gm12878_fitHiC_output<-read.csv(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_processed_fitHiC/chr1_Gm12878_fitHiC_output/output_subset.csv",header = T)

str(chr1_Gm12878_fitHiC_output)

chr1_Gm12878_fitHiC_output<-chr1_Gm12878_fitHiC_output[,c(2,3,4)]

str(chr1_Gm12878_fitHiC_output)

sample_length<-nrow(chr1_Gm12878_fitHiC_output)

# ids<-which(chr1_Gm12878_fitHiC_output[,6]<0.05)

#sample_length<-length(which(chr1_Gm12878_fitHiC_output[,6]<0.05))

#interaction_mat<-chr1_Gm12878_fitHiC_output[ids,]



#########################################################################################################

df_TAD<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/TAD/dp_raw_pen0.1_newestchr1.txt", header = FALSE, sep = "\t")
str(df_TAD)
TADs<-unique(sort(as.numeric(c(df_TAD[-1,1],df_TAD[-1,2]))))
TAD_locations<-TADs[-c(1,length(TADs))]
str(TAD_locations)

#########################################################################################################
design_mat_generator<-function(interaction_mat,TAD_locations){
  N<-nrow(interaction_mat)
  design_mat<-matrix(0,nrow = N,ncol = length(TAD_locations))
  for(i in 1:nrow(interaction_mat)){
    for (j in 1:length(TAD_locations)){
      if((TAD_locations[j]>=interaction_mat[i,1])&(TAD_locations[j]<=interaction_mat[i,2])){
        design_mat[i,j]<-1
      }
    }
    if((i%%1000)==0){
      print(i)
    }
  }
  return(design_mat)
}

design_mat<-design_mat_generator(interaction_mat = chr1_Gm12878_fitHiC_output,TAD_locations = TAD_locations)

save(design_mat,file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_processed_fitHiC/design_mat.RData")

df_chr<-data.frame("Intensity"=chr1_Gm12878_fitHiC_output[,3],"Distance"=abs(chr1_Gm12878_fitHiC_output[,2]-chr1_Gm12878_fitHiC_output[,1]),design_mat)

penFactor = c(0, rep(1,dim(df_chr)[2]-2))

lasso=glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=1,nlambda=100,upper.limits = rep(0,ncol(df_chr)),penalty.factor = penFactor)

cv.lasso = cv.glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=1,nfolds=10,upper.limits = rep(0,ncol(df_chr)), penalty.factor = penFactor)

## get lambda and best lasso fit
lambda.lasso=cv.lasso$lambda.1se
print(lambda.lasso)

## beta estimates for best lambda
betas.lasso=coef(cv.lasso,s=lambda.lasso)
beta_vec<-rep(0,length(betas.lasso@Dimnames[[1]]))
beta_vec[betas.lasso@i + 1]<-betas.lasso@x
#print(beta_vec[2])
#print(lasso$dev.ratio[which(lasso$lambda == lambda.lasso)])

df_coeff<-data.frame(name = betas.lasso@Dimnames[[1]], coefficient = beta_vec,stringsAsFactors = F)
df_coeff= df_coeff[-c(1,2),]
df_coeff$name = paste0("X",as.character(as.numeric(substring(df_coeff$name, 2))+start_ID - 1))
#str(df_coeff)
df_coeff$name <- factor(df_coeff$name, levels = df_coeff[["name"]])

predictedValues = predict(lasso, s = lambda.lasso, newx = as.matrix(df_chr[,-1]), type = "response")

corr_vec_generator<-function(chr1_Gm12878_fitHiC_output,TAD_locations){
  N<-max(TAD_locations)
  corr_vec<-rep(NA_real_,N)
  for(k in 1:N){
    ids<-which((chr1_Gm12878_fitHiC_output[,2]-chr1_Gm12878_fitHiC_output[,1])==k)
    if(length(ids)>10){
      intensity_val_k<-chr1_Gm12878_fitHiC_output[ids,3]
      intensity_val_k_predicted<-predictedValues[ids]
      corr_vec[k]<-cor(intensity_val_k,intensity_val_k_predicted)
    }
    if((k%%100)==0){
      print(k)
    }
  }
  return(corr_vec)
}

corr_vec<-corr_vec_generator(chr1_Gm12878_fitHiC_output = chr1_Gm12878_fitHiC_output, TAD_locations = TAD_locations)

non_NA_ids<-which(!is.na(corr_vec))

corr_vec_modified<-corr_vec[non_NA_ids]
N<-max(TAD_locations)
dist_vec<-(1:N)[non_NA_ids]


plot(x = dist_vec,y = corr_vec_modified)
## Overall correlation
cor(chr1_Gm12878_fitHiC_output[,3],predictedValues)

#########################################################################################################
## Linear model
lm_fit<-lm(Intensity~., data = df_chr)
summary(lm_fit)

predictedValues = predict(lm_fit, newx = as.matrix(df_chr[,-1]), type = "response")

corr_vec<-corr_vec_generator(chr1_Gm12878_fitHiC_output = chr1_Gm12878_fitHiC_output, TAD_locations = TAD_locations)

non_NA_ids<-which(!is.na(corr_vec))

corr_vec_modified<-corr_vec[non_NA_ids]
N<-max(TAD_locations)
dist_vec<-(1:N)[non_NA_ids]


plot(x = dist_vec,y = corr_vec_modified)

## Overall correlation
cor(chr1_Gm12878_fitHiC_output[,3],predictedValues)

#########################################################################################################
## Epigenetic data preprocessing
df_epigen<-fread(input = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/encodestate_chr1_10kb.txt",data.table = F)

str(df_epigen)

df_epigen[,1]

epigen_mat<-as.matrix(df_epigen)

epigen_mat_normalized<-t(apply(X = epigen_mat,MARGIN = 1,FUN = function(x){
  y<-x/(sum(x)+1e-4)
  return(y)
}))

df_epigen_modified<-data.frame(epigen_mat_normalized[,-1])
str(df_epigen_modified)
names(df_epigen_modified)

#########################################################################################################
## Linear model with epigenetic covariates
design_epigenetic<-matrix(NA_real_,nrow(chr1_Gm12878_fitHiC_output),2*ncol(df_epigen_modified))
for (k in 1:nrow(chr1_Gm12878_fitHiC_output)){
  design_epigenetic[k,1:35]<-as.numeric(df_epigen_modified[chr1_Gm12878_fitHiC_output[k,"fragmentMid1"],])
  design_epigenetic[k,36:70]<-as.numeric(df_epigen_modified[chr1_Gm12878_fitHiC_output[k,"fragmentMid2"],])
  if((k%%100)==0){
    print(k)
  }
}

df_chr_modified<-data.frame(df_chr,design_epigenetic)

str(df_chr_modified)

lm_fit<-lm(Intensity~., data = df_chr_modified)
summary(lm_fit)

predictedValues = predict(lm_fit, newx = as.matrix(df_chr_modified[,-1]), type = "response")

corr_vec<-corr_vec_generator(chr1_Gm12878_fitHiC_output = chr1_Gm12878_fitHiC_output, TAD_locations = TAD_locations)

non_NA_ids<-which(!is.na(corr_vec))

corr_vec_modified<-corr_vec[non_NA_ids]
N<-max(TAD_locations)
dist_vec<-(1:N)[non_NA_ids]


plot(x = dist_vec,y = corr_vec_modified)

## Overall correlation
cor(chr1_Gm12878_fitHiC_output[,3],predictedValues)

#########################################################################################################
## Lasso with epigenetic covariates
penFactor = c(0, rep(1,dim(df_chr_modified)[2]-2))

lasso=glmnet(x=as.matrix(df_chr_modified[,-1]),y=log(df_chr_modified[,1]+1),alpha=1,nlambda=100,upper.limits = rep(0,ncol(df_chr_modified)),penalty.factor = penFactor)

cv.lasso = cv.glmnet(x=as.matrix(df_chr_modified[,-1]),y=log(df_chr_modified[,1]+1),alpha=1,nfolds=10,upper.limits = rep(0,ncol(df_chr_modified)), penalty.factor = penFactor)

## get lambda and best lasso fit
lambda.lasso=cv.lasso$lambda.1se

predictedValues = predict(lasso, s = lambda.lasso, newx = as.matrix(df_chr_modified[,-1]), type = "response")

corr_vec<-corr_vec_generator(chr1_Gm12878_fitHiC_output = chr1_Gm12878_fitHiC_output, TAD_locations = TAD_locations)

non_NA_ids<-which(!is.na(corr_vec))

corr_vec_modified<-corr_vec[non_NA_ids]
N<-max(TAD_locations)
dist_vec<-(1:N)[non_NA_ids]


plot(x = dist_vec,y = corr_vec_modified)

## Overall correlation
cor(chr1_Gm12878_fitHiC_output[,3],predictedValues)

## beta estimates for best lambda
betas.lasso=coef(cv.lasso,s=lambda.lasso)
beta_vec<-rep(0,length(betas.lasso@Dimnames[[1]]))
beta_vec[betas.lasso@i + 1]<-betas.lasso@x
#print(beta_vec[2])
#print(lasso$dev.ratio[which(lasso$lambda == lambda.lasso)])

df_coeff<-data.frame(name = betas.lasso@Dimnames[[1]], coefficient = beta_vec,stringsAsFactors = F)
df_coeff= df_coeff[-c(1,2),]
df_coeff$name = paste0("X",as.character(as.numeric(substring(df_coeff$name, 2))))
#str(df_coeff)
df_coeff$name <- factor(df_coeff$name, levels = df_coeff[["name"]])

length(which(df_coeff[(nrow(df_coeff)-69):nrow(df_coeff),2]!=0))

length(which(df_coeff[1:(nrow(df_coeff)-69),2]!=0))




