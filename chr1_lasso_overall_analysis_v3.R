#########################################################################################################
## Loading required libraries
library(glmnet)
library(dplyr)

#########################################################################################################
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_mat_preprocessed/chr1_mat.RData")

# Defining a function to get the paired row IDs and column corrsponding the start and end IDs and bandwidth parameter (taking into account the symmetric matrix)
pairs_list_generator<-function(start_ID,end_ID,bandwidth_parameter){
  pairs_list<-list()
  for(i in start_ID:end_ID){
    col_ID_seq<-(i-(bandwidth_parameter-1)):(i+(bandwidth_parameter-1))
    col_ID_seq_positive<-col_ID_seq[which((col_ID_seq>=max(i,start_ID))&(col_ID_seq<=end_ID))]
    pairs_list[[i]]<-data.frame("row_ID"=rep(i,length(col_ID_seq_positive)),"col_ID"=col_ID_seq_positive)
  }
  return(pairs_list)
}

## Number of bins in the chromosome matrix
N<-dim(chr_mat)[1]
## Getting the chr pairs list
chr_pairs_list<-pairs_list_generator(start_ID = 1,end_ID = N,bandwidth_parameter = 500)

## Initializing the dataframe for the chromosome IDs and distance
df_chr_IDs_dist<-do.call(what = rbind,args = chr_pairs_list)
df_chr_IDs_dist$dist<-abs(df_chr_IDs_dist$row_ID-df_chr_IDs_dist$col_ID)
nrow(df_chr_IDs_dist)
#########################################################################################################
## Randomly sample from 5% of data
randomly_sampled_row_IDs<-sample(x = 1:nrow(df_chr_IDs_dist),size = 500000,replace = F)

df_chr_IDs_dist_subset<-df_chr_IDs_dist[randomly_sampled_row_IDs,]

#########################################################################################################
## Randomly sample 5% of data but stratified over sub-diagonals

summary(df_chr_IDs_dist[,2]-df_chr_IDs_dist[,1])

## Initializing the total row IDs vector for random sampling across sub-diagonals
row_IDs<-c()
for(k in 0:499){
  row_IDs_k<-which(df_chr_IDs_dist$dist==k)
  if(length(row_IDs_k)>=20){
    row_IDs_k_subset<-row_IDs_k[sample(x = 1:length(row_IDs_k),size = 0.05*length(row_IDs_k),replace = F)]
    row_IDs<-c(row_IDs,row_IDs_k_subset)
  }
}

df_chr_IDs_dist_subset<-df_chr_IDs_dist[row_IDs,]

summary(df_chr_IDs_dist_subset[,1])

#####################################################
## Getting the IDs matrix
chr_IDs_mat<-as.matrix(df_chr_IDs_dist_subset[,c(1,2)])
attr(chr_IDs_mat,"dimnames")<-NULL

#chr_IDs_mat_full<-as.matrix(df_chr_IDs_dist[,c(1,2)])
#attr(chr_IDs_mat_full,"dimnames")<-NULL

str(chr_IDs_mat)
#str(chr_IDs_mat_full)
#save(chr_IDs_mat,)

#####################################################
chr_intensity_extractR<-function(chr_mat, chr_IDs_mat){
  intensityVec<-rep(NA,nrow(chr_IDs_mat))
  for(i in 1:nrow(chr_IDs_mat)){
    intensityVec[i]<-chr_mat[chr_IDs_mat[i,1],chr_IDs_mat[i,2]]
    if((i%%100000)==0){
      print(i)
    }
  }
  return(intensityVec)
}

intensity<-chr_intensity_extractR(chr_mat = chr_mat,chr_IDs_mat = chr_IDs_mat)
#intensity_full<-chr_intensity_extractR(chr_mat = chr_mat,chr_IDs_mat = chr_IDs_mat_full)

#########################################################################################################
df_TAD<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/TAD/dp_raw_pen0.1_newestchr1.txt", header = FALSE, sep = "\t")
str(df_TAD)
TADs<-unique(sort(as.numeric(c(df_TAD[-1,1],df_TAD[-1,2]))))
TAD_locations<-TADs[-c(1,length(TADs))]
str(TAD_locations)

#save(TAD_locations,file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/TAD_locations.RData")

#####################################################
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

design_mat<-design_mat_generator(interaction_mat = chr_IDs_mat,TAD_locations = TAD_locations)
#design_mat_full<-design_mat_generator(interaction_mat = chr_IDs_mat_full,TAD_locations = TAD_locations)

#####################################################
df_chr<-data.frame("Intensity"=intensity,"Distance"=df_chr_IDs_dist_subset$dist,design_mat)

save(df_chr,file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/df_chr_stratified.RData")

penFactor<-c(0, rep(1,dim(df_chr)[2]-2))

lasso=glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=1,nlambda=100,upper.limits = rep(0,ncol(df_chr)),penalty.factor = penFactor)

save(lasso,file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/lasso/lasso.RData")

cv.lasso<-cv.glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=1,nfolds=10,upper.limits = rep(0,ncol(df_chr)), penalty.factor = penFactor)

## get lambda and best lasso fit
lambda.lasso<-cv.lasso$lambda.1se
print(lambda.lasso)

save(lambda.lasso,file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/lasso/lambda_lasso.RData")

## beta estimates for best lambda
betas.lasso=coef(cv.lasso,s=lambda.lasso)
beta_vec<-rep(0,length(betas.lasso@Dimnames[[1]]))
beta_vec[betas.lasso@i + 1]<-betas.lasso@x
#print(beta_vec[2])
#print(lasso$dev.ratio[which(lasso$lambda == lambda.lasso)])

df_coeff<-data.frame(name = as.character(TAD_locations), coefficient = beta_vec[-c(1,2)],stringsAsFactors = F)
#str(df_coeff)
df_coeff$name<-factor(df_coeff$name, levels = df_coeff[["name"]])

predictedValues<-predict(lasso, s = lambda.lasso, newx = as.matrix(df_chr[,-1]), type = "response")

#####################################################
## Calculating the correlation
corr_vec_generator<-function(chr_IDs_mat,TAD_locations){
  N<-max(TAD_locations)
  corr_vec<-rep(NA_real_,N)
  for(k in 1:N){
    ids<-which((chr_IDs_mat[,2]-chr_IDs_mat[,1])==k)
    if(length(ids)>10){
      intensity_val_k<-intensity[ids]
      intensity_val_k_predicted<-predictedValues[ids]
      corr_vec[k]<-cor(intensity_val_k,intensity_val_k_predicted)
    }
    if((k%%100)==0){
      print(k)
    }
  }
  return(corr_vec)
}

corr_vec<-corr_vec_generator(chr_IDs_mat = chr_IDs_mat, TAD_locations = TAD_locations)

non_NA_ids<-which(!is.na(corr_vec))

corr_vec_modified<-corr_vec[non_NA_ids]
N<-max(TAD_locations)
dist_vec<-(1:N)[non_NA_ids]


plot(x = dist_vec,y = corr_vec_modified)
## Overall correlation
cor(intensity,predictedValues)

#####################################################
#####################################################
## Sourcing plotloop to create to plot the original and predicted matrix
source('~/Box Sync/PSU/Fall 2018/BioStat_Research/Code/Yu_code.R')
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_mat_preprocessed/chr1_mat.RData")

matrix_generator_lasso<-function(start_ID=601,end_ID=800,bandwidth=200,chr_mat,TAD_locations,lasso_obj,lambda_lasso){
  plotRange<-c(start_ID,end_ID)
  # Defining a function to get the paired row IDs and column corrsponding the start and end IDs and bandwidth parameter (taking into account the symmetric matrix)
  pairs_list_generator<-function(start_ID,end_ID,bandwidth_parameter){
    pairs_list<-list()
    for(i in start_ID:end_ID){
      col_ID_seq<-(i-(bandwidth_parameter-1)):(i+(bandwidth_parameter-1))
      col_ID_seq_positive<-col_ID_seq[which((col_ID_seq>=max(i,start_ID))&(col_ID_seq<=end_ID))]
      pairs_list[[i]]<-data.frame("row_ID"=rep(i,length(col_ID_seq_positive)),"col_ID"=col_ID_seq_positive)
    }
    return(pairs_list)
  }
  chr_pairs_list<-pairs_list_generator(start_ID = 601,end_ID = 800,bandwidth_parameter = bandwidth)
  ## Initializing the dataframe for the chromosome IDs and distance
  df_chr_IDs_dist<-do.call(what = rbind,args = chr_pairs_list)
  df_chr_IDs_dist$dist<-abs(df_chr_IDs_dist$row_ID-df_chr_IDs_dist$col_ID)
  #####################################################
  ## Getting the IDs matrix
  chr_IDs_mat<-as.matrix(df_chr_IDs_dist[,c(1,2)])
  attr(chr_IDs_mat,"dimnames")<-NULL
  #####################################################
  chr_intensity_extractR<-function(chr_mat, chr_IDs_mat){
    intensityVec<-rep(NA,nrow(chr_IDs_mat))
    for(i in 1:nrow(chr_IDs_mat)){
      intensityVec[i]<-chr_mat[chr_IDs_mat[i,1],chr_IDs_mat[i,2]]
      if((i%%100000)==0){
        print(i)
      }
    }
    return(intensityVec)
  }
  intensity<-chr_intensity_extractR(chr_mat = chr_mat,chr_IDs_mat = chr_IDs_mat)
  #####################################################
  design_mat_generator<-function(interaction_mat,TAD_locations){
    N<-nrow(interaction_mat)
    design_mat<-matrix(0,nrow = N,ncol = length(TAD_locations))
    for(i in 1:nrow(interaction_mat)){
      for (j in 1:length(TAD_locations)){
        if((TAD_locations[j]>=interaction_mat[i,1])&(TAD_locations[j]<=interaction_mat[i,2])){
          design_mat[i,j]<-1
        }
      }
      # if((i%%1000)==0){
      #   print(i)
      # }
    }
    return(design_mat)
  }
  design_mat<-design_mat_generator(interaction_mat = chr_IDs_mat,TAD_locations = TAD_locations)
  #####################################################
  df_chr<-data.frame("Distance"=df_chr_IDs_dist$dist,design_mat)
  predictedValues<-predict(lasso_obj, s = lambda_lasso, newx = as.matrix(df_chr), type = "response")
  predicted_matrix<-matrix(0,length(unique(chr_IDs_mat[,1])),length(unique(chr_IDs_mat[,2])))
  for(i in 1:nrow(chr_IDs_mat)){
    predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]<-predictedValues[i]
    predicted_matrix[chr_IDs_mat[i,2]-start_ID+1,chr_IDs_mat[i,1]-start_ID+1]<-predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]
  }
  quantile_log_intensity<-quantile(c(log(intensity+1)),prob = c(0.0,0.99))
  quantile_predicted_log_intensity<-quantile(predictedValues,prob = c(0.0,0.99))
  Rstart<-max(quantile_log_intensity[1],quantile_predicted_log_intensity[1])
  Rend<-min(quantile_log_intensity[2],quantile_predicted_log_intensity[2])
  Rrange<-c(Rstart,Rend)
  original_matrix<-log(chr_mat[start_ID:end_ID,start_ID:end_ID]+1)
  return(list(original_matrix,predicted_matrix))
}

load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/TAD_locations.RData")
load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/lasso/lasso.RData")
load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/lasso/lambda_lasso.RData")

debug(matrix_generator_lasso)
matrix_output<-matrix_generator_lasso(chr_mat = chr_mat,TAD_locations = TAD_locations,lasso_obj = lasso,lambda_lasso = lambda.lasso)
plotLoop(x = matrix_output[[1]])
plotLoop(x = matrix_output[[2]])

########################################################################################################
########################################################################################################
########################################################################################################
## Fitting the linear model without epigenetic
load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/df_chr.RData")

df_chr[,-1]<--df_chr[,-1]

library(penalized)
lm_fit<-penalized(response = Intensity, penalized = as.matrix(df_chr[,-1]) , lambda1=0, lambda2=0, positive = rep(T,ncol(df_chr)-1), data=df_chr)
summary(lm_fit)

predictedValues<-predict(object = lm_fit,penalized=as.matrix(df_chr[,-1]))

#####################################################
## Calculating the correlation
corr_vec_generator<-function(chr_IDs_mat,TAD_locations){
  N<-max(TAD_locations)
  corr_vec<-rep(NA_real_,N)
  for(k in 1:N){
    ids<-which((chr_IDs_mat[,2]-chr_IDs_mat[,1])==k)
    if(length(ids)>10){
      intensity_val_k<-intensity[ids]
      intensity_val_k_predicted<-predictedValues[ids]
      corr_vec[k]<-cor(intensity_val_k,intensity_val_k_predicted)
    }
    if((k%%100)==0){
      print(k)
    }
  }
  return(corr_vec)
}

corr_vec<-corr_vec_generator(chr_IDs_mat = chr_IDs_mat, TAD_locations = TAD_locations)
non_NA_ids<-which(!is.na(corr_vec))
corr_vec_modified<-corr_vec[non_NA_ids]
N<-max(TAD_locations)
dist_vec<-(1:N)[non_NA_ids]
par(mfrow=c(1,1))
plot(x = dist_vec,y = corr_vec_modified)

## Overall correlation
cor(intensity,predictedValues)

save(lm_fit,file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/lm/lm_fit_stratified.RData")

########################################################################################################
source('~/Box Sync/PSU/Fall 2018/BioStat_Research/Code/Yu_code.R')
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_mat_preprocessed/chr1_mat.RData")
load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/TAD_locations.RData")
load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/lm/lm_fit.RData")
load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Sep25_Results/lm/lm_fit.RData.RData")

matrix_generator_lm<-function(start_ID=601,end_ID=800,bandwidth=200,chr_mat,TAD_locations,lm_obj){
  plotRange<-c(start_ID,end_ID)
  # Defining a function to get the paired row IDs and column corrsponding the start and end IDs and bandwidth parameter (taking into account the symmetric matrix)
  pairs_list_generator<-function(start_ID,end_ID,bandwidth_parameter){
    pairs_list<-list()
    for(i in start_ID:end_ID){
      col_ID_seq<-(i-(bandwidth_parameter-1)):(i+(bandwidth_parameter-1))
      col_ID_seq_positive<-col_ID_seq[which((col_ID_seq>=max(i,start_ID))&(col_ID_seq<=end_ID))]
      pairs_list[[i]]<-data.frame("row_ID"=rep(i,length(col_ID_seq_positive)),"col_ID"=col_ID_seq_positive)
    }
    return(pairs_list)
  }
  chr_pairs_list<-pairs_list_generator(start_ID = 601,end_ID = 800,bandwidth_parameter = bandwidth)
  ## Initializing the dataframe for the chromosome IDs and distance
  df_chr_IDs_dist<-do.call(what = rbind,args = chr_pairs_list)
  df_chr_IDs_dist$dist<-abs(df_chr_IDs_dist$row_ID-df_chr_IDs_dist$col_ID)
  #####################################################
  ## Getting the IDs matrix
  chr_IDs_mat<-as.matrix(df_chr_IDs_dist[,c(1,2)])
  attr(chr_IDs_mat,"dimnames")<-NULL
  #####################################################
  chr_intensity_extractR<-function(chr_mat, chr_IDs_mat){
    intensityVec<-rep(NA,nrow(chr_IDs_mat))
    for(i in 1:nrow(chr_IDs_mat)){
      intensityVec[i]<-chr_mat[chr_IDs_mat[i,1],chr_IDs_mat[i,2]]
      if((i%%100000)==0){
        print(i)
      }
    }
    return(intensityVec)
  }
  intensity<-chr_intensity_extractR(chr_mat = chr_mat,chr_IDs_mat = chr_IDs_mat)
  #####################################################
  design_mat_generator<-function(interaction_mat,TAD_locations){
    N<-nrow(interaction_mat)
    design_mat<-matrix(0,nrow = N,ncol = length(TAD_locations))
    for(i in 1:nrow(interaction_mat)){
      for (j in 1:length(TAD_locations)){
        if((TAD_locations[j]>=interaction_mat[i,1])&(TAD_locations[j]<=interaction_mat[i,2])){
          design_mat[i,j]<-1
        }
      }
      # if((i%%1000)==0){
      #   print(i)
      # }
    }
    return(design_mat)
  }
  design_mat<-design_mat_generator(interaction_mat = chr_IDs_mat,TAD_locations = TAD_locations)
  #####################################################
  df_chr<-data.frame("Distance"=-df_chr_IDs_dist$dist,-design_mat)
  predictedValues<-as.numeric(predict(object = lm_obj, penalized = as.matrix(df_chr))[,1])
  predicted_matrix<-matrix(0,length(unique(chr_IDs_mat[,1])),length(unique(chr_IDs_mat[,2])))
  for(i in 1:nrow(chr_IDs_mat)){
    predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]<-predictedValues[i]
    predicted_matrix[chr_IDs_mat[i,2]-start_ID+1,chr_IDs_mat[i,1]-start_ID+1]<-predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]
  }
  quantile_log_intensity<-quantile(c(log(intensity+1)),prob = c(0.0,0.99))
  quantile_predicted_log_intensity<-quantile(predictedValues,prob = c(0.0,0.99))
  Rstart<-max(quantile_log_intensity[1],quantile_predicted_log_intensity[1])
  Rend<-min(quantile_log_intensity[2],quantile_predicted_log_intensity[2])
  Rrange<-c(Rstart,Rend)
  original_matrix<-log(chr_mat[start_ID:end_ID,start_ID:end_ID]+1)
  return(list(original_matrix,predicted_matrix,Rrange))
}

debug(matrix_generator_lm)
matrix_output<-matrix_generator_lm(chr_mat = chr_mat,TAD_locations = TAD_locations,lm_obj = lm_fit_2)
par(mfrow=c(1,2))
plotLoop(x = matrix_output[[1]])
plotLoop(x = matrix_output[[2]])



########################################################################################################
predictedValues = predict(lm_fit, newx = as.matrix(df_chr[,-1]), type = "response")


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







