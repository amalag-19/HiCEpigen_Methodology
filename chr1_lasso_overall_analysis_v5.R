#########################################################################################################
## Building linear model with constraints on first 500 window with TAD locations and epigenetic variables as predictors
#########################################################################################################
## Loading required libraries
library(glmnet)
library(dplyr)
library(penalized)
library(biglm)
library(ff)
library(ffbase)
library(data.table)

#########################################################################################################
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_mat_preprocessed/chr1_mat.RData")

chr_mat_full<-chr_mat

window<-c(1,500)
chr_mat<-chr_mat_full[window[1]:window[2],window[1]:window[2]]

#########################################################################################################
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
chr_pairs_list<-pairs_list_generator(start_ID = window[1],end_ID = window[2],bandwidth_parameter = 500)

## Initializing the dataframe for the chromosome IDs and distance
df_chr_IDs_dist<-do.call(what = rbind,args = chr_pairs_list)
df_chr_IDs_dist$dist<-abs(df_chr_IDs_dist$row_ID-df_chr_IDs_dist$col_ID)
nrow(df_chr_IDs_dist)

#########################################################################################################
## Randomly sample from 5% of data
df_chr_IDs_dist_subset<-df_chr_IDs_dist

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

intensity<-chr_intensity_extractR(chr_mat = chr_mat_full,chr_IDs_mat = chr_IDs_mat)
#intensity_full<-chr_intensity_extractR(chr_mat = chr_mat,chr_IDs_mat = chr_IDs_mat_full)

#########################################################################################################
df_TAD<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/TAD/dp_raw_pen0.1_newestchr1.txt", header = FALSE, sep = "\t")
str(df_TAD)
TADs<-unique(sort(as.numeric(c(df_TAD[-1,1],df_TAD[-1,2]))))
TAD_locations<-TADs[-c(1,length(TADs))]
str(TAD_locations)

TAD_locations<-TAD_locations[which((TAD_locations>=window[1])&(TAD_locations<=window[2]))]

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

design_mat_TAD<-design_mat_generator(interaction_mat = chr_IDs_mat,TAD_locations = TAD_locations)
head(design_mat_TAD)

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

#####################################################
design_mat_epigenetic<-matrix(NA_real_,nrow(chr_IDs_mat),2*ncol(df_epigen_modified))
for (k in 1:nrow(chr_IDs_mat)){
  design_mat_epigenetic[k,1:35]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,1],])
  design_mat_epigenetic[k,36:70]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,2],])
  if((k%%1000)==0){
    print(k)
  }
}

#####################################################
df_chr<-data.frame("Intensity"=intensity,"Distance"=-df_chr_IDs_dist_subset$dist,-design_mat_TAD,design_mat_epigenetic)

lm_fit<-penalized(response = log(Intensity+1), penalized = as.matrix(df_chr[,-1]) , lambda1=0, lambda2=0, positive = c(rep(T,ncol(design_mat_TAD)+1),rep(F,ncol(design_mat_epigenetic))), data=df_chr)

#lm_fit<-penalized(response = log(Intensity+1), penalized = as.matrix(df_chr[,-1]) , lambda1=0, lambda2=0, positive = c(rep(T,ncol(design_mat_TAD)+1),rep(F,ncol(design_mat_epigenetic))), data=df_chr,maxiter = 1000)

## Extracting all coefficients
coeff_all<-as.numeric(coef(lm_fit,which="all"))
length(coeff_all)

## Number of TAD coefficients
length(TAD_locations)

## Number of non-zero epigenetic coefficients
length(which(coeff_all[(length(coeff_all)):(length(coeff_all)-69)]!=0))

## Number of non-zero TAD coefficients
length(which(coeff_all[(3:(length(TAD_locations)+2))]!=0))

## Number of total non-zero coefficients
length(as.numeric(coef(lm_fit)))

predictedValues<-predict(object = lm_fit,penalized=as.matrix(df_chr[,-1]))

#####################################################
## Calculating the correlation
corr_vec_generator<-function(chr_IDs_mat,TAD_locations){
  N<-max(TAD_locations)
  corr_vec<-rep(NA_real_,N)
  for(k in 1:N){
    ids<-which((chr_IDs_mat[,2]-chr_IDs_mat[,1])==k)
    if(length(ids)>10){
      intensity_val_k<-log(intensity[ids]+1)
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
cor(log(intensity+1),as.numeric(predictedValues[,1]))

#########################################################################################################
#########################################################################################################
source('~/Box Sync/PSU/Fall 2018/BioStat_Research/Code/Yu_code.R')
matrix_generator_lm<-function(start_ID=1,end_ID=500,bandwidth=500,chr_mat,TAD_locations,lm_obj){
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
  chr_pairs_list<-pairs_list_generator(start_ID = start_ID,end_ID = end_ID,bandwidth_parameter = bandwidth)
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
  design_mat_TAD<-design_mat_generator(interaction_mat = chr_IDs_mat,TAD_locations = TAD_locations)
  
  #####################################################
  design_mat_epigenetic<-matrix(NA_real_,nrow(chr_IDs_mat),2*ncol(df_epigen_modified))
  for (k in 1:nrow(chr_IDs_mat)){
    design_mat_epigenetic[k,1:35]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,1],])
    design_mat_epigenetic[k,36:70]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,2],])
    if((k%%1000)==0){
      print(k)
    }
  }
  
  #####################################################
  df_chr<-data.frame("Distance"=-df_chr_IDs_dist$dist,-design_mat_TAD,design_mat_epigenetic)
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

#undebug(matrix_generator_lm)
matrix_output<-matrix_generator_lm(start_ID=window[1], end_ID= window[2], chr_mat = chr_mat_full,TAD_locations = TAD_locations,lm_obj = lm_fit)
par(mfrow=c(1,2))

predicted_matrix<-matrix_output[[2]]
common_scale<-matrix_output[[3]]

#plotLoop(x = matrix_output[[1]])
#plotLoop(x = matrix_output[[2]])
plotLoop(x = matrix_output[[1]],myrange = matrix_output[[3]])
plotLoop(x = matrix_output[[2]],myrange = matrix_output[[3]])

