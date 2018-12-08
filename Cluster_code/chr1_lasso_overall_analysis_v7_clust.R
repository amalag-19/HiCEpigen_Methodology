#########################################################################################################
## Building linear model with constraints on first 500 window with TAD locations and epigenetic variables as predictors
#########################################################################################################
## Loading required libraries
library(glmnet)
library(dplyr)
library(penalized)
library(data.table)
#library(DAAG) 

#########################################################################################################
file_path<-"/storage/home/a/aua257/work/Bio_Stat/penalized_full_models/"
load(file = paste0(file_path,"data/chr1_mat.RData"))

chr_mat_full<-chr_mat

window<-c(1,nrow(chr_mat_full))
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
set.seed(1)
randomly_sampled_row_IDs<-sample(x = 1:nrow(df_chr_IDs_dist),size = nrow(df_chr_IDs_dist)/100,replace = F) 
head(randomly_sampled_row_IDs)
df_chr_IDs_dist_subset<-df_chr_IDs_dist[sort(randomly_sampled_row_IDs),]

#####################################################
## Getting the IDs matrix
chr_IDs_mat<-as.matrix(df_chr_IDs_dist_subset[,c(1,2)])
attr(chr_IDs_mat,"dimnames")<-NULL

#chr_IDs_mat_full<-as.matrix(df_chr_IDs_dist[,c(1,2)])
#attr(chr_IDs_mat_full,"dimnames")<-NULL

str(chr_IDs_mat)
head(chr_IDs_mat)
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
df_TAD<-read.table(file = paste0(file_path,"data/dp_raw_pen0.1_newestchr1.txt"), header = FALSE, sep = "\t")
str(df_TAD)
TADs<-unique(sort(as.numeric(c(df_TAD[-1,1],df_TAD[-1,2]))))
TAD_locations<-TADs[-c(1,length(TADs))]
str(TAD_locations)
#save(TAD_locations,file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Apps/HiCEpigen/TAD_locations.RData")

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

print("design_mat_TAD is done")

save(design_mat_TAD,file = paste0(file_path,"results/design_mat_TAD.RData"))

#########################################################################################################
## Epigenetic data preprocessing
df_epigen<-fread(input = paste0(file_path,"data/encodestate_chr1_10kb.txt"),data.table = F)

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

#save(df_epigen_modified,file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Apps/HiCEpigen/df_epigen_modified.RData")

#####################################################
design_mat_epigenetic<-matrix(NA_real_,nrow(chr_IDs_mat),2*ncol(df_epigen_modified))
for (k in 1:nrow(chr_IDs_mat)){
  design_mat_epigenetic[k,1:35]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,1],])
  design_mat_epigenetic[k,36:70]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,2],])
  if((k%%1000)==0){
    print(k)
  }
}
print("design_mat_epigenetic is done")

#####################################################
df_chr<-data.frame(intensity,df_chr_IDs_dist_subset$dist,design_mat_TAD,design_mat_epigenetic)

## Processing column names of the dataframe
TAD_numbers<-as.character(TAD_locations) 
TAD_names<-paste("TAD_",TAD_numbers, sep="")
Epi_numbers<-as.character(1:ncol(design_mat_epigenetic)) 
Epi_names<-paste("Epi_",Epi_numbers, sep="")
df_chr_names<-c("Intensity","Distance",TAD_names,Epi_names)
names(df_chr)<-df_chr_names

## Checking vif
#lm_fit<-lm(log(Intensity+1)~.,data = df_chr)
#vif(lm_fit)
#save(lm_fit,file = paste0(file_path,"results/lm_fit.RData"))

## Defining the penalty factor to impose no penalty on distance and TAD locations but L2 penalty on epigenetic variables
penFactor<-c(rep(0,ncol(design_mat_TAD)+1),rep(1,ncol(design_mat_epigenetic)))
#penFactor<-c(rep(1,ncol(design_mat_TAD)+1),rep(0,ncol(design_mat_epigenetic)))

model_fit<-glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=0,nlambda=100,upper.limits = c(rep(0,ncol(design_mat_TAD)+1),rep(Inf,ncol(design_mat_epigenetic))),penalty.factor = penFactor)

cvl_fit<-cv.glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=0,nfolds=10,upper.limits = c(rep(0,ncol(design_mat_TAD)+1),rep(Inf,ncol(design_mat_epigenetic))),penalty.factor = penFactor)

# library(glmnetUtils)
# formula<-log(Intensity + 1) ~ Distance + TAD_126 + TAD_149 + 
#   TAD_156 + TAD_170 + TAD_185 + TAD_207 + TAD_214 + TAD_236 + 
#   TAD_246 + TAD_258 + TAD_370 + TAD_379 + Epi_1 + Epi_2 + Epi_4 + 
#   Epi_5 + Epi_6 + Epi_7 + Epi_8 + Epi_9 + Epi_10 + Epi_12 + 
#   Epi_13 + Epi_14 + Epi_15 + Epi_16 + Epi_17 + Epi_18 + Epi_19 + 
#   Epi_20 + Epi_21 + Epi_22 + Epi_23 + Epi_24 + Epi_25 + Epi_26 + 
#   Epi_27 + Epi_28 + Epi_29 + Epi_30 + Epi_31 + Epi_32 + Epi_33 + 
#   Epi_34 + Epi_35 + Epi_36 + Epi_37 + Epi_39 + Epi_40 + Epi_41 + 
#   Epi_42 + Epi_43 + Epi_44 + Epi_45 + Epi_46 + Epi_47 + Epi_48 + 
#   Epi_49 + Epi_50 + Epi_51 + Epi_52 + Epi_53 + Epi_54 + Epi_55 + 
#   Epi_56 + Epi_57 + Epi_58 + Epi_59 + Epi_60 + Epi_61 + Epi_62 + 
#   Epi_63 + Epi_64 + Epi_65 + Epi_66 + Epi_67 + Epi_68 + Epi_69 + 
#   Epi_70
# model_fit<-glmnet(formula,alpha=0,nlambda=100,upper.limits = c(rep(0,ncol(design_mat_TAD)+1),rep(Inf,ncol(design_mat_epigenetic))),penalty.factor = penFactor,data=df_chr)
# cvl_fit<-cv.glmnet(formula,alpha=0,nfolds=10,upper.limits = c(rep(0,ncol(design_mat_TAD)+1),rep(Inf,ncol(design_mat_epigenetic))),penalty.factor = penFactor,data = df_chr)

## Extracting all coefficients
tuning_param<-cvl_fit$lambda.min
print(tuning_param)

## beta estimates for best lambda
model_coeff_obj<-coef(model_fit,s=tuning_param)
model_coeff_obj

model_coeff_vec<-rep(0,length(model_coeff_obj@Dimnames[[1]]))
model_coeff_vec[model_coeff_obj@i + 1]<-model_coeff_obj@x

df_coeff<-data.frame(name = model_coeff_obj@Dimnames[[1]], coefficient = model_coeff_vec,stringsAsFactors = F)

save(model_fit,file = paste0(file_path,"results/model_fit.RData"))
save(tuning_param,file = paste0(file_path,"results/tuning_param.RData"))

#########################################################################################################
#########################################################################################################
# source('~/Box Sync/PSU/Fall 2018/BioStat_Research/Code/Yu_code.R')
# predict_output_gen<-function(start_ID=1,end_ID=500,bandwidth=500,chr_mat,TAD_locations,model_fit_obj,model_tuning_param=0,df_chr_names){
#   plotRange<-c(start_ID,end_ID)
#   # Defining a function to get the paired row IDs and column corrsponding the start and end IDs and bandwidth parameter (taking into account the symmetric matrix)
#   pairs_list_generator<-function(start_ID,end_ID,bandwidth_parameter){
#     pairs_list<-list()
#     for(i in start_ID:end_ID){
#       col_ID_seq<-(i-(bandwidth_parameter-1)):(i+(bandwidth_parameter-1))
#       col_ID_seq_positive<-col_ID_seq[which((col_ID_seq>=max(i,start_ID))&(col_ID_seq<=end_ID))]
#       pairs_list[[i]]<-data.frame("row_ID"=rep(i,length(col_ID_seq_positive)),"col_ID"=col_ID_seq_positive)
#     }
#     return(pairs_list)
#   }
#   chr_pairs_list<-pairs_list_generator(start_ID = start_ID,end_ID = end_ID,bandwidth_parameter = bandwidth)
#   ## Initializing the dataframe for the chromosome IDs and distance
#   df_chr_IDs_dist<-do.call(what = rbind,args = chr_pairs_list)
#   df_chr_IDs_dist$dist<-abs(df_chr_IDs_dist$row_ID-df_chr_IDs_dist$col_ID)
#   #####################################################
#   ## Getting the IDs matrix
#   chr_IDs_mat<-as.matrix(df_chr_IDs_dist[,c(1,2)])
#   attr(chr_IDs_mat,"dimnames")<-NULL
#   #####################################################
#   chr_intensity_extractR<-function(chr_mat, chr_IDs_mat){
#     intensityVec<-rep(NA,nrow(chr_IDs_mat))
#     for(i in 1:nrow(chr_IDs_mat)){
#       intensityVec[i]<-chr_mat[chr_IDs_mat[i,1],chr_IDs_mat[i,2]]
#       if((i%%100000)==0){
#         print(i)
#       }
#     }
#     return(intensityVec)
#   }
#   intensity<-chr_intensity_extractR(chr_mat = chr_mat,chr_IDs_mat = chr_IDs_mat)
#   #####################################################
#   design_mat_generator<-function(interaction_mat,TAD_locations){
#     N<-nrow(interaction_mat)
#     design_mat<-matrix(0,nrow = N,ncol = length(TAD_locations))
#     for(i in 1:nrow(interaction_mat)){
#       for (j in 1:length(TAD_locations)){
#         if((TAD_locations[j]>=interaction_mat[i,1])&(TAD_locations[j]<=interaction_mat[i,2])){
#           design_mat[i,j]<-1
#         }
#       }
#       # if((i%%1000)==0){
#       #   print(i)
#       # }
#     }
#     return(design_mat)
#   }
#   design_mat_TAD<-design_mat_generator(interaction_mat = chr_IDs_mat,TAD_locations = TAD_locations)
#   
#   #####################################################
#   design_mat_epigenetic<-matrix(NA_real_,nrow(chr_IDs_mat),2*ncol(df_epigen_modified))
#   for (k in 1:nrow(chr_IDs_mat)){
#     design_mat_epigenetic[k,1:35]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,1],])
#     design_mat_epigenetic[k,36:70]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,2],])
#     if((k%%1000)==0){
#       print(k)
#     }
#   }
#   
#   #####################################################
#   df_chr<-data.frame(df_chr_IDs_dist$dist,design_mat_TAD,design_mat_epigenetic)
#   names(df_chr)<-df_chr_names
#   #predictedValues<-as.numeric(predict(object = model_fit_obj, s = model_tuning_param, newx = as.matrix(df_chr), type = "response"))
#   predictedValues<-as.numeric(predict(object = model_fit_obj, s = model_tuning_param, newx = as.matrix(df_chr), type = "response"))
#   #predictedValues<-as.numeric(predict(object = lm_obj, penalized = as.matrix(df_chr))[,1])
#   predicted_matrix<-matrix(0,length(unique(chr_IDs_mat[,1])),length(unique(chr_IDs_mat[,2])))
#   for(i in 1:nrow(chr_IDs_mat)){
#     predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]<-predictedValues[i]
#     predicted_matrix[chr_IDs_mat[i,2]-start_ID+1,chr_IDs_mat[i,1]-start_ID+1]<-predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]
#   }
#   quantile_log_intensity<-quantile(c(log(intensity+1)),prob = c(0.0,0.99))
#   quantile_predicted_log_intensity<-quantile(predictedValues,prob = c(0.0,0.99))
#   Rstart<-max(quantile_log_intensity[1],quantile_predicted_log_intensity[1])
#   Rend<-min(quantile_log_intensity[2],quantile_predicted_log_intensity[2])
#   Rrange<-c(Rstart,Rend)
#   ## Calculating the original matrix for specified window
#   original_matrix<-log(chr_mat[start_ID:end_ID,start_ID:end_ID]+1)
#   ## Correlation analysis for the specified window
#   corr_vec_generator<-function(df_chr_IDs_dist){
#     N<-max(df_chr_IDs_dist$dist)
#     corr_vec<-rep(NA_real_,N)
#     for(k in 1:N){
#       ids<-which(df_chr_IDs_dist$dist==k)
#       if(length(ids)>10){
#         intensity_val_k<-log(intensity[ids]+1)
#         intensity_val_k_predicted<-predictedValues[ids]
#         corr_vec[k]<-cor(intensity_val_k,intensity_val_k_predicted)
#       }
#       if((k%%100)==0){
#         print(k)
#       }
#     }
#     return(corr_vec)
#   }
#   corr_vec<-corr_vec_generator(df_chr_IDs_dist = df_chr_IDs_dist)
#   non_NA_ids<-which(!is.na(corr_vec))
#   corr_vec_modified<-corr_vec[non_NA_ids]
#   dist_vec<-(1:max(df_chr_IDs_dist$dist))[non_NA_ids]
#   ## Overall correlation
#   corr_overall<-cor(log(intensity+1),as.numeric(predictedValues))
#   return(list(original_matrix,predicted_matrix,Rrange,dist_vec,corr_vec_modified,corr_overall))
# }
# 
# debug(predict_output_gen)
# predict_output<-predict_output_gen(start_ID=101, end_ID= 300, chr_mat = chr_mat_full,TAD_locations = TAD_locations,model_fit_obj = model_fit,model_tuning_param = tuning_param,df_chr_names = df_chr_names[-1])
# 
# par(mfrow=c(1,2))
# plotLoop(x = predict_output[[1]],myrange = predict_output[[3]])
# plotLoop(x = predict_output[[2]],myrange = predict_output[[3]])
# 
# par(mfrow=c(1,1))
# plot(x = predict_output[[4]],y = predict_output[[5]])
# predict_output[[6]]
# 
# #########################################################################################################
# ## Regression of log(gamma_TAD+1) vs. epigenetic
# 
# load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Apps/HiCEpigen/model2/lm_fit.RData")
# 
# load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Apps/HiCEpigen/TAD_locations.RData")
# 
# load(file="/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Apps/HiCEpigen/df_epigen_modified.RData")
# 
# ## Extracting all coefficients
# coeff_all<-as.numeric(coef(lm_fit,which="all"))
# 
# ## Extracting the TAD coefficients
# gamma_TAD_estimated<-coeff_all[(3:(length(TAD_locations)+2))]
# 
# df_epigen




