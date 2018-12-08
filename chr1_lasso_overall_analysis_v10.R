#########################################################################################################
## Building linear model with constraints on 1-500 window (clust version is 1% subsampling of whole chromosome) with TAD locations, and L2 penalty over epigenetic variables as predictors (including the step of removing zeros before training the model and adding them back after.)
#########################################################################################################

#########################################################################################################
## Getting final dataframe for GM12878
#########################################################################################################
source(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Code/penalized_full_models/common_func.R")

## Defining the window size
window<-c(1,500)

#########################################################################################################
## Getting the chr object for Gm12878
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_mat_preprocessed/chr1_mat.RData")
chr_obj_Gm12878<-chr_obj_generator(chr_mat = chr_mat,window_start_ID = window[1], window_end_ID = window[2], bandwidth_parameter = 500, sampling_fraction = 1)
remove(chr_mat)

## Getting the chr object for K562
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/K562/10kb/chr_mat_preprocessed/chr1_mat.RData")
chr_obj_K562<-chr_obj_generator(chr_mat = chr_mat, window_start_ID = window[1], window_end_ID = window[2], bandwidth_parameter = 500, sampling_fraction = 1)
remove(chr_mat)

#########################################################################################################
## Getting TADs for Gm12878 and K562 within the specified window
TAD_locations_Gm12878_window<-TAD_locations_Gm12878[which(((TAD_locations_Gm12878 >= window[1]) & (TAD_locations_Gm12878 <= window[2])) & (!TAD_locations_Gm12878 %in% chr_obj_Gm12878[[4]]))]

TAD_locations_K562_window<-TAD_locations_K562[which(((TAD_locations_K562 >= window[1]) & (TAD_locations_K562 <= window[2])) & (!TAD_locations_K562 %in% chr_obj_K562[[4]]))]

## Taking union of TADs for Gm12878 and K562 within the specified window
TAD_locations_union_window<-unique(sort(union(TAD_locations_Gm12878_window,TAD_locations_K562_window)))

#####################################################
## Getting the design matrix for TAD boundaries for Gm12878 and K562
design_mat_TAD_Gm12878 <- design_mat_TAD_generator(interaction_mat = chr_obj_Gm12878[[6]],TAD_locations = TAD_locations_union_window)
head(design_mat_TAD_Gm12878)

design_mat_TAD_K562 <- design_mat_TAD_generator(interaction_mat = chr_obj_K562[[6]],TAD_locations = TAD_locations_union_window)
head(design_mat_TAD_K562)

## Getting the combined design_mat_TAD
design_mat_TAD <- rbind(cbind(design_mat_TAD_Gm12878, matrix(0, nrow = nrow(design_mat_TAD_Gm12878), ncol = length(TAD_locations_union_window))), cbind(matrix(0, nrow = nrow(design_mat_TAD_K562), ncol = length(TAD_locations_union_window)), design_mat_TAD_K562)) 

str(design_mat_TAD)

#########################################################################################################
## Getting the design matrix for epigenetic for Gm12878 and K562
design_mat_epigenetic_Gm12878 <- design_mat_epigenetic_generator(df_epigen_modified = df_epigen_modified_Gm12878, chr_IDs_mat = chr_obj_Gm12878[[6]])
head(design_mat_epigenetic_Gm12878)

design_mat_epigenetic_K562 <- design_mat_epigenetic_generator(df_epigen_modified = df_epigen_modified_K562, chr_IDs_mat = chr_obj_K562[[6]])
head(design_mat_epigenetic_K562)

## Getting the combined design_mat_epigenetic
design_mat_epigenetic <- rbind(design_mat_epigenetic_Gm12878, design_mat_epigenetic_K562)

str(design_mat_epigenetic)

#####################################################
## Building the combined dataframe
df_chr<-data.frame(c(chr_obj_Gm12878[[7]],chr_obj_K562[[7]]), c(chr_obj_Gm12878[[5]]$dist, chr_obj_K562[[5]]$dist), design_mat_TAD, design_mat_epigenetic)

## Processing column names of the dataframe
TAD_numbers<-as.character(TAD_locations_union_window) 
TAD_names<-paste("TAD_",TAD_numbers, sep="")
Epi_names<-states_Gm12878[-which(states_Gm12878=="Zero")]
df_chr_names<-c("Intensity","Distance",paste0(TAD_names,"_Gm12878"),paste0(TAD_names,"_K562"),paste0(Epi_names,"_index1"),paste0(Epi_names,"_index2"))
names(df_chr)<-df_chr_names

##########################################################################################################
## Defining the penalty factor to impose no penalty on distance and TAD locations but L2 penalty on epigenetic variables
penFactor<-c(rep(0,ncol(design_mat_TAD)+1),rep(1,ncol(design_mat_epigenetic)))
#penFactor<-c(rep(1,ncol(design_mat_TAD)+1),rep(0,ncol(design_mat_epigenetic)))

model_fit<-glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=0,nlambda=100,upper.limits = c(rep(0,ncol(design_mat_TAD)+1),rep(Inf,ncol(design_mat_epigenetic))),penalty.factor = penFactor)

cvl_fit<-cv.glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=0,nfolds=10,upper.limits = c(rep(0,ncol(design_mat_TAD)+1),rep(Inf,ncol(design_mat_epigenetic))),penalty.factor = penFactor)

## Extracting all coefficients
tuning_param<-cvl_fit$lambda.min
print(tuning_param)

## beta estimates for best lambda
model_coeff_obj<-coef(model_fit,s=tuning_param)
model_coeff_obj

model_coeff_vec<-rep(0,length(model_coeff_obj@Dimnames[[1]]))
model_coeff_vec[model_coeff_obj@i + 1]<-model_coeff_obj@x

df_coeff<-data.frame(name = model_coeff_obj@Dimnames[[1]], coefficient = model_coeff_vec,stringsAsFactors = F)

##########################################################################################################
source('~/Box Sync/PSU/Fall 2018/BioStat_Research/Code/Yu_code.R')
predict_output_gen<-function(window_start_ID, window_end_ID, bandwidth_parameter, chr_obj_Gm12878, chr_obj_K562, TAD_locations_union_window, model_fit_obj, model_tuning_param = 0, df_chr_names){
  window_size<-window_end_ID-window_start_ID+1
  plotRange<-c(start_ID,end_ID)
  # Getting chr objects for the two cells for local windows
  chr_obj_Gm12878_window <- chr_obj_generator(chr_mat = chr_obj_Gm12878[[1]], window_start_ID = window_start_ID, window_end_ID = window_end_ID, bandwidth_parameter = bandwidth_parameter, sampling_fraction = 1)
  chr_obj_K562_window <- chr_obj_generator(chr_mat = chr_obj_K562[[1]], window_start_ID = window_start_ID, window_end_ID = window_end_ID, bandwidth_parameter = bandwidth_parameter, sampling_fraction = 1)
  
  #####################################################
  ## Getting the local design matrices for TAD boundaries for Gm12878 and K562
  design_mat_TAD_Gm12878_window <- design_mat_TAD_generator(interaction_mat = chr_obj_Gm12878_window[[6]], TAD_locations = TAD_locations_union_window)
  
  design_mat_TAD_K562_window <- design_mat_TAD_generator(interaction_mat = chr_obj_K562_window[[6]], TAD_locations = TAD_locations_union_window)
  
  ## Getting the combined design_mat_TAD
  design_mat_TAD_window <- rbind(cbind(design_mat_TAD_Gm12878_window, matrix(0, nrow = nrow(design_mat_TAD_Gm12878_window), ncol = length(TAD_locations_union_window))), cbind(matrix(0, nrow = nrow(design_mat_TAD_K562_window), ncol = length(TAD_locations_union_window)), design_mat_TAD_K562_window)) 
  
  str(design_mat_TAD_window)
  
  #####################################################
  ## Getting the design matrix for epigenetic for Gm12878 and K562
  design_mat_epigenetic_Gm12878_window <- design_mat_epigenetic_generator(df_epigen_modified = df_epigen_modified_Gm12878, chr_IDs_mat = chr_obj_Gm12878_window[[6]])
  
  design_mat_epigenetic_K562_window <- design_mat_epigenetic_generator(df_epigen_modified = df_epigen_modified_K562, chr_IDs_mat = chr_obj_K562_window[[6]])
  
  ## Getting the combined design_mat_epigenetic
  design_mat_epigenetic_window <- rbind(design_mat_epigenetic_Gm12878_window, design_mat_epigenetic_K562_window)
  
  str(design_mat_epigenetic_window)
  
  #####################################################
  ## Building the combined dataframe
  df_chr<-data.frame(c(chr_obj_Gm12878_window[[7]], chr_obj_K562_window[[7]]), c(chr_obj_Gm12878_window[[5]]$dist, chr_obj_K562_window[[5]]$dist), design_mat_TAD, design_mat_epigenetic)
  names(df_chr)<-df_chr_names
  
  ## Prediction
  predictedValues<-as.numeric(predict(object = model_fit_obj, s = model_tuning_param, newx = as.matrix(df_chr[, -1]), type = "response"))
  
  ## Original values
  originalValues<-log(df_chr[, 1]+1)
  
  quantile_normalisation <- function(df){
    df_rank <- apply(df,2,rank,ties.method="min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    index_to_mean <- function(my_index, my_mean){
      return(my_mean[my_index])
    }
    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
    rownames(df_final) <- rownames(df)
    return(df_final)
  }
  
  ## Putting original and predicted in one dataframe
  #df_org_pred<-data.frame("original"=originalValues,"predicted"=predictedValues)
  # ## Quantile normalization
  # df_org_pred_normalized<-quantile_normalisation(df_org_pred)
  # predictedValues_quantile_normalized<-df_org_pred_normalized[,2]
  # originalValues_quantile_normalized<-df_org_pred_normalized[,1]
  
  ## Function to do Quantile normalization with sub-diagonal
  quantileNorm_output_generator<-function(df_chr_IDs_dist, originalValues, predictedValues){
    N<-max(df_chr_IDs_dist$dist)
    predictedValues_quantile_normalized_subdiagonal<-rep(NA_real_,length(originalValues))
    originalValues_quantile_normalized_subdiagonal<-rep(NA_real_,length(predictedValues))
    for(k in 0:N){
      ids<-which(df_chr_IDs_dist$dist==k)
      if(length(ids)>=2){
        df_org_pred_k<-data.frame("original"=originalValues[ids],"predicted"=predictedValues[ids])
        ## Quantile normalization
        df_org_pred_normalized_k<-quantile_normalisation(df_org_pred_k)
        new_predictedValues_k<-df_org_pred_normalized_k[,2]
        new_originalValues_k<-df_org_pred_normalized_k[,1]
      }else{
        new_predictedValues_k<-predictedValues[ids]
        new_originalValues_k<-originalValues[ids]
      }
      predictedValues_quantile_normalized_subdiagonal[ids]<-new_predictedValues_k
      originalValues_quantile_normalized_subdiagonal[ids]<-new_originalValues_k
      #print(k)
    }
    return(list(predictedValues_quantile_normalized_subdiagonal,originalValues_quantile_normalized_subdiagonal))
  }
  
  quantileNorm_output <- quantileNorm_output_generator(df_chr_IDs_dist = chr_obj_Gm12878_window[[4]], originalValues = df_chr[,1], predictedValues = predictedValues)
  
  
  ## Function to get original and predicted matrices, common range, and correlation output
  result_generator<-function(predictedVal, originalVal, df_chr_IDs_dist, chr_IDs_mat){
    #predictedValues<-as.numeric(predict(object = lm_obj, penalized = as.matrix(df_chr))[,1])
    predicted_matrix<-matrix(0,window_size,window_size)
    original_matrix<-matrix(0,window_size,window_size)
    for(i in 1:nrow(chr_IDs_mat)){
      predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]<-predictedVal[i]
      predicted_matrix[chr_IDs_mat[i,2]-start_ID+1,chr_IDs_mat[i,1]-start_ID+1]<-predicted_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]
      original_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]<-originalVal[i]
      original_matrix[chr_IDs_mat[i,2]-start_ID+1,chr_IDs_mat[i,1]-start_ID+1]<-original_matrix[chr_IDs_mat[i,1]-start_ID+1,chr_IDs_mat[i,2]-start_ID+1]
    }
    quantile_log_intensity<-quantile(originalVal,prob = c(0.0,0.99))
    quantile_predicted_log_intensity<-quantile(predictedVal,prob = c(0.0,0.99))
    Rstart<-max(quantile_log_intensity[1],quantile_predicted_log_intensity[1])
    Rend<-min(quantile_log_intensity[2],quantile_predicted_log_intensity[2])
    Rrange<-c(Rstart,Rend)
    ## Correlation analysis for the specified window
    ## Correlation analysis for the specified window
    corr_vec_generator<-function(df_chr_IDs_dist,org_log_intensity,pred_log_intensity){
      N<-max(df_chr_IDs_dist$dist)
      corr_vec<-rep(NA_real_,N)
      for(k in 0:N){
        ids<-which(df_chr_IDs_dist$dist==k)
        if(length(ids)>10){
          intensity_val_k<-org_log_intensity[ids]
          intensity_val_k_predicted<-pred_log_intensity[ids]
          corr_vec[k]<-cor(intensity_val_k,intensity_val_k_predicted)
        }
        if((k%%100)==0){
          print(k)
        }
      }
      return(corr_vec)
    }
    corr_vec<-corr_vec_generator(df_chr_IDs_dist = df_chr_IDs_dist, org_log_intensity = originalVal, pred_log_intensity = predictedVal)
    non_NA_ids<-which(!is.na(corr_vec))
    corr_vec_modified<-corr_vec[non_NA_ids]
    dist_vec<-(1:max(df_chr_IDs_dist$dist))[non_NA_ids]
    ## Overall correlation
    corr_overall<-cor(originalVal,as.numeric(predictedVal))
    return(list(original_matrix,predicted_matrix,Rrange,dist_vec,corr_vec_modified,corr_overall))
  }
  
  ## Getting results for two cells
  ## For Gm12878
  result_list[[1]]<-result_generator(predictedValues = quantileNorm_output[[1]][1:length(chr_obj_Gm12878_window[[7]])], originalValues = quantileNorm_output[[1]][1:length(chr_obj_Gm12878_window[[7]])])
  ## For K562
  result_list[[2]]<-result_generator(predictedValues = quantileNorm_output[[1]][(length(chr_obj_Gm12878_window[[7]])+1):(length(chr_obj_Gm12878_window[[7]])+length(chr_obj_K562_window[[7]]))], originalValues = quantileNorm_output[[1]][(length(chr_obj_Gm12878_window[[7]])+1):(length(chr_obj_Gm12878_window[[7]])+length(chr_obj_K562_window[[7]]))])
  
  return(result_list)
}

#debug(predict_output_gen)
predict_output<-predict_output_gen(start_ID=201, end_ID= 400, remove_IDs = remove_IDs, chr_mat = chr_mat_full,TAD_locations = TAD_locations,model_fit_obj = model_fit,model_tuning_param = tuning_param,df_chr_names = df_chr_names[-1])

par(mfrow=c(1,2))
plotLoop(x = predict_output[[1]][[1]],myrange = predict_output[[1]][[3]])
plotLoop(x = predict_output[[1]][[2]],myrange = predict_output[[1]][[3]])




