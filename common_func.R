## Defining common functions
## Loading required libraries
library(glmnet)
library(dplyr)
library(penalized)
library(data.table)

#########################################################################################################
# Function for removing zero IDs
remove_IDs_func <- function(chr_mat_window, signal_threshold = 100){
  non_zeros_IDs_vec <- rep(NA_real_, nrow(chr_mat_window))
  signal_strength_vec <- rep(NA_real_, nrow(chr_mat_window))
  for (i in 1:nrow(chr_mat_window)){
    non_zeros_IDs_vec[i] <- sum(chr_mat_window[i, ]!=0)
    signal_strength_vec[i] <- sum(chr_mat_window[i, ])
  }
  remove_IDs <- which(signal_strength_vec < signal_threshold)
  ## Training the model after removing zeros.
  chr_mat <- chr_mat_window[-remove_IDs, -remove_IDs]
  return(list(chr_mat, remove_IDs))
}

#########################################################################################################
# Defining a function to get the paired row IDs and column corrsponding the start and end IDs and bandwidth parameter (taking into account the symmetric matrix)
pairs_list_generator<-function(start_ID,end_ID,remove_IDs,bandwidth_parameter){
  pairs_list<-list()
  for(i in start_ID:end_ID){
    if(!(i%in%remove_IDs)){
      col_ID_seq<-(i-(bandwidth_parameter-1)):(i+(bandwidth_parameter-1))
      col_ID_seq_positive<-col_ID_seq[which((col_ID_seq>=max(i,start_ID))&(col_ID_seq<=end_ID))]
      pairs_list[[i]]<-data.frame("row_ID"=rep(i,length(col_ID_seq_positive)),"col_ID"=col_ID_seq_positive)
    }
  }
  return(pairs_list)
}

#########################################################################################################
# Defining a function to extract the intensities corresponding to given chr_mat and the pair ids inside chr_IDs_mat
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

#########################################################################################################
## Defining a function to generate the chr object
chr_obj_generator<-function(chr_mat,window_start_ID,window_end_ID,bandwidth_parameter,sampling_fraction=1){
  chr_mat_full<-chr_mat
  window<-c(window_start_ID,window_end_ID)
  chr_mat_window<-chr_mat_full[window[1]:window[2],window[1]:window[2]]
  remove(chr_mat)
  
  #########################################################################################################
  remove_IDs_output<-remove_IDs_func(chr_mat_window = chr_mat_window)
  
  chr_mat<-remove_IDs_output[[1]]
  remove_IDs<-remove_IDs_output[[2]]
  remove(remove_IDs_output)
  
  #########################################################################################################
  ## Number of bins in the chromosome matrix
  N<-dim(chr_mat_window)[1]-length(remove_IDs)
  ## Getting the chr pairs list
  chr_pairs_list<-pairs_list_generator(start_ID = window[1], end_ID = window[2], remove_IDs = remove_IDs, bandwidth_parameter = bandwidth_parameter)
  
  ## Initializing the dataframe for the chromosome IDs and distance
  df_chr_IDs_dist<-do.call(what = rbind,args = chr_pairs_list)
  df_chr_IDs_dist$dist<-abs(df_chr_IDs_dist$row_ID-df_chr_IDs_dist$col_ID)
  print(nrow(df_chr_IDs_dist))
  
  #########################################################################################################
  ## Randomly sample from 5% of data
  set.seed(1)
  randomly_sampled_row_IDs<-sample(x = 1:nrow(df_chr_IDs_dist),size = sampling_fraction*nrow(df_chr_IDs_dist),replace = F) 
  #head(randomly_sampled_row_IDs)
  df_chr_IDs_dist_subset<-df_chr_IDs_dist[sort(randomly_sampled_row_IDs),]
  
  #####################################################
  ## Getting the IDs matrix
  chr_IDs_mat<-as.matrix(df_chr_IDs_dist_subset[,c(1,2)])
  attr(chr_IDs_mat,"dimnames")<-NULL
  
  #chr_IDs_mat_full<-as.matrix(df_chr_IDs_dist[,c(1,2)])
  #attr(chr_IDs_mat_full,"dimnames")<-NULL
  
  #str(chr_IDs_mat)
  #head(chr_IDs_mat)
  #str(chr_IDs_mat_full)
  #save(chr_IDs_mat,)
  
  #####################################################
  intensity<-chr_intensity_extractR(chr_mat = chr_mat_full,chr_IDs_mat = chr_IDs_mat)
  #intensity_full<-chr_intensity_extractR(chr_mat = chr_mat,chr_IDs_mat = chr_IDs_mat_full)
  
  return(list(chr_mat_full,chr_mat_window,chr_mat,remove_IDs,df_chr_IDs_dist_subset,chr_IDs_mat,intensity))
}

#########################################################################################################
# Defining a function to get the design matrix for TAD locations and the given chr_IDs_mat obtained after removing zero IDs
design_mat_TAD_generator<-function(interaction_mat,TAD_locations){
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

#########################################################################################################
# Defining a function to get the design matrix for epigenetic states
design_mat_epigenetic_generator<-function(df_epigen_modified,chr_IDs_mat){
  design_mat_epigenetic<-matrix(NA_real_,nrow(chr_IDs_mat),2*ncol(df_epigen_modified))
  for (k in 1:nrow(chr_IDs_mat)){
    design_mat_epigenetic[k,1:35]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,1],])
    design_mat_epigenetic[k,36:70]<-as.numeric(df_epigen_modified[chr_IDs_mat[k,2],])
    if((k%%1000)==0){
      print(k)
    }
  }
  return(design_mat_epigenetic)
}

# #########################################################################################################
# ## Preprocessing and getting TAD locations for Gm12878
# df_TAD_Gm12878<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/TAD/dp_raw_pen0.1_newestchr1.txt", header = FALSE, sep = "\t")
# str(df_TAD_Gm12878)
# TADs_Gm12878<-unique(sort(as.numeric(c(df_TAD_Gm12878[-1,1],df_TAD_Gm12878[-1,2]))))
# TAD_locations_Gm12878<-TADs_Gm12878[-c(1,length(TADs_Gm12878))]
# str(TAD_locations_Gm12878)
# save(TAD_locations_Gm12878,file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Apps/HiCEpigen/TAD_locations_Gm12878.RData")
# 
# ## Preprocessing and getting TAD locations for K562
# df_TAD_K562<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/K562/10kb/tadk562.txt", header = FALSE, sep = " ")
# str(df_TAD_K562)
# head(df_TAD_K562)
# TADs_K562<-unique(sort(as.numeric(c(df_TAD_K562[which(df_TAD_K562$V1=="chr1"),2],df_TAD_K562[which(df_TAD_K562$V1=="chr1"),3]))))
# TAD_locations_K562<-TADs_K562[-c(1,length(TADs_K562))]
# str(TAD_locations_K562)
# save(TAD_locations_K562,file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Apps/HiCEpigen/TAD_locations_K562.RData")
# 
# intersect(TAD_locations_Gm12878,TAD_locations_K562)

# #########################################################################################################
# ## Preprocessing and getting TAD locations for Gm12878
# df_TAD_Gm12878<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/TAD_boundaries/Gm12878_TAD/dp_raw_pen0.1_newest.chr1.txt", header = FALSE, sep = "\t")
# str(df_TAD_Gm12878)
# TADs_Gm12878<-unique(sort(as.numeric(c(df_TAD_Gm12878[-1,1],df_TAD_Gm12878[-1,2]))))
# TAD_locations_Gm12878<-TADs_Gm12878[-c(1,length(TADs_Gm12878))]
# str(TAD_locations_Gm12878)
# 
# df_TAD_K562<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/TAD_boundaries/K562_TAD/dp_raw_pen0.1_newest.chr1.txt", header = FALSE, sep = "\t")
# str(df_TAD_K562)
# TADs_K562<-unique(sort(as.numeric(c(df_TAD_K562[-1,1],df_TAD_K562[-1,2]))))
# TAD_locations_K562<-TADs_K562[-c(1,length(TADs_K562))]
# str(TAD_locations_K562)
# 
# intersect(TAD_locations_Gm12878,TAD_locations_K562)

#########################################################################################################
## Preprocessing and getting TAD locations for Gm12878 (joint call)
df_TAD_Gm12878<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/TAD_boundaries/Joint_call/OnTADraw_pen0.1_max200_meannorm_chr1_data1.tad.txt", header = FALSE, sep = "\t")
str(df_TAD_Gm12878)
TADs_Gm12878<-unique(sort(as.numeric(c(df_TAD_Gm12878[-1,1],df_TAD_Gm12878[-1,2]))))
TAD_locations_Gm12878<-TADs_Gm12878[-c(1,length(TADs_Gm12878))]
str(TAD_locations_Gm12878)

df_TAD_K562<-read.table(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/TAD_boundaries/Joint_call/OnTADraw_pen0.1_max200_meannorm_chr1_data2.tad.txt", header = FALSE, sep = "\t")
str(df_TAD_K562)
TADs_K562<-unique(sort(as.numeric(c(df_TAD_K562[-1,1],df_TAD_K562[-1,2]))))
TAD_locations_K562<-TADs_K562[-c(1,length(TADs_K562))]
str(TAD_locations_K562)

intersect(TAD_locations_Gm12878,TAD_locations_K562)

#########################################################################################################
## Epigenetic data preprocessing
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/df_epigen_output_Gm12878.RData")
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/states_Gm12878.RData")
epigen_mat_Gm12878<-df_epigen_output
states_Gm12878<-states

load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/df_epigen_output_K562.RData")
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/states_K562.RData")
epigen_mat_K562_raw<-df_epigen_output
states_K562<-states

states_K562_permute_IDs<-rep(NA_real_,36)
for(i in 1:length(states_Gm12878)){
  states_K562_permute_IDs[i]<-which(states_K562==states_Gm12878[i])
}
prod(states_Gm12878==states_K562[states_K562_permute_IDs])

epigen_mat_K562<-epigen_mat_K562_raw[,states_K562_permute_IDs]

## Defining a function for epigenetic normalization
epiNormFunc <- function(epigen_mat,states){
  epigen_mat_normalized<-t(apply(X = epigen_mat,MARGIN = 1,FUN = function(x){
    y<-x/(sum(x)+1e-4)
    return(y)
  }))
  df_epigen_modified<-data.frame(epigen_mat_normalized[,-which(states=="Zero")])
  names(df_epigen_modified)<-states[-which(states=="Zero")]
  return(df_epigen_modified)
}

## Getting final epigenetic dataframe for Gm12878
df_epigen_modified_Gm12878<-epiNormFunc(epigen_mat = epigen_mat_Gm12878, states = states_Gm12878)
str(df_epigen_modified_Gm12878)

## Getting final epigenetic dataframe for K562
df_epigen_modified_K562<-epiNormFunc(epigen_mat = epigen_mat_K562, states = states_Gm12878)
str(df_epigen_modified_K562)


#########################################################################################################
# ## Epigenetic data preprocessing
# df_epigen_Gm12878_v1<-fread(input = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/encodestate_chr1_10kb.txt",data.table = F)
# 
# str(df_epigen_Gm12878_v1)
# 
# df_epigen_Gm12878_v1[,1]
# 
# epigen_mat_Gm12878_v1<-as.matrix(df_epigen_Gm12878_v1)
# 
# epigen_mat_normalized_Gm12878_v1<-t(apply(X = epigen_mat_Gm12878_v1,MARGIN = 1,FUN = function(x){
#   y<-x/(sum(x)+1e-4)
#   return(y)
# }))
# 
# df_epigen_modified_Gm12878_v1<-data.frame(epigen_mat_normalized_Gm12878_v1[,-1])
# str(df_epigen_modified_Gm12878_v1)
# names(df_epigen_modified_Gm12878_v1)
# 
# 
# ## Epigenetic data preprocessing
# load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/df_epigen_output_Gm12878.RData")
# 
# str(df_epigen_output)
# 
# df_epigen_Gm12878_v2<-df_epigen_output
# 
# permute_IDs<-matrix(c(1:36,rep(NA,36)),36,2)
# for(i in 1:36){
#   for(j in 1:36){
#     if(prod(df_epigen_Gm12878_v2[,j]==df_epigen_Gm12878_v1[,i])==1){
#       permute_IDs[i,2]<-j
#     }
#   }
# }
# 
# 
# prod(df_epigen_Gm12878_v1==df_epigen_Gm12878_v2[,permute_IDs[,2]])


