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
chr_pairs_list<-pairs_list_generator(start_ID = 1,end_ID = N,bandwidth_parameter = 200)

## Initializing the dataframe for the chromosome IDs and distance
df_chr_IDs_dist<-do.call(what = rbind,args = chr_pairs_list)
df_chr_IDs_dist$dist<-abs(df_chr_IDs_dist$row_ID-df_chr_IDs_dist$col_ID)
#str(df_chr_IDs_dist)

#####################################################
## Getting the IDs matrix
chr_IDs_mat<-as.matrix(df_chr_IDs_dist[,c(1,2)])
attr(chr_IDs_mat,"dimnames")<-NULL

str(chr_IDs_mat)

#########################################################################################################
## Defining a function to generate fragment matrix
fitHiC_fragment_mat_generator<-function(chr_mat,chromosome_number,cell){
  fragment_mat<-matrix(0,nrow = N,ncol = 5)
  fragment_mat[,1]<-rep(1,N)
  fragment_mat[,3]<-(1:N)
  fragment_mat[,4]<-(rowSums(chr_mat))
  write.table(fragment_mat, file = paste0("/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_processed_fitHiC/chr",chromosome_number,"_",cell,"_fragment_mat.txt"), sep = "\t",row.names = FALSE,col.names = FALSE)
  return(fragment_mat)
}

## Defining a function to generate interaction matrix
fitHiC_interaction_mat_generator<-function(chr_mat,chr_IDs_mat,chromosome_number,cell){
  n_samples<-nrow(chr_IDs_mat)
  interaction_mat<-matrix(0,nrow = n_samples,ncol = 5)
  interaction_mat[,1]<-rep(1,n_samples)
  interaction_mat[,2]<-chr_IDs_mat[,1]
  interaction_mat[,3]<-rep(1,n_samples)
  interaction_mat[,4]<-chr_IDs_mat[,2]
  for (k in 1:n_samples){
    interaction_mat[k,5]<-chr_mat[chr_IDs_mat[k,1],chr_IDs_mat[k,2]]
    if((k%%100000)==0){
      print(k)
    }
  }
  interaction_mat_modified<-interaction_mat[which(interaction_mat[,5]>=100),]
  write.table(interaction_mat_modified, file = paste0("/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_processed_fitHiC/chr",chromosome_number,"_",cell,"_interaction_mat.txt"), sep = "\t",row.names = FALSE,col.names = FALSE)
  return(interaction_mat_modified)
}

fitHiC_fragment_mat<-fitHiC_fragment_mat_generator(chr_mat = chr_mat, chromosome_number = 1,cell = "Gm12878")

fitHiC_interaction_mat<-fitHiC_interaction_mat_generator(chr_mat = chr_mat, chr_IDs_mat=chr_IDs_mat, chromosome_number = 1,cell = "Gm12878")

interaction_mat<-fitHiC_output[[2]]
chromosome_number<-1
cell<-"Gm12878"
N<-dim(chr_mat)[1]
chunks<-c(seq(1,N^2,by = 1000000),((N^2)+1))
for(i in 1:(length(chunks)-1)){
  write.table(interaction_mat[chunks[i]:(chunks[i+1]-1),], file = paste0("/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_processed_fitHiC/chr",chromosome_number,"_",cell,"_interaction_mat.txt"), sep = "\t",row.names = FALSE,col.names = FALSE,append = T)
  if((i%%100)==0){
    print(i)
  }
}


#########################################################################################################
## Loading the p value matrix from fitHiC output
load(file = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/HiC/Gm12878/Primary/10kb/chr_processed_fitHiC/chr1_Gm12878_fitHiC_output/")

ids<-which(chr1_Gm12878_fitHiC_output[,6]<0.05)

sample_length<-length(which(chr1_Gm12878_fitHiC_output[,6]<0.05))

interaction_mat<-chr1_Gm12878_fitHiC_output[ids,]



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
    if(interaction_mat[i,4]>=interaction_mat[i,2]){
      for (j in 1:length(TAD_locations)){
        if((TAD_locations[j]>=interaction_mat[i,2])&(TAD_locations[j]<=interaction_mat[i,4])){
          design_mat[i,j]<-1
        }
      }
    }
  }
  return(design_mat)
}

df_chr<-data.frame("Intensity"=interaction_mat[,5],"Distance"=abs(interaction_mat[,4]-interaction_mat[,2]),design_mat)

penFactor = c(0, rep(1,dim(df_chr)[2]-2))

lasso=glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=1,nlambda=100,upper.limits = rep(0,ncol(df_chr)),penalty.factor = penFactor)

cv.lasso = cv.glmnet(x=as.matrix(df_chr[,-1]),y=log(df_chr[,1]+1),alpha=alpha_input,nfolds=10,upper.limits = rep(0,ncol(df_chr)), penalty.factor = penFactor)

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



#########################################################################################################
