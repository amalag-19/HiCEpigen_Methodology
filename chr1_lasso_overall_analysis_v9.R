## Script for preprocessing the epigenetic data from bigbed format

library(data.table)

## Reading the epigenetic raw data for the Gm12878 cell
df_epigen<-fread(input = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/GM12878_ENCODE36.bigBed",data.table = F)

#df_epigen<-fread(input = "/Users/Amal/Box Sync/PSU/Fall 2018/BioStat_Research/Data/Epigenetic/K562_ENCODE36.bigBed",data.table = F)

## Checking the df_epigen dataframe
str(df_epigen)
head(df_epigen)

## Subsetting the whole dataframe for chromosome 1
df_epigen_chr1<-df_epigen[which(df_epigen$chrom=="chr1"),]
str(df_epigen_chr1)

## Defining the bin sequence. Here I add 2 to match the HiC data chromosome 1 matrix dimensions.
bin_seq<-seq(1,(max(df_epigen_chr1$chromEnd)/10000)+2)
length(bin_seq)

## Extracting all unique 36 states
states<-unique(df_epigen_chr1$name)

## Initializing the window start and end IDs
window_start_IDs<-rep(NA_integer_,length(bin_seq))
window_end_IDs<-rep(NA_integer_,length(bin_seq))

## Function (used later) to give state count vector (length of 36 corresponding to state vector) for the current window.
state_counter<-function(df,window_start_ID,window_end_ID,states,max_val,min_val){
  total_count_vec<-floor((min_val-max_val+5)/200)
  state_count_vec<-rep(0,length(states))
  for(j in 1:length(states)){
    state_j_row_IDs<-which(df$name==states[j])
    if(length(state_j_row_IDs)>=1){
      state_count_vec[j]<-state_count_vec[j]+sum(total_count_vec[state_j_row_IDs])
    }
  }
  return(state_count_vec)
}

## Initializing the output of epigenetic data
df_epigen_output<-matrix(NA_integer_,length(bin_seq),length(states))
ptm<-proc.time()
for(i in 1:length(bin_seq)){
  ## Extracting the window start and end IDs assuming 10 kb resolution.
  window_start_IDs[i]<-(i-1)*10000+1
  window_end_IDs[i]<-i*10000
  ## Getting the maximum values of chromStart and current window start ID
  max_val<-apply(X = matrix(df_epigen_chr1$chromStart),MARGIN = 1,FUN = function(x){
    max(x,window_start_IDs[i])
  })
  ## Getting the minimum values of chromEnd and current window end ID
  min_val<-apply(X = matrix(df_epigen_chr1$chromEnd),MARGIN = 1,FUN = function(x){
    min(x,window_end_IDs[i])
  })
  ## Extracting the row IDs when the max_val is less than equal to min_val. These are the rows we will care about for current window. Note that this covers cases for intervals spanning more than 1 10kb window.
  subset_row_IDs<-which(max_val<=min_val)
  ## Subsetting the dataframe and counting states for each interval in units of 200 b.p.
  if(length(subset_row_IDs)>=1){
    df_epigen_chr1_subset<-df_epigen_chr1[subset_row_IDs,]
    max_val_subset<-max_val[subset_row_IDs]
    min_val_subset<-min_val[subset_row_IDs]
    state_count_vec<-state_counter(df = df_epigen_chr1_subset, window_start_ID = window_start_IDs[i], window_end_ID = window_end_IDs[i], states = states, max_val = max_val_subset, min_val = min_val_subset)
  }else{
    state_count_vec<-rep(0,length(states))
  }
  df_epigen_output[i,]<-state_count_vec
  ## Printing out time and bin ID every 100 bins to keep track of program (total estimated time is ~7 hrs.)
  if((i%%100)==0){
    print(proc.time()-ptm)
    print(i)
    ptm<-proc.time()
  }
}
