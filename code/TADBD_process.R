library(tidyverse)
library(TADBD)
library(vroom)
library(Matrix)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%
           filter(!(is.na(X3)))%>%filter(X1!=X2))
}

full_f_mat<-function(cl_mat,res){
  
  
  range_5kb<-range(unique(c(cl_mat$X1,cl_mat$X2)))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  cl_mat<-cl_mat %>% 
    bind_rows(.,tibble(X1=cl_mat$X2,X2=cl_mat$X1,X3=cl_mat$X3)) %>% 
    bind_rows(.,tibble(X1=bin_5kb,X2=bin_5kb,X3=0))
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$X1)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$X2)]
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=cl_mat$X3,dimnames = list(bin_5kb,bin_5kb))
  return(chr_mat)
}
#-----------------------------------------
HiC_data_folder<-"~/Documents/multires_bhicect/data/HMEC/"
tmp_res<-"10kb" # resolution at which Rao-arrowhead detected their domains in original paper

chr_set<-str_split_fixed(list.files(paste0(HiC_data_folder,tmp_res)),"\\.",2)[,1]

TADBD_res_tbl<-do.call(bind_rows,lapply(chr_set,function(chromo){
  message(chromo)
  chr_dat<-hic_dat_in(HiC_data_folder,tmp_res,chromo)
  tmp_count<-as.matrix(full_f_mat(chr_dat,res_num[tmp_res]))
  rm(chr_dat)
  gc(full = T)
  df_result <- TADBD(tmp_count)
  gc(full = T)
  return(df_result$domains %>% 
    as_tibble %>% 
    mutate(chr=chromo,
           start = as.numeric(rownames(tmp_count)[from.id2+1]),
           end = as.numeric(rownames(tmp_count)[to.id2+1])) %>% 
      dplyr::select(chr,start,end))
  rm(tmp_count)
  
}))
