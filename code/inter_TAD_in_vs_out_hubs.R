library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------------
##Utils. Fn

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#--------------------------------
TAD_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_GM12878_TAD_pval_tbl.Rda"
top_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda"
res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
hub_tbl<-get_tbl_in_fn(top_hub_file)
TAD_tbl<-get_tbl_in_fn(TAD_file)

hub_chr_set<-unique(hub_tbl$chr)
map(hub_chr_set,function(chromo){
  chr_spec_res<-get_tbl_in_fn(paste0(res_file,chromo"_spec_res.Rda"))
  
})