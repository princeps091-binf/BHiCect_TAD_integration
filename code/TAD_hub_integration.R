library(tidyverse)
library(parallel)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#---------------------------------------------

inter_file<-"./data/intersection_out_tbl/HMEC_TAD_cl_intersect.Rda"
top_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_trans_res_dagger_tbl.Rda"

TAD_cl_inter_tbl<-get_tbl_in_fn(inter_file)
hub_tbl<-get_tbl_in_fn(top_hub_file)

tmp_TAD_Grange<-lapply(TAD_cl_inter_tbl$GRange,function(x)narrow(x,start=1,end=-2))
tmp_cl_Grange<-TAD_cl_inter_tbl$inter.cl.GRange
length(tmp_TAD_Grange)
cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(dplyr)
  library(GenomicRanges)
  print("node ready")
})
clusterExport(cl, c("tmp_TAD_Grange", "tmp_cl_Grange"))
tmp_stat_l <- parLapply(cl, 1:5, function(x) {
  tmp_cl_set<-GRangesList(tmp_cl_Grange[[x]])
  tmp_l<-GRangesList(lapply(1:length(tmp_cl_set),function(y)return(tmp_TAD_Grange[[x]])))
  
  
  tad_cl_inter<-intersect(tmp_cl_set,tmp_l)
  inter_w<-unlist(lapply(width(IRanges::intersect(tmp_cl_set,tmp_l)),sum))
  return(tibble(TAD.w=width(tmp_TAD_Grange[[x]]),cl.w=unlist(lapply(width(tmp_cl_set),sum)),inter.w=inter_w))
  })
stopCluster(cl)
rm(cl)
