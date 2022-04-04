library(tidyverse)
library(GenomicRanges)
library(furrr)
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
inter_file<-"./data/intersection_out_tbl/HMEC_TAD_cl_intersect.Rda"


TAD_cl_inter_tbl<-get_tbl_in_fn(inter_file)

TAD_cl_inter_tbl %>% 
  ggplot(.,aes(max.inter))+
  geom_histogram()


TAD_cl_inter_tbl<-TAD_cl_inter_tbl %>% 
  mutate(max.GRange=pmap(list(inter.cl.GRange,inter.stat),function(inter.cl.GRange,inter.stat){
    inter.cl.GRange[[which.max(inter.stat)]]
  })) 

max_cl_GRange_l<-GRangesList(TAD_cl_inter_tbl$max.GRange)
TAD_Grange_l<-GRangesList(TAD_cl_inter_tbl$GRange)
tad_cover<-unlist(lapply(width(IRanges::intersect(TAD_Grange_l,max_cl_GRange_l)),sum))/(as.numeric(width(TAD_Grange_l)))
cl_cover<-unlist(lapply(width(IRanges::intersect(TAD_Grange_l,max_cl_GRange_l)),sum))/as.numeric(unlist(lapply(width(max_cl_GRange_l),sum)))
tibble(tad=tad_cover,cl=cl_cover) %>% 
  ggplot(.,aes(tad,cl))+
  geom_density_2d_filled()
