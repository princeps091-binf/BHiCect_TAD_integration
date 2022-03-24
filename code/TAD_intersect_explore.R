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
inter_file<-"./data/intersection_out_tbl/GM12878_TAD_cl_intersect.Rda"


TAD_cl_inter_tbl<-get_tbl_in_fn(inter_file)

TAD_cl_inter_tbl %>% 
  ggplot(.,aes(max.inter))+
  geom_histogram()

TAD_cl_inter_tbl %>% 
  mutate(max.inter.res=str_split_fixed(max.cl,"_",2)[,1]) %>% 
  ggplot(.,aes(max.inter))+
  geom_density()+
  facet_wrap(max.inter.res~.,scales='free')

TAD_cl_inter_tbl<-TAD_cl_inter_tbl %>% 
  mutate(max.GRange=pmap(list(inter.cl.GRange,inter.stat),function(inter.cl.GRange,inter.stat){
    inter.cl.GRange[[which.max(inter.stat)]]
  })) 

max_cl_GRange_l<-GRangesList(TAD_cl_inter_summary_tbl$max.GRange)
TAD_Grange_l<-GRangesList(TAD_cl_inter_summary_tbl$GRange)
plot(unlist(lapply(width(IRanges::intersect(TAD_Grange_l,max_cl_GRange_l)),sum))/as.numeric(width(TAD_Grange_l)),TAD_cl_inter_summary_tbl$max.inter)
