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
TAD_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_GM12878_TAD_pval_tbl.Rda"
top_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda"

hub_tbl<-get_tbl_in_fn(top_hub_file)

TAD_tbl<-get_tbl_in_fn(TAD_file)%>% 
  group_by(chr) %>% 
  mutate(FDR=p.adjust(emp.pval,method="fdr")) %>% 
  filter(FDR<=0.01)

TAD_cl_inter_tbl<-get_tbl_in_fn(inter_file)

TAD_cl_inter_tbl %>% 
  ggplot(.,aes(max.inter))+
  geom_histogram()

TAD_cl_inter_tbl %>%
  mutate(ID=paste(chr,start,end,sep="_")) %>% 
  mutate(cage.rich=ifelse(ID %in% TAD_tbl$ID,"rich",'other')) %>% 
  ggplot(.,aes(max.inter,color=cage.rich))+
  geom_density()


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
  geom_density_2d_filled(show.legend = F)
ggsave("~/Documents/multires_bhicect/weeklies/weekly55/img/GM12878_TAD_cl_intersect.png")

TAD_cl_inter_tbl %>% 
  mutate(ID=paste(chr,start,end,sep="_")) %>% 
  mutate(tad.cover=tad_cover,cl.cover=cl_cover) %>% 
  mutate(cage.rich=ifelse(ID %in% TAD_tbl$ID,"rich",'other')) %>% 
  ggplot(.,aes(cl.cover,color=cage.rich))+
  geom_density()
  

TAD_cl_inter_tbl %>% 
  dplyr::select(chr,start,end,inter.cl,inter.stat) %>% 
  mutate(ID=1:n()) %>% 
  unnest(cols=c(inter.cl,inter.stat)) %>% 
  left_join(.,hub_tbl %>% 
              dplyr::select(chr,node) %>% 
              dplyr::rename(inter.cl=node) %>% 
              mutate(io="hub")) %>% 
  mutate(io=ifelse(is.na(io),"out",io)) %>% 
  group_by(ID) %>% 
  summarise(n=sum(io == "hub")) %>% 
  filter(n>0)
  

tmp_tad_stat.tbl %>% 
  ggplot(.,aes(inter.stat,color=io))+
  geom_density()
