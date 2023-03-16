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

cage_rich_TAD_tbl<-get_tbl_in_fn(TAD_file)%>% 
  group_by(chr) %>% 
  mutate(FDR=p.adjust(emp.pval,method="fdr")) %>% 
  filter(FDR<=0.01)

TAD_cl_inter_tbl<-get_tbl_in_fn(inter_file)


TAD_hub_inter_tbl<-TAD_cl_inter_tbl %>% 
  mutate(ID=1:n()) %>% 
  dplyr::select(ID,chr,start,end) %>% 
  left_join(.,
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
  )

TAD_hub_inter_tbl<-TAD_hub_inter_tbl %>% 
  mutate(ID2=paste(chr,start,end,sep = "_")) %>% 
  left_join(.,cage_rich_TAD_tbl %>% 
              ungroup %>% 
              mutate(rich="rich") %>% 
              dplyr::select(ID,rich),by=c("ID2"="ID"))

io_mat<-TAD_hub_inter_tbl %>% 
  mutate(all.hub.io=ifelse(is.na(n),"out","in")) %>% 
  group_by(all.hub.io) %>% 
  summarise(n.all=n()) %>% 
  left_join(.,
            TAD_hub_inter_tbl %>% 
              filter(!(is.na(rich))) %>% 
              mutate(all.hub.io=ifelse(is.na(n),"out","in")) %>% 
              group_by(all.hub.io) %>% 
              summarise(n.cage.rich=n())
  ) %>% 
  dplyr::select(-all.hub.io) %>% as.matrix
chisq.test(io_mat[,2],p=io_mat[,1]/sum(io_mat[,1]))$obs
