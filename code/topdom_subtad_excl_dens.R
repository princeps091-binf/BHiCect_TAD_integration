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
jaccard_res_file<-"./data/TopDom_data/HMEC/TopDom_res/HMEC_topdom_tads_jaccard.Rda"
tad_file<-"./data/TopDom_data/HMEC/TopDom_res/topdom_tads.Rda"
#--------------------------------
topdom_tad_tbl<-get_tbl_in_fn(tad_file) %>% 
  filter(name=="domain") %>% 
  mutate(TAD.ID=paste(chrom,chromStart,chromEnd,sep="_"))

chr_res_tbl<-get_tbl_in_fn(jaccard_res_file)

chr_res_tbl<-chr_res_tbl %>% 
  filter(TAD.ID %in% topdom_tad_tbl$TAD.ID)

chr_res_tbl %>% 
  filter(inter < 1) %>% 
  group_by(chr,cl) %>% 
  summarise(n.TAD=n()) %>% 
  ggplot(.,aes(n.TAD))+
  geom_density()+
  scale_x_log10()
