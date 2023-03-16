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

inter_file<-"~/data_transfer/candidate_trans_DAGGER_hub/TAD_tbl/topdom/HMEC_topdom_tads_inter_tbl.Rda"
spec_res_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/HMEC/HMEC/spec_res/"
#--------------------------------
inter_tbl<-get_tbl_in_fn(inter_file)

inter_tbl %>% 
  mutate(chr=fct_relevel(chr,paste0("chr",1:22))) %>% 
  group_by(chr,TAD.ID) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(n,color=chr))+
  geom_density()+
  scale_color_viridis_d()

chr_set<-unique(inter_tbl$chr)
chr_res_tbl<-do.call(bind_rows,lapply(chr_set,function(chromo){
  message(chromo)
  chr_inter_tbl<-inter_tbl %>% 
    filter(chr==chromo)
  chr_spec_res<-get_tbl_in_fn(paste0(spec_res_folder,chromo,"_spec_res.Rda"))
  chr_inter_tbl<-chr_inter_tbl %>% 
    mutate(cl.start=chr_spec_res$cl_member[cl],
           res=str_split_fixed(cl,"_",2)[,1])
  plan(multisession, workers=5)
  res_tbl<-chr_inter_tbl %>%
    mutate(jaccard=future_pmap_dbl(list(chr,cl.start,res,TAD.ID),function(chr,cl.start,res,TAD.ID){
      cl_Grange<- IRanges::reduce(GRanges(seqnames=chr,
                          ranges = IRanges(start=as.numeric(cl.start),
                                           end=as.numeric(cl.start) + res_num[res] -1
                          )))
      
      tad_split<-str_split_fixed(TAD.ID,"_",3)
      
      tad_Grange<-GRanges(seqnames=tad_split[1,1],
                          ranges = IRanges(start=as.numeric(tad_split[1,2]),
                                           end=as.numeric(tad_split[1,3])-1
                          ))
      return(sum(width(intersect(cl_Grange,tad_Grange)))/sum(width(union(cl_Grange,tad_Grange))))

      
    }))
  plan(sequential)
  return(res_tbl)
}))

res_tbl %>% 
  group_by(TAD.ID) %>% 
  slice_max(inter) %>% 
  arrange(inter)
  ggplot(.,aes(jaccard))+
  geom_density()
save(chr_res_tbl,file="~/data_transfer/candidate_trans_DAGGER_hub/TAD_tbl/topdom/HMEC_topdom_tads_jaccard.Rda")
#---------------------------------------------------------------
jaccard_res_file<-"./data/TopDom_data/HMEC/TopDom_res/HMEC_topdom_tads_jaccard.Rda"
chr_res_tbl<-get_tbl_in_fn(jaccard_res_file)
chr_res_tbl %>% 
  group_by(TAD.ID) %>% 
  slice_max(jaccard) %>% 
  ggplot(.,aes(jaccard))+
  geom_density()
