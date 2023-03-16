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
top_fit_tad_cl_tbl<-chr_res_tbl %>% 
  group_by(TAD.ID) %>% 
  slice_max(jaccard) %>% 
  ungroup

plan(multisession, workers=5)

top_fit_tad_cl_tbl<-top_fit_tad_cl_tbl %>%
    mutate(inter=future_pmap_int(list(TAD.ID,cl.start,res,chr),function(TAD.ID,cl.start,res,chr){
      cl_Grange<- IRanges::reduce(GRanges(seqnames=chr,
                                          ranges = IRanges(start=as.numeric(cl.start),
                                                           end=as.numeric(cl.start) + res_num[res] -1
                                          )))
      tad_split<-str_split_fixed(TAD.ID,"_",3)
      
      tad_Grange<-GRanges(seqnames=tad_split[1,1],
                          ranges = IRanges(start=as.numeric(tad_split[1,2]),
                                           end=as.numeric(tad_split[1,3])-1
                          ))
      return(sum(width(intersect(cl_Grange,tad_Grange))))
    }))
plan(sequential)

plan(multisession, workers=5)
top_fit_tad_cl_tbl<-top_fit_tad_cl_tbl %>% 
  mutate(tad.inter=future_pmap_dbl(list(TAD.ID,inter),function(TAD.ID,inter){
    
    tad_split<-str_split_fixed(TAD.ID,"_",3)
    
    return(inter/(as.numeric(tad_split[1,3]) - as.numeric(tad_split[1,2])))
  }))
plan(sequential)

plan(multisession, workers=5)
top_fit_tad_cl_tbl<-top_fit_tad_cl_tbl %>% 
  mutate(cl.inter=future_pmap_dbl(list(cl.start,res,chr,inter),function(cl.start,res,chr,inter){
    
    cl_Grange<- IRanges::reduce(GRanges(seqnames=chr,
                                        ranges = IRanges(start=as.numeric(cl.start),
                                                         end=as.numeric(cl.start) + res_num[res] -1
                                        )))
    
    return(inter/sum(width(cl_Grange)))
  }))
plan(sequential)

top_fit_tad_cl_tbl %>% 
  ggplot(.,aes(tad.inter,cl.inter))+
  geom_density_2d_filled()
