library(tidyverse)
library(GenomicRanges)
library(furrr)
library(igraph)
library(vroom)
library(mgcv)
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
TAD_file<-"~/data_transfer/candidate_trans_DAGGER_hub/TAD_tbl/CAGE_union_HMEC_TAD_pval_tbl.Rda"

TAD_tbl<-get_tbl_in_fn(TAD_file)
# Eliminate sub-TADs
TAD_tbl<-TAD_tbl %>% 
  dplyr::select(-GRange) %>% 
  mutate(start=str_split_fixed(ID,"_",3)[,2],
         end=str_split_fixed(ID,"_",3)[,3])
TAD_GRange<-  GRanges(seqnames=TAD_tbl$chr,
                      ranges = IRanges(start=as.numeric(TAD_tbl$start),
                                       end=as.numeric(TAD_tbl$end)-1
                      ))
plan(multisession,workers=10)
TAD_tbl<-TAD_tbl %>%
  mutate(GRange=future_pmap(list(chr,start,end),function(chr,start,end){
    GRanges(seqnames=chr,
            ranges = IRanges(start=as.numeric(start),
                             end=as.numeric(end)-1
            ))
  }))
plan(sequential)


TAD_GRangeL<-GRangesList(TAD_tbl$GRange)

subTAD_tbl<-as_tibble(findOverlaps(TAD_GRangeL,TAD_GRangeL)) %>%
  group_by(queryHits) %>% 
  summarise(tad.inter=list(unique(subjectHits)))

subTAD_inter_tbl<-as_tibble(findOverlaps(TAD_GRangeL,TAD_GRangeL))
detect_TAD_set<-components(graph_from_data_frame(subTAD_inter_tbl))

plan(multisession,workers=4)
TAD_set_GRange<-future_map(1:detect_TAD_set$no,function(x){
  do.call("c",TAD_tbl$GRange[which(detect_TAD_set$membership == x)])
})
plan(sequential)

plan(multisession,workers=10)
top_TAD_GRange<-future_map(TAD_set_GRange,function(x){
  IRanges::reduce(x)
})
plan(sequential)

top_TAD_tbl<-tibble(chr=as.character(unlist(lapply(top_TAD_GRange,function(x)as.character(seqnames(x))))),
                    start=map_int(top_TAD_GRange,start),
                    end=map_int(top_TAD_GRange,end),
                    GRange=top_TAD_GRange)
save(top_TAD_tbl,file="~/data_transfer/candidate_trans_DAGGER_hub/TAD_tbl/HMEC_top_TAD_tbl.Rda")