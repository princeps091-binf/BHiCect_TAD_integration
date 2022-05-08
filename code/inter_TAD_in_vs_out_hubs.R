library(tidyverse)
library(GenomicRanges)
library(furrr)
library(igraph)
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
hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_trans_res_dagger_tbl.Rda"
res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
hub_tbl<-get_tbl_in_fn(hub_file)
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
plan(multisession,workers=4)
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

plan(multisession,workers=4)
top_TAD_GRange<-future_map(TAD_set_GRange,function(x){
  IRanges::reduce(x)
})
plan(sequential)
hist(log10(unlist(lapply(top_TAD_GRange,width))))
top_TAD_tbl<-tibble(chr=as.character(unlist(lapply(top_TAD_GRange,function(x)as.character(seqnames(x))))),
       start=map_int(top_TAD_GRange,start),
       end=map_int(top_TAD_GRange,end),
       GRange=top_TAD_GRange)

hub_chr_set<-unique(hub_tbl$chr)

plan(multisession,workers=4)
hub_tbl<-do.call(bind_rows,future_map(hub_chr_set,function(chromo){
  chr_spec_res<-get_tbl_in_fn(paste0(res_file,chromo,"_spec_res.Rda"))
  tmp_hub_set<-hub_tbl %>% 
    filter(chr==chromo) %>% 
    dplyr::select(node) %>% unlist
  return(tibble(chr=chromo,node=tmp_hub_set,bins=chr_spec_res$cl_member[tmp_hub_set]) %>% 
    mutate(res=str_split_fixed(node,"_",2)[,1]))
  
}))
plan(sequential)

plan(multisession,workers=4)
hub_tbl<-hub_tbl %>% 
  mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    IRanges::reduce(GRanges(seqnames=chr,
                   ranges = IRanges(start=as.numeric(bins),
                                    end=as.numeric(bins) + res_num[res] -1
                   )))
    
  }))
plan(sequential)
hub_Grange<-GRangesList(hub_tbl$GRange)
plan(multisession,workers=4)
TAD_hub_inter<-future_map(TAD_tbl$GRange,function(x){
  unique(subjectHits(findOverlaps(x,hub_Grange)))
})
plan(sequential)
#---------------------------------------------------------------
chromo<-"chr22"
chr_hub_tbl<-hub_tbl %>% 
  filter(chr==chromo)
chr_hub_Grange<-GRangesList(chr_hub_tbl$GRange)
chr_TAD_tbl<-TAD_tbl %>% 
  filter(chr==chromo)

plan(multisession,workers=4)
chr_TAD_hub_inter<-future_map(chr_TAD_tbl$GRange,function(x){
  chr_hub_tbl$node[unique(subjectHits(findOverlaps(x,chr_hub_Grange)))]
})
plan(sequential)

TAD_pair_combo<-t(combn(chr_TAD_tbl$ID,2))

tmp_hub<-chr_hub_tbl$node[23]
chr_TAD_tbl$ID[which(unlist(lapply(chr_TAD_hub_inter,function(x){
  tmp_hub %in% x
})))]
