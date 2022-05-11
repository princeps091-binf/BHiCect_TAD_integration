library(tidyverse)
library(GenomicRanges)
library(furrr)
library(data.tree)
library(igraph)
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
res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/chr1_spec_res.Rda"
#HiC_dat_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/"
HiC_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"

TAD_cl_inter_tbl<-get_tbl_in_fn(inter_file)

TAD_best_cl_fit_tbl<-TAD_cl_inter_tbl %>% 
  dplyr::select(chr,start,end,GRange,max.inter,max.cl)
rm(TAD_cl_inter_tbl)

chromo<-"chr22"

chr_TAD_best_cl_fit_tbl<-TAD_best_cl_fit_tbl %>% 
  filter(chr==chromo)

chr_spec_res<-get_tbl_in_fn(paste0("~/Documents/multires_bhicect/data/GM12878/spec_res/",chromo,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
node_ancestor<-lapply(node_ancestor,'[',-1)

TAD_node<-unique(chr_TAD_best_cl_fit_tbl%>%dplyr::select(max.cl) %>% unlist)
TAD_set<-unique(c(TAD_node,unique(unlist(node_ancestor[TAD_node])))) 
Prune(chr_bpt, function(x) x$name %in% TAD_set)

g_bpt<-as.igraph.Node(chr_bpt,directed = T,direction = 'climb')

tmp_d<-distances(g_bpt,TAD_node,TAD_node)

TAD_cl_combo<-t(combn(TAD_node,2)) 
TAD_cl_combo_tbl<-TAD_cl_combo %>% 
  as_tibble %>% 
  mutate(bpt.dist=tmp_d[TAD_cl_combo])

chr_TAD_best_cl_fit_tbl %>% 
  mutate(TAD.ID=paste(chr,start,end,sep="_"))


tmp_res_set<-c("50kb","10kb","5kb")



chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)
tmp_bins<-unique(c(chr_dat$X1,chr_dat$X2))
chr_bin_GRange<-GRanges(seqnames=chromo,
                        ranges = IRanges(start=tmp_bins,
                                         end=tmp_bins + res_num[tmp_res] -1
                        ))

bin_TAD_inter_tbl<-findOverlaps(top_TAD_GRangeL,chr_bin_GRange) %>% 
  as_tibble %>% 
  mutate(TAD.ID=TAD_ID[queryHits],bin=tmp_bins[subjectHits]) 
intersect_size<-unlist(lapply(pintersect(top_TAD_GRangeL[bin_TAD_inter_tbl$queryHits],chr_bin_GRange[bin_TAD_inter_tbl$subjectHits]),width))

bin_TAD_inter_tbl<-bin_TAD_inter_tbl %>% 
  mutate(inter.size=intersect_size)

chr_TAD_bin_content<-bin_TAD_inter_tbl %>% 
  filter(inter.size>=res_num[tmp_res]) %>% 
  group_by(TAD.ID) %>% 
  summarise(bin.content=list(unique(bin)))

TAD_pair_combo<-t(combn(chr_TAD_bin_content$TAD.ID,2))

inter_TAD_hic_dat<-chr_dat %>% 
  filter(X1 %in% unique(unlist(chr_TAD_bin_content$bin.content)) & X2 %in% unique(unlist(chr_TAD_bin_content$bin.content)))
