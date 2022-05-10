library(tidyverse)
library(GenomicRanges)
library(furrr)
library(mgcv)
library(vroom)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
options(future.globals.maxSize=20971520000)

#--------------------------------
##Utils. Fn

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.na(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo,res_num){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat,cluster = 10)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}  

#--------------------------------

hub_file<-"~/data_transfer/candidate_trans_DAGGER_hub/GM12878_union_trans_res_dagger_tbl.Rda"
res_file<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/spec_res/"
TAD_file<-"~/data_transfer/candidate_trans_DAGGER_hub/TAD_tbl/GM12878_top_TAD_tbl.Rda"

hub_tbl<-get_tbl_in_fn(hub_file)
top_TAD_tbl<-get_tbl_in_fn(TAD_file)

hub_chr_set<-unique(hub_tbl$chr)

plan(multisession,workers=10)
hub_tbl<-do.call(bind_rows,future_map(hub_chr_set,function(chromo){
  chr_spec_res<-get_tbl_in_fn(paste0(res_file,chromo,"_spec_res.Rda"))
  tmp_hub_set<-hub_tbl %>% 
    filter(chr==chromo) %>% 
    dplyr::select(node) %>% unlist
  return(tibble(chr=chromo,node=tmp_hub_set,bins=chr_spec_res$cl_member[tmp_hub_set]) %>% 
           mutate(res=str_split_fixed(node,"_",2)[,1]))
  
}))
plan(sequential)

plan(multisession,workers=10)
hub_tbl<-hub_tbl %>% 
  mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    IRanges::reduce(GRanges(seqnames=chr,
                            ranges = IRanges(start=as.numeric(bins),
                                             end=as.numeric(bins) + res_num[res] -1
                            )))
    
  }))
plan(sequential)
hub_Grange<-GRangesList(hub_tbl$GRange)

top_TAD_GRangeL<-GRangesList(top_TAD_tbl$GRange)
TAD_ID<-top_TAD_tbl %>% 
  mutate(ID=paste(chr,start,end,sep="_")) %>% 
  dplyr::select(ID) %>% unlist
hub_ID<-hub_tbl %>% 
  mutate(ID=paste(chr,node,sep="_")) %>% 
  dplyr::select(ID) %>% unlist

hub_inter_TAD_tbl<-findOverlaps(top_TAD_GRangeL,hub_Grange) %>% 
  as_tibble %>% 
  mutate(TAD.ID=TAD_ID[queryHits],hub.ID=hub_ID[subjectHits]) %>% 
  group_by(hub.ID) %>% 
  summarise(TAD.content=list(unique(TAD.ID)),n.TAD=n()) %>% 
  filter(n.TAD>1) %>% 
  arrange(desc(n.TAD))

#---------------------------------------------------------------
tmp_res_set<-c("50kb","10kb","5kb")
HiC_dat_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/"



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

plan(multisession,workers=10)
inter_TAD_hic_dat_l<-future_map(1:nrow(TAD_pair_combo),function(i){
  
  bins_a<-chr_TAD_bin_content %>% 
    filter(TAD.ID==TAD_pair_combo[i,1]) %>% 
    unnest(cols=c(bin.content)) %>% 
    dplyr::select(bin.content) %>% unlist
  
  bins_b<-chr_TAD_bin_content %>% 
    filter(TAD.ID==TAD_pair_combo[i,2]) %>% 
    unnest(cols=c(bin.content)) %>% 
    dplyr::select(bin.content) %>% unlist
  
  return(inter_TAD_hic_dat %>% 
           filter(X1 %in% bins_a & X2 %in% bins_b | X1 %in% bins_b & X2 %in% bins_a) %>% 
           dplyr::select(X1,X2,X3,zscore))
  
})
plan(sequential)

TAD_pair_combo<-TAD_pair_combo %>% 
  as_tibble %>% 
  mutate(hic.dat=inter_TAD_hic_dat_l)

chr_hub_tbl<-hub_tbl %>% 
  filter(chr==chromo)
chr_hub_Grange<-GRangesList(chr_hub_tbl$GRange)

hub_TAD_content<-findOverlaps(chr_hub_Grange,top_TAD_GRangeL) %>% 
  as_tibble %>% 
  mutate(hub=chr_hub_tbl$node[queryHits],TAD.ID=TAD_ID[subjectHits]) %>% 
  group_by(hub) %>% 
  summarise(TAD.content=list(unique(TAD.ID)),n.TAD=n()) %>% 
  filter(n.TAD>1)

plan(multisession,workers=10)
hub_subset_tbl<-do.call(bind_rows,future_map(1:nrow(hub_TAD_content),function(i){
  tmp_TAD<- hub_TAD_content$TAD.content[[i]]
  TAD_pair_combo %>% 
    filter(V1 %in% tmp_TAD & V2 %in% tmp_TAD | V2 %in% tmp_TAD & V1 %in% tmp_TAD) %>% 
    dplyr::select(hic.dat) %>% 
    unnest(cols=c(hic.dat)) %>% 
    distinct
})) %>% 
  mutate(hub.io="hub")
plan(sequential)

TAD_in_hub<-findOverlaps(chr_hub_Grange,top_TAD_GRangeL) %>% 
  as_tibble %>% 
  mutate(hub=chr_hub_tbl$node[queryHits],TAD.ID=TAD_ID[subjectHits]) %>% 
  distinct(TAD.ID) %>% 
  unlist

gg_dens<-TAD_pair_combo %>% 
  filter(!(V1 %in% TAD_in_hub | V2 %in% TAD_in_hub )) %>% 
  dplyr::select(hic.dat) %>% 
  unnest(cols=c(hic.dat)) %>% 
  distinct %>%
  mutate(hub.io="out") %>% 
  bind_rows(.,hub_subset_tbl) %>% 
  ggplot(.,aes(zscore,color=hub.io))+
  geom_density()
