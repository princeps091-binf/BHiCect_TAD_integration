library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
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
jaccard_res_file<-"./data/TopDom_data/HMEC/TopDom_res/HMEC_topdom_tads_jaccard.Rda"
res_file<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/spec_res/"
HiC_dat_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/"

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
## Determine maximum bin-size fitting within smallest TAD
min_TAD_size<-top_fit_tad_cl_tbl %>% 
  mutate(w=map_dbl(TAD.ID,function(x){
    tad_coord<-str_split_fixed(x,"_",3)
    return(as.numeric(tad_coord[3])-as.numeric(tad_coord[2]))
  })) %>% 
  summarise(m=min(w)) %>% 
  unlist

tmp_res<-names(which.max(res_num[which(res_num<=min_TAD_size)]))

chr_set<-grep("^chr",str_split_fixed(list.files(res_file),"_",2)[,1],value=T)

chr_res_l<-vector("list",length(chr_set))
names(chr_res_l)<-chr_set

for(chromo in chr_set){
  message(chromo)
  chr_TAD_best_cl_fit_tbl<-top_fit_tad_cl_tbl %>% 
    filter(chr==chromo)

  message("compute BPT distance for ", chromo)
  chr_spec_res<-get_tbl_in_fn(paste0(res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  
  TAD_node<-unique(chr_TAD_best_cl_fit_tbl%>%dplyr::select(cl) %>% unlist)
  TAD_set<-unique(c(TAD_node,unique(unlist(node_ancestor[TAD_node])))) 
  Prune(chr_bpt, function(x) x$name %in% TAD_set)
  
  g_bpt<-as.igraph.Node(chr_bpt,directed = T,direction = 'climb')
  
  tmp_d<-distances(g_bpt,TAD_node,TAD_node)
  
  
  message("compute zscore for ", chromo)
  chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)
  tmp_bins<-unique(c(chr_dat$X1,chr_dat$X2))
  
  chr_bin_GRange<-GRanges(seqnames=chromo,
                          ranges = IRanges(start=tmp_bins,
                                           end=tmp_bins + res_num[tmp_res] -1
                          ))
  TAD_ID<-chr_TAD_best_cl_fit_tbl$TAD.ID
  top_TAD_GRangeL<-GRangesList(chr_TAD_best_cl_fit_tbl$GRange)
  
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
  
  message("Subset inter-TAD interaction for", chromo)
  
  plan(multisession,workers=15)
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
  
  best_cl_vec<-chr_TAD_best_cl_fit_tbl$max.cl
  names(best_cl_vec)<-chr_TAD_best_cl_fit_tbl$TAD.ID
  
  TAD_pair_combo_tbl<-TAD_pair_combo %>% 
    as_tibble %>% 
    mutate(cl.A=best_cl_vec[TAD_pair_combo[,1]],cl.B=best_cl_vec[TAD_pair_combo[,2]])
  
  TAD_pair_combo_tbl<-TAD_pair_combo_tbl %>% 
    mutate(bpt.d=tmp_d[as.matrix(TAD_pair_combo_tbl[,3:4])])
  
  TAD_pair_combo_tbl<-TAD_pair_combo_tbl %>% 
    mutate(med.z=inter_TAD_hic_dat_l)
  
  chr_res_l[[chromo]]<-TAD_pair_combo_tbl
}
