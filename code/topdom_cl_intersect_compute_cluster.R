library(tidyverse)
library(vroom)
library(GenomicRanges)
library(parallel)
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
BHiCect_res_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/HMEC/HMEC/spec_res/"
TAD_file<-"~/data_transfer/candidate_trans_DAGGER_hub/TAD_tbl/topdom/HMEC_topdom_tads.Rda"
#--------------------------------
TAD_cl_inter_tbl<-get_tbl_in_fn(TAD_file) %>%
  mutate(TAD.ID=paste(chrom,chromStart,chromEnd,sep="_"))

tad_Grange<-   GRanges(seqnames=TAD_cl_inter_tbl$chrom,
                       ranges = IRanges(start=TAD_cl_inter_tbl$chromStart,
                                        end=TAD_cl_inter_tbl$chromEnd - 1
                       ))
mcols(tad_Grange)<-tibble(TAD.ID=TAD_cl_inter_tbl$TAD.ID)

chr_set<-str_split_fixed(grep("^chr",list.files(BHiCect_res_folder),value = T),"_",2)[,1]
chr_res_l<-vector("list",length(chr_set))
names(chr_res_l)<-chr_set

for(chromo in chr_set){
  message(chromo)
  load(paste0(BHiCect_res_folder,chromo,"_spec_res.Rda"))
  chr_cl_tbl<-tibble(chr=chromo,res=str_split_fixed(names(chr_spec_res$cl_member),"_",2)[,1],cl=names(chr_spec_res$cl_member),bin=chr_spec_res$cl_member)
  cl <- makeCluster(20)
  clusterEvalQ(cl, {
    library(dplyr)
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl, c("chr_cl_tbl", "tad_Grange","res_num","chromo"))
  tmp_stat_l <- parLapply(cl, 1:nrow(chr_cl_tbl), function(x) {
    
    cl_Grange<- GRanges(seqnames=chr_cl_tbl$chr[x],
                        ranges = IRanges(start=as.numeric(chr_cl_tbl$bin[[x]]),
                                         end=as.numeric(chr_cl_tbl$bin[[x]]) + res_num[chr_cl_tbl$res[x]] -1
                        ))
    tad_in_idx<-unique(subjectHits(findOverlaps(cl_Grange,tad_Grange)))
    
    inter_tad_Grange_l<-GRangesList(lapply(tad_in_idx,function(i){
      tad_Grange[i]
    }))
    inter_cl_Grange_l<-GRangesList(lapply(tad_in_idx,function(i){
      cl_Grange
    }))
    
    tad_cl_inter<-intersect(inter_tad_Grange_l,inter_cl_Grange_l)
    
    tad_cl_inter_stat<-unlist(lapply(seq_along(tad_in_idx),function(i){
      sum(width(tad_cl_inter[[i]]))/width(tad_Grange[tad_in_idx[i]])
    }))
    
    return(tibble(chr=chromo,cl=chr_cl_tbl$cl[x],TAD.ID=mcols(tad_Grange)$TAD.ID[tad_in_idx],inter=tad_cl_inter_stat))
  })
  stopCluster(cl)
  rm(cl)
  tmp_stat_tbl<-do.call(bind_rows,tmp_stat_l) %>%
    mutate(chr=chromo)
  chr_res_l[[chromo]]<-tmp_stat_tbl
}
chr_res_tbl<-do.call(bind_rows,chr_res_l)
