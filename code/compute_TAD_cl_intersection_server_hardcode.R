library(GenomicRanges)
library(tidyverse)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
# Utils Fn.
data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}


#-----------------------------------------
TAD_file<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
BHiCect_res_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/spec_res/"
out_file<-"~/data_transfer/GM12878_TAD_cl_intersect.Rda"
#-----------------------------------------
TAD_tbl<-read_delim(TAD_file,delim = "\t",col_names = T) %>% 
  mutate(chr=paste0("chr",chr1)) %>% 
  dplyr::select(chr,x1,x2) %>% 
  dplyr::rename(start=x1,end=x2)
nworker<-15
chr_set<-str_split_fixed(grep("^chr",list.files(BHiCect_res_folder),value=T),pattern = "_",n=2)[,1]

TAD_cl_inter_summary_l<-vector('list',length(chr_set))
names(TAD_cl_inter_summary_l)<-chr_set
for(chromo in chr_set){
  message(chromo)
  
  tmp_spec_res<-data_tbl_load_fn(paste0(BHiCect_res_folder,chromo,"_spec_res.Rda"))
  tmp_res<-str_split_fixed(names(tmp_spec_res$cl_member),pattern = "_",2)[,1]
  cl_tbl<-tibble(chr=chromo,res=tmp_res,cl=names(tmp_spec_res$cl_member),bins=lapply(tmp_spec_res$cl_member,as.numeric)) %>% 
    mutate(ends=pmap(list(res,bins),function(res,bins){
      bins + res_num[res]-1
    }))
  
  message(chromo,": BHiCect GRange build")
  tmp_bins<-cl_tbl$bins
  tmp_ends<-cl_tbl$ends
  tmp_chr<-cl_tbl$chr
  
  cl <- makeCluster(nworker)
  clusterEvalQ(cl, {
    library(dplyr)
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl, c("tmp_bins", "tmp_ends", "tmp_chr"))
  
  tmp_Grange <- parLapply(cl, 1:length(tmp_bins), function(x) {
    reduce(GRanges(seqnames=tmp_chr[x],
                   ranges = IRanges(start=tmp_bins[[x]],
                                    end=tmp_ends[[x]]
                   )))
    
  })
  stopCluster(cl)
  rm(cl)
  
  cl_tbl<-cl_tbl %>% 
    mutate(GRange=tmp_Grange)
  rm(tmp_Grange,tmp_ends,tmp_bins,tmp_chr)
  
  message(chromo,": TAD GRange build")
  
  chr_TAD_tbl<-TAD_tbl %>% 
    filter(chr==chromo) 
  
  tmp_bins<-chr_TAD_tbl$start
  tmp_ends<-chr_TAD_tbl$end
  tmp_chr<-chr_TAD_tbl$chr
  
  cl <- makeCluster(nworker)
  clusterEvalQ(cl, {
    library(dplyr)
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl, c("tmp_bins", "tmp_ends", "tmp_chr"))
  
  tmp_Grange <- parLapply(cl, 1:length(tmp_bins), function(x) {
    reduce(GRanges(seqnames=tmp_chr[x],
                   ranges = IRanges(start=tmp_bins[x],
                                    end=tmp_ends[x]
                   )))
    
  })
  stopCluster(cl)
  rm(cl)
  
  chr_TAD_tbl<-chr_TAD_tbl %>% 
    mutate(GRange=tmp_Grange)
  rm(tmp_Grange,tmp_ends,tmp_bins,tmp_chr)
  message(chromo,": Build GRangeList")
  
  cl_GRange_l<-GRangesList(cl_tbl$GRange)
  names(cl_GRange_l)<-cl_tbl$cl
  
  message(chromo,": TAD-BHiCect intersection")
  
  chr_TAD_tbl<-chr_TAD_tbl %>% 
    mutate(inter.cl.GRange=map(GRange,function(x){
      cl_tbl$GRange[which(countOverlaps(cl_GRange_l,x)>0)]
    }),inter.cl=map(GRange,function(x){
      cl_tbl$cl[which(countOverlaps(cl_GRange_l,x)>0)]
    }))
  message(chromo,": TAD-BHiCect intersection level metric")
  
  tmp_tad<-chr_TAD_tbl$GRange
  tmp_inter_cl_GRange<-chr_TAD_tbl$inter.cl.GRange
  cl <- makeCluster(nworker)
  clusterEvalQ(cl, {
    library(dplyr)
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl, c("tmp_tad", "tmp_inter_cl_GRange"))
  
  tmp_stat_l <- parLapply(cl, 1:length(tmp_tad), function(x) {
    tmp_l<-GRangesList(lapply(1:length(tmp_inter_cl_GRange[[x]]),function(y)return(tmp_tad[[x]])))
    tmp_inter_l<-GRangesList(tmp_inter_cl_GRange[[x]])
    return(unlist(lapply(width(IRanges::intersect(tmp_l,tmp_inter_l)),sum))/sum(width(IRanges::union(tmp_l,tmp_inter_l))))
    
  })
  stopCluster(cl)
  rm(cl)
  chr_TAD_tbl<-chr_TAD_tbl %>% 
    mutate(inter.stat=tmp_stat_l)
  rm(tmp_stat_l,tmp_tad,tmp_inter_cl_GRange)
  
  
  chr_TAD_tbl<-chr_TAD_tbl %>% 
    mutate(max.inter=map_dbl(inter.stat,max),max.cl=pmap_chr(list(inter.cl,inter.stat),function(inter.cl,inter.stat){
      inter.cl[which.max(inter.stat)]
    }))
  
  TAD_cl_inter_summary_l[[chromo]]<-chr_TAD_tbl
  
}

TAD_cl_inter_summary_tbl<-do.call(bind_rows,TAD_cl_inter_summary_l)

save(TAD_cl_inter_summary_tbl,file = out_file)
