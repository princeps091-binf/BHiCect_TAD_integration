library(GenomicRanges)
library(tidyverse)
library(furrr)
library(optparse)

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

option_list <- list(
  make_option(c("-t", "--tad"),
              type = "character", default = NULL,
              help = "The TAD files.", metavar = "character"
  ),
  make_option(c("-c", "--cluster"),
              type = "character", default = NULL,
              help = "The cluster folder", metavar = "character"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = "The outputfile", metavar = "character"
  ),
  make_option(c("-n", "--workers"),
              type = "character", default = NULL,
              help = "The number of workers", metavar = "integer"
  )
)

#-----------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

TAD_file<-opt$tad
BHiCect_res_folder<-opt$cluster
nworker<-opt$workers
out_file<-opt$output

TAD_tbl<-read_delim(TAD_file) %>% 
  mutate(chr=paste0("chr",chr1)) %>% 
  dplyr::select(chr,x1,x2) %>% 
  dplyr::rename(start=x1,end=x2)

chr_set<-str_split_fixed(grep("^chr",list.files(BHiCect_res_folder),value=T),pattern = "_",n=2)[,1]
TAD_cl_inter_summary_l<-map(chr_set,function(chromo){
  message(chromo)
  
  tmp_spec_res<-data_tbl_load_fn(paste0(BHiCect_res_folder,chromo,"_spec_res.Rda"))
  tmp_res<-str_split_fixed(names(tmp_spec_res$cl_member),pattern = "_",2)[,1]
  cl_tbl<-tibble(chr=chromo,res=tmp_res,cl=names(tmp_spec_res$cl_member),bins=lapply(tmp_spec_res$cl_member,as.numeric)) %>% 
    mutate(ends=pmap(list(res,bins),function(res,bins){
      bins + res_num[res]-1
    }))
  
  message(chromo,": BHiCect GRange build")
  
  plan(multisession,workers=nworker)
  cl_tbl<-cl_tbl %>% 
    mutate(GRange=future_pmap(list(chr,bins,ends),function(chr,bins,ends){
      reduce(GRanges(seqnames=chr,
                     ranges = IRanges(start=bins,
                                      end=ends
                     )))
    }))
  plan(sequential)
  message(chromo,": TAD GRange build")
  
  plan(multisession,workers=nworker)
  chr_TAD_tbl<-TAD_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(GRange=future_pmap(list(chr,start,end),function(chr,start,end){
      reduce(GRanges(seqnames=chr,
                     ranges = IRanges(start=start,
                                      end=end
                     )))
    }))
  plan(sequential)
  message(chromo,": Build GRangeList")
  
  cl_GRange_l<-GRangesList(cl_tbl$GRange)
  names(cl_GRange_l)<-cl_tbl$cl
  
  chr_TAD_tbl<-chr_TAD_tbl %>% 
    mutate(inter.cl.GRange=map(GRange,function(x){
      cl_tbl$GRange[which(countOverlaps(cl_GRange_l,x)>0)]
    }),inter.cl=map(GRange,function(x){
      cl_tbl$cl[which(countOverlaps(cl_GRange_l,x)>0)]
    }))
  message(chromo,": TAD-BHiCect intersection")
  
  plan(multisession,workers=nworker)
  chr_TAD_tbl<-chr_TAD_tbl %>% 
    mutate(inter.stat=future_pmap(list(GRange,inter.cl.GRange),function(GRange,inter.cl.GRange){
      tmp_l<-GRangesList(lapply(1:length(inter.cl.GRange),function(x)return(GRange)))
      tmp_inter_l<-GRangesList(inter.cl.GRange)
      return(unlist(lapply(width(IRanges::intersect(tmp_l,tmp_inter_l)),sum))/sum(width(IRanges::union(tmp_l,tmp_inter_l))))
    }))
  plan(sequential)
  chr_TAD_tbl<-chr_TAD_tbl %>% 
    mutate(max.inter=map_dbl(inter.stat,max),max.cl=pmap_chr(list(inter.cl,inter.stat),function(inter.cl,inter.stat){
      inter.cl[which.max(inter.stat)]
    }))
  
  return(chr_TAD_tbl)
  
})
TAD_cl_inter_summary_tbl<-do.call(bind_rows,TAD_cl_inter_summary_l)

save(TAD_cl_inter_summary_tbl,file = out_file)
