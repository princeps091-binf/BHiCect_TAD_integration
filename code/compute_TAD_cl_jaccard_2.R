library(tidyverse)
library(readr)
library(purrr)
library(GenomicRanges)
library(furrr)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6L,5e5L,1e5L,5e4L,1e4L,5e3L)
names(res_num)<-res_set

#--------------------------
# Utils Fn.
get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

spec_res_to_tbl_fn<-function(chr_spec_res,chromo){
  tibble(chr=chromo,cl=names(chr_spec_res$cl_member),bins=map(chr_spec_res$cl_member,as.integer))
}

bin_to_GRange<-function(chr,cl,bins){
  tmp_res<-str_split_fixed(cl,"_",2)[,1]
  return(GRanges(seqnames=chr,
                          ranges = IRanges(start=bins,
                                           end=bins + res_num[tmp_res] -1
                          )))
  
}
#--------------------------
# data files
## TAD
TAD_folder<-"~/Documents/multires_bhicect/data/HMEC/GSE63525_HMEC_Arrowhead_domainlist.txt"
## BHiCect
BHiCect_cl_folder<-"~/Documents/multires_bhicect/data/HMEC/determinate_spec_res/"
#--------------------------
# Load data
tad_tbl<-readr::read_delim(TAD_folder,delim = '\t')
## Loop genome-wide
chromo<-'chr1'
chr_spec_res<-get_tbl_in_fn(paste0(BHiCect_cl_folder,chromo,"_spec_res.Rda"))
chr_cl_tbl<-spec_res_to_tbl_fn(chr_spec_res,chromo)
#--------------------------
# Convert to GRange object
plan(multisession, workers=4)
chr_cl_tbl<-chr_cl_tbl %>% 
  mutate(Grange=future_pmap(list(chr,cl,bins),function(chr,cl,bins){
    bin_to_GRange(chr,cl,bins)
  }))
plan(sequential)
chr_tad_tbl<-tad_tbl %>% 
  mutate(chr=paste0('chr',chr1)) %>% 
  filter(chr==chromo) %>% 
  dplyr::select(chr,x1,x2) %>% 
  dplyr::mutate(Grange=pmap(list(chr,x1,x2),function(chr,x1,x2){
    return(GRanges(seqnames=chr,
                   ranges = IRanges(start=x1,
                                    end=x2-1
                   )))
    
  }))

#--------------------------
# Intersect and find best fit
tad_cl_inter_tbl<-as_tibble(findOverlaps(GRangesList(chr_tad_tbl$Grange),GRangesList(chr_cl_tbl$Grange)))
# Compute Jaccard
tmp_ego<-GRangesList(chr_tad_tbl$Grange[tad_cl_inter_tbl$queryHits])
tmp_alter<-GRangesList(chr_cl_tbl$Grange[tad_cl_inter_tbl$subjectHits])
tad_cl_inter_tbl<-tad_cl_inter_tbl %>% 
  mutate(inter.w=sum(width(intersect(tmp_ego,tmp_alter))),
         union.w=sum(width(union(tmp_ego,tmp_alter))),
         tad.w=sum(width(tmp_ego)),
         cl.w=sum(width(tmp_alter))) %>% 
  mutate(jaccard=inter.w/union.w)

#---------------------------------------------------------
# Recovery of TADs
tad_cl_inter_tbl1 %>% 
  group_by(queryHits) %>% 
  slice_max(jaccard) %>% 
  ggplot(.,aes(jaccard))+
  geom_density()

# Recovery of sub-TADs
tad_cl_inter_tbl1 %>% 
  filter(cl.w < tad.w) %>% 
  group_by(subjectHits) %>% 
  count %>%
  ggplot(.,aes(n))+
  geom_density()