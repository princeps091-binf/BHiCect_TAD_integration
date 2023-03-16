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
topdom_res_file<-"./data/TopDom_data/HMEC/TopDom_res/topdom_tads.Rda"
BHiCect_cl_folder<-"~/Documents/multires_bhicect/data/HMEC/determinate_spec_res/"
#--------------------------
# Load data
tad_tbl<-get_tbl_in_fn(topdom_res_file)
#--------------------------
#to loop genome-wide
chromo<-'chr3'

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

plan(multisession, workers=4)
chr_tad_tbl<-tad_tbl %>% 
  filter(chrom == chromo) %>% 
  filter(name=="domain") %>% 
  mutate(ID=paste(chrom,1:n(),"topdom",sep="_")) %>% 
  mutate(Grange=future_pmap(list(chrom, chromStart, chromEnd),function(chrom,chromStart,chromEnd){
    IRanges::reduce(GRanges(seqnames=chrom,
                            ranges = IRanges(start=chromStart,
                                             end=chromEnd-1
                            )))
  }))
plan(sequential)

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
  mutate(jaccard=inter.w/union.w) %>% 
  mutate(tad_ID=paste0("tad_",queryHits,"_",chromo),
         cl_ID=paste0("cl_",subjectHits,"_",chromo)) %>% 
  dplyr::select(-c(queryHits,subjectHits))
#---------------------------------------------------------
# Recovery of TADs
tad_cl_inter_tbl %>% 
  group_by(tad_ID) %>% 
  slice_max(jaccard) %>% 
  ggplot(.,aes(jaccard))+
  geom_density()

tad_cl_inter_tbl %>% 
  group_by(tad_ID) %>% 
  slice_max(jaccard) %>% 
  ungroup() %>% 
  summarise(m=median(jaccard))

# Recovery of sub-TADs
tad_cl_inter_tbl %>% 
  filter(cl.w < tad.w) %>% 
  group_by(subjectHits) %>% 
  count %>%
  ggplot(.,aes(n))+
  geom_density()
