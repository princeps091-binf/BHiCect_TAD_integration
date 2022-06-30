library(GenomicRanges)
library(tidyverse)
library(vroom)
library(furrr)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
options(scipen = 999999999)

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#---------------------------------------------
tad_file<-"~/Documents/multires_bhicect/data/H1/Dekker/TAD.txt"

hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/H1_union_trans_res_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"
# Build GRange object for TADs and Clusters
tad_tbl<-vroom(tad_file,col_select = c(1:3)) %>% 
  mutate(chr1=paste0("chr",chr1)) %>% 
  dplyr::rename(chr=chr1,start=x1,end=x2)
hub_tbl<-get_tbl_in_fn(hub_file)

## Loop through chromosome
chr_set<-unique(hub_tbl$chr)

res_l<-vector('list',length(chr_set))
names(res_l)<-chr_set
for(chromo in chr_set){
  message(chromo)
  chr_tad_tbl<-tad_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(GRange=pmap(list(chr,start,end),function(chr,start,end){
      IRanges::reduce(GRanges(seqnames=chr,
                              ranges = IRanges(start=start,
                                               end=end-1
                              )))
      
    })) %>% 
    mutate(ID=paste(chr,start,end,sep="_"))
  
  chr_hub_tbl<-hub_tbl %>% 
    filter(chr==chromo)
  chr_spec_res<-get_tbl_in_fn(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  
  chr_hub_tbl<-chr_hub_tbl %>% 
    mutate(bins=chr_spec_res$cl_member[node]) %>%
    mutate(bins=map(bins,as.numeric)) %>% 
    mutate(ends=pmap(list(res,bins),function(res,bins){
      bins + res_num[res] - 1
    })) %>% 
    mutate(GRange=pmap(list(chr,bins,ends),function(chr,bins,ends){
      IRanges::reduce(GRanges(seqnames=chr,
                              ranges = IRanges(start=bins,
                                               end=ends
                              )))
      
    }))
  message("Build Grange: ",chromo)
  chr_TAD_GRangeL<-GRangesList(chr_tad_tbl$GRange)
  chr_cl_GRangeL<-GRangesList(chr_hub_tbl$GRange)
  chr_tad_w<-unlist(width(chr_TAD_GRangeL))
  chr_cl_w<-unlist(lapply(width(chr_cl_GRangeL),sum))
  
  message("Compute intersection: ",chromo)
  plan(multisession,workers=5)
  res_l[[chromo]]<-findOverlaps(chr_TAD_GRangeL,chr_cl_GRangeL) %>% 
    as_tibble %>% 
    mutate(inter.w=future_pmap_dbl(list(queryHits,subjectHits),function(queryHits,subjectHits){
      sum(width(IRanges::intersect(chr_TAD_GRangeL[[queryHits]],chr_cl_GRangeL[[subjectHits]])))
    })) %>% 
    mutate(tad.w=chr_tad_w[queryHits],
           cl.w=chr_cl_w[subjectHits],
           TAD.id=chr_tad_tbl$ID[queryHits],
           cl.id= chr_hub_tbl$node[subjectHits],
           chr=chromo) %>% 
    dplyr::select(chr,TAD.id,cl.id,tad.w,cl.w,inter.w)
  plan(sequential)
  
}

res_tbl<-do.call(bind_rows,res_l)

res_tbl %>% 
  filter(cl.w<tad.w) %>% 
  group_by(cl.id) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(n))+
  geom_density()

res_tbl %>% 
  filter(cl.w>tad.w) %>% 
  mutate(n=inter.w/tad.w) %>% 
  ggplot(.,aes(n))+
  geom_density()

unique(res_tbl$TAD.id)

tad_tbl %>% 
  mutate(ID=paste(chr,start,end,sep="_")) %>% 
  mutate(io=ifelse(ID %in% res_tbl$TAD.id,"in","out")) %>% 
  group_by(io) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(x="TADs",n,fill=io))+
  geom_bar(stat="identity",position="fill")

