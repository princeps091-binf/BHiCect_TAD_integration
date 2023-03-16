library(GenomicRanges)
library(tidyverse)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#---------------------------------------------
topdom_res_file<-"./data/TopDom_data/HMEC/TopDom_res/topdom_tads.Rda"
tad_file<-"~/Documents/multires_bhicect/data/HMEC/GSE63525_HMEC_Arrowhead_domainlist.txt"
#---------------------------------------------
topdom_tbl<-get_tbl_in_fn(topdom_res_file)
arrow_tbl<-read_tsv(tad_file)

topdom_Grange<-GRanges(seqnames=topdom_tbl$chrom,
        ranges = IRanges(start=topdom_tbl$chromStart,
                         end=topdom_tbl$chromEnd-1
        ))

arrow_Grange<-GRanges(seqnames=paste0('chr',arrow_tbl$chr1),
                       ranges = IRanges(start=arrow_tbl$x1,
                                        end=arrow_tbl$x2-1
                       ))

topdom_arrow_inter_tbl <- findOverlaps(topdom_Grange,arrow_Grange) %>% 
  as_tibble() %>% 
  mutate(inter.w=width(pintersect(topdom_Grange[queryHits],arrow_Grange[subjectHits])),
         union.w=width(punion(topdom_Grange[queryHits],arrow_Grange[subjectHits]))) %>% 
  mutate(jaccard=inter.w/union.w)

topdom_arrow_inter_tbl %>% 
  group_by(queryHits) %>% 
  slice_max(jaccard) %>% 
  ggplot(.,aes(jaccard))+
  geom_density()

topdom_arrow_inter_tbl %>% 
  group_by(queryHits) %>% 
  slice_max(jaccard) %>% 
  ungroup() %>% 
  summarise(m=median(jaccard))

topdom_tbl %>% 
  mutate(w=chromEnd-chromStart) %>% 
  mutate(set="TopDom") %>% 
  dplyr::select(set,w) %>% 
  bind_rows(.,
            arrow_tbl %>% 
              mutate(w=end-start,
                     set="Arrowhead") %>% 
              dplyr::select(w,set)) %>% 
  ggplot(.,aes(w,color=set))+
  geom_density()+
  scale_color_brewer(palette="Dark2")

topdom_tbl %>% 
  group_by(chrom) %>% 
  mutate(chr.n=1:n()) %>% 
  mutate(set="TopDom") %>%
  filter(chrom=="chr1") %>% 
  dplyr::select(chrom,chromStart,chromEnd,ID,chr.n,set) %>% 
  bind_rows(.,arrow_tbl %>% 
              dplyr::rename(
                chrom=chr,
                chromStart=start,
                chromEnd=end
              ) %>% 
              mutate(set="Arrowhead") %>%
              group_by(chrom) %>% 
              arrange(chromStart) %>% 
              mutate(chr.n=1:n()) %>% 
              filter(chrom=="chr1") %>%
              dplyr::select(chrom,chromStart,chromEnd,ID,chr.n,set)) %>% 
  ggplot(.,aes(x=chromStart,xend=chromEnd,y=chr.n,yend=chr.n,group=ID,color=ID))+
  geom_segment()+facet_grid(set~chrom)+
  theme(legend.position = "none")  

