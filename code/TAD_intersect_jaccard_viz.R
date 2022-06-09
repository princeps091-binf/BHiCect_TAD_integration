library(tidyverse)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#---------------------------------------------

tad_inter_hmec_tbl<-get_tbl_in_fn("~/data_transfer/HMEC_top_TAD_cl_intersect.Rda")
tad_inter_h1_tbl<-get_tbl_in_fn("~/data_transfer/H1_top_TAD_cl_intersect.Rda")
tad_inter_gm12878_tbl<-get_tbl_in_fn("~/data_transfer/GM12878_top_TAD_cl_intersect.Rda")

gm12878_res_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/spec_res/"
hmec_res_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/HMEC/HMEC/spec_res/"
h1_res_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/H1/Dekker/spec_res/"

gg_tmp<-tad_inter_hmec_tbl %>% 
  dplyr::select(chr,start,end,max.inter,max.cl) %>% 
  mutate(cell.line="HMEC") %>% 
  bind_rows(.,
            tad_inter_h1_tbl %>% 
              dplyr::select(chr,start,end,max.inter,max.cl) %>% 
              mutate(cell.line="H1")) %>% 
  bind_rows(.,
            tad_inter_gm12878_tbl %>% 
              dplyr::select(chr,start,end,max.inter,max.cl) %>% 
              mutate(cell.line="GM12878")) %>% 
  ggplot(.,aes(max.inter,color=cell.line))+
  geom_density()

ggsave("~/data_transfer/TAD_inter_mcell_dens.png",gg_tmp)

cl_Grange_build_fn<-function(tad_inter_hmec_tbl,hmec_res_folder){
  
  tmp_chr_set<-unique(tad_inter_hmec_tbl$chr)
  
  out_l<-lapply(tmp_chr_set,function(chromo){
    message(chromo)
    tmp_cl<-tad_inter_hmec_tbl %>% 
      filter(chr==chromo) %>% 
      dplyr::select(max.cl) %>% 
      unlist
    load(paste0(hmec_res_folder,chromo,"_spec_res.Rda"))
    tmp_tbl<-tibble(chr=chromo,max.cl=tmp_cl,bin=chr_spec_res$cl_member[tmp_cl])%>% 
      mutate(res=str_split_fixed(max.cl,"_",2)[,1]) %>% 
      mutate(ends=pmap(list(bin,res),function(bin,res){
        as.numeric(bin) + res_num[res]-1
      })) %>% 
      mutate(max.GRange=pmap(list(chr,bin,ends),function(chr,bin,ends){
        IRanges::reduce(GRanges(seqnames=chr,
                       ranges = IRanges(start=as.numeric(bin),
                                        end=ends
                       )))
        
      }))
    
    return(chr_tad_inter_hmec_tbl<-tad_inter_hmec_tbl %>% 
      filter(chr==chromo) %>% 
      dplyr::select(chr,start,end,GRange,max.inter,max.cl) %>% 
      left_join(.,tmp_tbl %>% 
                  dplyr::select(chr,max.cl,max.GRange)))
  })
  
  return(do.call(bind_rows,out_l))
}
tad_cl_inter_fn<-function(tad_inter_hmec_tbl,cell_line){
  
  max_cl_GRange_l<-GRangesList(tad_inter_hmec_tbl$max.GRange)
  TAD_Grange_l<-GRangesList(tad_inter_hmec_tbl$GRange)
  tad_cover<-unlist(lapply(width(IRanges::intersect(TAD_Grange_l,max_cl_GRange_l)),sum))/(as.numeric(width(TAD_Grange_l)))
  cl_cover<-unlist(lapply(width(IRanges::intersect(TAD_Grange_l,max_cl_GRange_l)),sum))/as.numeric(unlist(lapply(width(max_cl_GRange_l),sum)))
  
  return(tibble(tad=tad_cover,cl=cl_cover,cell.line=cell_line))
  
}

tad_max_inter_h1_tbl<-cl_Grange_build_fn(tad_inter_h1_tbl,h1_res_folder) %>% 
  tad_cl_inter_fn(.,"H1")
tad_max_inter_hmec_tbl<-cl_Grange_build_fn(tad_inter_hmec_tbl,hmec_res_folder) %>% 
  tad_cl_inter_fn(.,"HMEC")
tad_max_inter_gm12878_tbl<-cl_Grange_build_fn(tad_inter_gm12878_tbl,gm12878_res_folder) %>% 
  tad_cl_inter_fn(.,"GM12878")

gg_tmp<-tad_max_inter_h1_tbl %>%
  bind_rows(.,tad_max_inter_hmec_tbl) %>% 
  bind_rows(.,tad_max_inter_gm12878_tbl) %>% 
  ggplot(.,aes(tad,cl))+
  geom_density_2d_filled(show.legend = F)+facet_wrap(cell.line~.)
ggsave("~/data_transfer/tad_cl_inter_heat_mcell.png",gg_tmp)
