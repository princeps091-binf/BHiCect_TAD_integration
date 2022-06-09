library(tidyverse)
library(vroom)
library(Matrix)
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

hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.na(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

full_f_mat<-function(cl_mat,res,var){
  
  range_5kb<-range(as.numeric(unique(c(cl_mat$bin.A,cl_mat$bin.B))))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$bin.A)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$bin.B)]
  
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat[[var]]))
  return(chr_mat)
}

#---------------------------------------------
hmec_res_folder<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
hmec_dat_folder<-"~/Documents/multires_bhicect/data/HMEC/"
tad_inter_hmec_tbl<-get_tbl_in_fn("./data/intersection_out_tbl/HMEC_top_TAD_cl_intersect.Rda")


tad_inter_cl_hmec_tbl<-cl_Grange_build_fn(tad_inter_hmec_tbl,hmec_res_folder)

tad_inter_cl_hmec_tbl %>% 
  mutate(res=str_split_fixed(max.cl,"_",2)[,1]) %>% 
  group_by(chr,res) %>% 
  summarise(n=n()) %>% 
  filter(chr=="chr20")

#---------------------------------------------
tmp_res<-"10kb"
chromo<-"chr20"

tad_chr_tbl<-tad_inter_cl_hmec_tbl %>% 
  filter(chr==chromo)

chr_dat<-hic_dat_in(hmec_dat_folder,tmp_res,chromo)

chr_bin_tbl<-tibble(chr=chromo,res=tmp_res,start=unique(c(chr_dat$X1,chr_dat$X2))) %>% 
  mutate(ends=start + res_num[res] - 1)
  
bin_Grange<-GRanges(seqnames=chr_bin_tbl$chr,
                    ranges = IRanges(start=chr_bin_tbl$start,
                                     end=chr_bin_tbl$ends
                    ))
tad_chr_tbl<-tad_chr_tbl %>% 
  mutate(bins=lapply(tad_chr_tbl$GRange,function(x){
    chr_bin_tbl$start[unique(subjectHits(findOverlaps(x,bin_Grange)))]
  })
  )

tad_inter_mat<-do.call(bind_rows,lapply(seq_along(tad_chr_tbl$bins),function(i){
  tmp_bin<-tad_chr_tbl$bins[[i]]
  return(expand_grid(bin.A=tmp_bin,bin.B=tmp_bin) %>% 
    mutate(value=i))
}))
tad_mat_tbl<-tad_mat_tbl %>% 
  left_join(.,tad_inter_mat)

tad_mat<-full_f_mat(tad_inter_mat,res_num[tmp_res],"value")

image(log10(as.matrix(tad_mat)))

png(paste0('~/Documents/multires_bhicect/Poster/img/F1/HMEC_',chromo,"_TAD","_mat",'.png'), width =40,height = 40,units = 'mm',type='cairo',res=5000)
par(mar = c(0, 0, 0,0))
plot.new()
image(log10(as.matrix(tad_mat)))
dev.off()

tad_chr_tbl<-tad_chr_tbl %>% 
  mutate(cl.bins=lapply(tad_chr_tbl$max.GRange,function(x){
    chr_bin_tbl$start[unique(subjectHits(findOverlaps(x,bin_Grange)))]
  })
  )

cl_inter_mat<-do.call(bind_rows,lapply(seq_along(tad_chr_tbl$cl.bins),function(i){
  tmp_bin<-tad_chr_tbl$cl.bins[[i]]
  return(expand_grid(bin.A=tmp_bin,bin.B=tmp_bin) %>% 
           mutate(value=i))
}))

tad_mat<-full_f_mat(cl_inter_mat,res_num[tmp_res],"value")

image(log10(as.matrix(tad_mat)))

png(paste0('~/Documents/multires_bhicect/Poster/img/F1/HMEC_',chromo,"_TAD_cl","_mat",'.png'), width =40,height = 40,units = 'mm',type='cairo',res=5000)
par(mar = c(0, 0, 0,0))
plot.new()
image(log10(as.matrix(tad_mat)))
dev.off()
