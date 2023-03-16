library(tidyverse)
library(TopDom)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
topDom_folder<-"./data/TopDom_data/HMEC/"
topDom_res_file<-"./data/TopDom_data/HMEC/TopDom_res/topdom_tads.Rda"
chr_set<-grep("^chr",str_split_fixed(list.files(topDom_folder),"_",2)[,1],value=T)
topDom_bed_tbl<-do.call(bind_rows,lapply(chr_set,function(chromo){
  message(chromo)
  topdom_file<-paste0("./data/TopDom_data/HMEC/",chromo,"_topdom_data.tsv")
  fit <- TopDom(topdom_file, window.size = 5L,debug=T)
  return(fit$bed %>% as_tibble)
}))

save(topDom_bed_tbl,file=topDom_res_file)
