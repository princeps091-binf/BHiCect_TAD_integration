# optparse test
library(optparse)

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
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

TAD_file<-opt$tad
BHiCect_res_folder<-opt$cluster
nworker<-as.numeric(opt$workers)
out_file<-opt$output

message("input: ", class(nworker))