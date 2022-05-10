# renv setup
library(renv)

renv::init(bioconductor = TRUE)
renv::install("tidyverse")
renv::install("furrr")
renv::install("valr")
renv::install("vroom")
renv::install("parallel")
renv::install("data.tree")
renv::install("igraph")
renv::install("crayon")
renv::install("bedr")
renv::install("UpSetR")
renv::install("svglite")
renv::install("optparse")
renv::install("DiagrammeR")
renv::install("bioc::GenomicRanges")
renv::install("bioc::rtracklayer")


renv::snapshot()
