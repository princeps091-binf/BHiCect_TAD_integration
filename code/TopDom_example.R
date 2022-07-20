library(tidyverse)
library(TopDom)
path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
chr <- "chr19"
pathname <- file.path(path, sprintf("nij.%s.gz", chr))
data <- readHiC(pathname, chr = chr, binSize = 40e3)
print(data)
str(data)

fit <- TopDom(tmp_data, window.size = 5L,debug=T)

## Subset TopDomData object
td <- subset(subset(fit$domain, tag == "domain"), size == max(size))
td<-fit$domain %>% 
  slice(8)
data_s <- subsetByRegion(data, region = td, margin = 0.9999)
vp <- grid::viewport(angle = -45, width = 0.7, y = 0.3)
gg <- ggCountHeatmap(data_s)
gg <- gg + ggDomain(td, color = "#cccc00") + ggDomainLabel(td)
print(gg, newpage = TRUE, vp = vp)

gg <- ggCountHeatmap(data_s, colors = list(mid = "white", high = "black"))
gg_td <- ggDomain(td, delta = 0.08)
dx <- attr(gg_td, "gg_params")$dx
gg <- gg + gg_td + ggDomainLabel(td, vjust = 2.5)
print(gg, newpage = TRUE, vp = vp)
