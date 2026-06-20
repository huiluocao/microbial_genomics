#https://github.com/gtonkinhill/rhierbaps
#https://github.com/ocbe-uio/rBAPS
#https://github.com/gtonkinhill/fastbaps

install.packages("rhierbaps")
library(rhierbaps)
library(ggtree)
library(ape)
fasta.file.name <- "data/phylogeny/hku_pm109_parsnp.snp"
snp.matrix <- load_fasta(fasta.file.name)
hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 20)
hb.results$partition.df
