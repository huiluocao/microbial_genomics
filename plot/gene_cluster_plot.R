
library(ggplot2)
library(gggenes)
library(magrittr)
library(xlsx)
extrafont::loadfonts()
rm(list = ls())
options(stringsAsFactors = F)

setwd('~/Daily/jindundun/20220710/gene-cluster_plot')
sam <- c('M23_bin.5', 'M2_bin.4', 'M3_bin.14', 'M9_bin.5', 'M22_bin.5')

# BGC0001017 bed
bed <- read.delim('input/BGC0001017.bed', header = F)
colnames(bed) <- c('start', 'end', 'strand', 'gene')
bed$genome <- 'BGC0001017'

# gene info
gene_symbol <- read.xlsx('input/microsyste_bgc_gene.xlsx',
                         sheetName = 'bgc_gene')

df <- lapply(sam, function(x) {
  
  genes <- read.delim(
    file = paste0('input/microcystis_bgc_blastp/', x, '_micro_bgc_blastp.txt'),
    header = F) %>%
    subset(V3 > 70 & V13 > 70) %>%
    .[grep('BGC0001017', .$V2), 'V2'] %>%
    plyr::mapvalues(from = gene_symbol$bgc_id, to = gene_symbol$gene)
 
  target <- bed
  target[!(target$gene %in% genes), 'gene'] <- 'Not Found'
  target$genome <- x
  return(target)

})  %>% do.call(rbind, .)
df <- rbind(bed, df)


df$direction <- ifelse(df$strand == '+', 1, -1)
df$genome <- factor(df$genome, levels = c('BGC0001017', sam))


## colors
cols <- setNames(c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99',
                   '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A',
                   '#FFFF99', '#B15928', '#6B9DF8', '#984EA3', '#7B967F',
                   '#C85080', '#666666', 'white'),
                   c(bed$gene, 'Not Found')
)


pdf('gene-cluster_plot.pdf', family="Times New Roman", 9, 5.3)
ggplot(df, aes(xmin = start, xmax = end, y = genome, fill = gene,
               forward = direction)) +
  geom_gene_arrow(arrow_body_height = grid::unit(5, "mm"),
                  arrowhead_width = grid::unit(6, "mm"),
                  arrowhead_height = grid::unit(7, "mm")) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  scale_fill_manual(values = cols) +
  theme_genes() +
  theme(axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_blank())
dev.off()


df_filtered <- subset(df, gene != 'Not Found') %>% .[c(5,1,2,4,3)]

file.copy('input/microsyste_bgc_gene_20220707_1b.xlsx', './microsyste_bgc_gene_20220707_1b.xlsx',
          overwrite = T)
write.xlsx(df_filtered, 'microsyste_bgc_gene_20220707_1b.xlsx', sheetName = 'Microcystin_20220718', 
           row.names = F, append = T)

