#### reading trees
mt147_tre <- read.tree('./data/01_all_bins/all_bins/mt_vhq_drep95_50c_bins_concatenated.aln.contree') %>%
  phytools::midpoint.root()

mt147_p <- ggtree(mt147_tre, layout = "rectangular", open.angle=15, size=.2) + ## fan
  geom_nodepoint(aes(size = as.numeric(label)), stroke = 0) +
  xlim(-.1, 4)+
  geom_tiplab(size=2, align = TRUE, family = 'Times')+
  scale_size_continuous(name = 'Bootstrap', range = c(0, 1.5), 
                        guide = guide_legend(order = 1)) +
  theme(text = element_text(family = 'Times'))
mt147_p

ar_tre <- read.tree('./data/04_phylogenomics/mt_vhq_drep95_50c/mt_vhq_drep95_50c_gtdbtk214.ar53.user_msa.fasta.contree') %>%
  phytools::midpoint.root()

ar_p <- ggtree(ar_tre, layout = "rectangular", open.angle=15, size=.5) + ## fan
  geom_nodepoint(aes(size = as.numeric(label)), stroke = 0) +
  geom_tiplab()+
  xlim(-.1, 4)+
  scale_size_continuous(name = 'Bootstrap', range = c(0, 1.5), 
                        guide = guide_legend(order = 1)) +
  theme(text = element_text(family = 'Times'))
ar_p

bac_tre <- read.tree('./data/04_phylogenomics/mt_vhq_drep95_50c/mt_vhq_drep95_50c_gtdbtk214.bac120.user_msa.fasta.contree') %>%
  phytools::midpoint.root()

bac_p <- ggtree(bac_tre, layout = "rectangular", open.angle=15, size=.2) + ## fan
  geom_nodepoint(aes(size = as.numeric(label)), stroke = 0) +
  scale_size_continuous(name = 'Bootstrap', range = c(0, 1.5), 
                        guide = guide_legend(order = 1)) +
  #geom_tiplab(size = 3)+
  xlim(-.1, 4)+
  theme(text = element_text(family = 'Times'))
bac_p

d = fortify(ar_tre)
d = fortify(bac_tre)
d = fortify(mt147_tre)
d = subset(d, isTip)
d
ar_bins=with(d, label[order(y, decreasing=T)])
bac_bins=with(d, label[order(y, decreasing=T)])
mt_147bins=with(d, label[order(y, decreasing=T)])
data.frame(c(ar_bins,bac_bins))

#### add phylum info
dat1 <- mt_bins147_new[c('bin_id', 'new_phylum')]
#cluster_all <- c(paste0('Cluster', 1:50), 'OG') 
#dat1$Clusters <- factor(dat1$`new cluster`,
#                        levels = intersect(cluster_all, dat1$`new cluster`))
colnames(dat1) <- c('ID','phylumn')
head(dat1)
dim(dat1)
cols <- c('#807f00', '#f97e00', '#70ffcb', '#6695a0', '#71fe5e', '#fdcb61', 
          '#f86565', '#abcf95', '#6bceff', '#a395ce', '#f78dff', '#f7f7a1', '#d2e6c8',
          '#f2a9d0', '#eebfc1', '#d0d6de', '#c7e3f7', '#dff99c', '#f1d5ae', '#b6b7f5',
          '#7e3f00', '#656aff')

library(randomcoloR)
cols <-distinctColorPalette(43)
p1 <- bac_p +
  new_scale_fill() +
  geom_fruit(
    data = dat1,
    geom=geom_tile, 
    mapping=aes(y=ID, fill=phylumn),
    offset = .65, #.15,
    width = .19,
    alpha = .5) + 
  scale_fill_manual(name = 'phylumn', values = cols, guide=guide_legend(order=2)) +
  geom_tiplab(size=2, align = TRUE, linetype = '8f', color = 'grey40', family='Times',
              linesize = .05, offset = .01)+
  #coord_flip()+
  theme(legend.position = 'left')
p1

### add GC
dat2 <- mt_bins147_new[c('bin_id', 'GC')]
colnames(dat2) <- c('ID','GC(%)')
#dat4 <- df[c('ID', 'GC (%)')]
dat2$GC <- as.numeric(dat2$`GC(%)`) #- 50

p2 <- p1 +
  new_scale_fill() +
  geom_fruit(
    data = dat2,
    geom=geom_tile,
    mapping=aes(y=ID, fill = GC),
    color = 'white',
    width = .05,
    offset = .12
  ) +
  scale_fill_gradient(name = 'GC content (%)', low = 'white', high = '#4d9220',
                      guide=guide_colorbar(order = 5))
p2
#### add number of genomes of each cluster
dat3 <- mt_bins147_new[c('bin_id', 'total')]
colnames(dat3) <- c('ID','Number')
#dat4 <- df[c('ID', 'GC (%)')]
dat3$Number <- as.numeric(dat3$Number) #- 50

p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = dat3,
    geom=geom_tile,
    mapping=aes(y=ID, fill = Number),
    color = 'white',
    width = .05,
    offset = .12
  ) +
  scale_fill_gradient(name = 'Number of genomes', low = 'white', high = '#7e3f00',
                      guide=guide_colorbar(order = 6), )

## layer4: heatmap
### BGC count for GCC
dat4 <- table(df$GCC) %>% data.frame %>% { colnames(.) <- c('GCC', 'BGC_count'); .}
head(hq95cluster_cov)
dim(hq95cluster_cov)
dat3 <- mt_bins147_new[c('bin_id', 'total')]

p4 <- p3 +
  new_scale_fill() +
  geom_fruit(
    data = dat4,
    geom=geom_tile, 
    mapping=aes(y=GCC, fill=log2(BGC_count+1)),
    color = 'white',
    offset=.15,
    width=4
  ) + 
  scale_fill_distiller(palette = 'Greys', direction = 1,
                       breaks = c(log2(1+1), log2(100+1), log2(1000+1)),
                       labels = c(1, 100, 1000)) +
  guides(fill = guide_colourbar(title = 'BGC count', direction = "horizontal"))


ar_p1 <- ar_p +
  new_scale_fill() +
  geom_fruit(
    data = dat1,
    geom=geom_tile, 
    mapping=aes(y=ID, fill=phylumn),
    offset = 1.5, #.15,
    width = .19,
    alpha = .5) + 
  scale_fill_manual(name = 'phylumn', values = cols, guide=guide_legend(order=2)) +
  #geom_tiplab(size=2, align = TRUE, linetype = '8f', color = 'grey40', family='Times',
  #            linesize = .05, offset = .01)+
  #coord_flip()+
  theme(legend.position = 'left')
ar_p1


mt147_p1 <- mt147_p +
  new_scale_fill() +
  geom_fruit(
    data = dat1,
    geom=geom_tile, 
    mapping=aes(y=ID, fill=phylumn),
    offset = .65, #.15,
    width = .19,
    alpha = .8
    ) + 
  #scale_fill_manual(name = 'phylumn', values = cols, guide=guide_legend(order=2)) +
  #geom_tiplab(size=1, align = TRUE, linetype = '8f', color = 'grey40', family='Times',
  #            linesize = .05, offset = .01)+
  #coord_flip()+
  scale_fill_manual(name = 'phylumn', values = c("p__Actinobacteriota"="#a6cee3",
                                                 "p__Bacteroidota" ="#1f78b4",
                                                 "p__Chloroflexota"="#b2df8a",
                                                 "p__Elusimicrobiota" ="#33a02c",
                                                 "p__Gemmatimonadota" ="#fb9a99",
                                                 "p__Marinisomatota" ="#e31a1c",
                                                 "p__Nitrospirota"= "#fdbf6f",
                                                 "p__Planctomycetota" = "#ff7f00",
                                                 "p__Proteobacteria" ="#cab2d6",
                                                 "p__Thermoproteota" = "#6a3d9a",
                                                 "others"="#bdbdbd"
                                                 ), guide=guide_legend(order=2)) +
  
  theme(legend.position = 'left')
mt147_p1

"#a6cee3"
"#1f78b4"
"#b2df8a"
"#33a02c"
"#fb9a99"
"#e31a1c"
"#fdbf6f"
"#ff7f00"
"#cab2d6"
"#6a3d9a"
"#ffff99"
"#b15928"


### add GC
dat2 <- mt_bins147_new[c('bin_id', 'GC')]
colnames(dat2) <- c('ID','GC(%)')
#dat4 <- df[c('ID', 'GC (%)')]
dat2$GC <- as.numeric(dat2$`GC(%)`) #- 50

mt147_p2 <- mt147_p1 +
  new_scale_fill() +
  geom_fruit(
    data = dat2,
    geom=geom_tile,
    mapping=aes(y=ID, fill = GC),
    color = 'white',
    width = .05,
    offset = .22
  ) +
  scale_fill_gradient(name = 'GC content (%)', low = 'white', high = '#4d9220',
                      guide=guide_colorbar(order = 5))
mt147_p2
#### add number of genomes of each cluster
dat3 <- mt_bins147_new[c('bin_id', 'total')]
colnames(dat3) <- c('ID','Number')
#dat4 <- df[c('ID', 'GC (%)')]
dat3$Number <- as.numeric(dat3$Number) #- 50

p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = dat3,
    geom=geom_tile,
    mapping=aes(y=ID, fill = Number),
    color = 'white',
    width = .05,
    offset = .12
  ) +
  scale_fill_gradient(name = 'Number of genomes', low = 'white', high = '#7e3f00',
                      guide=guide_colorbar(order = 6), )
