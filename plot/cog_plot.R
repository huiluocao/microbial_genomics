
library(ggplot2)
library(randomcoloR)
library(dplyr)
library(ggsci)
library(scales)
show_col(pal_npg('nrc')(10)[c(1,4)])

cog4raw <- read.table("./cog4genomes.txt", 
                        head = TRUE, sep = "\t",check.names = F)
dim(cog4raw)
head(cog4raw)
s<-lapply(cog4raw$class, function (x) unique(unlist(strsplit(as.character(x), ""))))
s
length(s)
cog4raw_rep<-data.frame(cog = rep(cog4raw$cog, sapply(s, length)), 
                        GCA_003149185.1 = rep(cog4raw$GCA_003149185.1, sapply(s, length)),
                        GCA_012844305.1 = rep(cog4raw$GCA_012844305.1, sapply(s, length)),
                        GCA_000379765.1 = rep(cog4raw$GCA_000379765.1, sapply(s, length)),
                        WSW3_B12= rep(cog4raw$`WSW3-B12`, sapply(s, length)),
                        class = unlist(s))

dim(cog4raw_rep)
head(cog4raw_rep)
cog4raw_rep[,-c(1)] %>% 
  group_by(class) %>% 
  summarise_each(list(sum))

cog4_class <-cog4raw_rep[,-c(1)]%>%
  group_by(class) %>% 
  summarise_each(list(sum))%>%
  pivot_longer(!class, names_to = "species", values_to = "count")%>%
  data.frame()
head(cog4_class)

cog_cat<-read.table('../../LMK/cog_fun-20.tab',head=F, sep='\t', 
                    stringsAsFactors = FALSE,check.names = F)
head(cog_cat)

cog4_category<-merge(cog4_class,cog_cat,by.x='class',by.y = 'V1', all.x=TRUE)
head(cog4_category)

ggplot(data = cog4_category, aes(V3, count, fill = species)) + 
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  scale_fill_manual(name='Species',values=pal_npg(c("nrc"))(10)[c(1,2,3,4)],
                    labels = c("Flexithrix dorotheae DSM 6795", 
                               "Sediminitomix flava DSM 28229",
                               "Flammeovirga aprica JL-4",
                               "WSW3-B12"))+
  labs(y="Number of Genes", x='Categories of COG')+
  #scale_fill_manual(name="Categories of COG",values =distinctColorPalette(25))+
  theme(axis.text.y = element_text(size = 12,family = "Times"),
        axis.text.x = element_text(size = 12,family = "Times"),
        axis.title.x = element_text(size = 16,family = "Times"),
        axis.title.y = element_text(size = 16,family = "Times"),
        legend.position = 'bottom',
  legend.text=element_text(size=14,family = "Times"),
  legend.title=element_text(size=16,family = "Times"))+
  guides(fill=guide_legend(ncol=2))+
  coord_flip()
