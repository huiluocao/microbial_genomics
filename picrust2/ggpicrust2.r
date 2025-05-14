###### picrust
pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools",
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq",
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

library(ggpicrust2)
library(tidyverse)

devtools::install_github("cafferychen777/MicrobiomeStat")
library(MicrobiomeStat)

metadata1<-column_to_rownames(metadata,var='Sample')
metadata1<-rownames_to_column(metadata1,var='sample_name')

ko_abundance <- read.delim("KO_metagenome_out/pred_metagenome_unstrat.tsv")
ko_abundance_filtered <- ko_abundance[rowSums(ko_abundance > 0) >= 3, ]
metadata_G <- metadata1[metadata1$type == "G", ]
ko_abundance_G<-ko_abundance_filtered%>%
  subset(select=c('X.NAME',metadata_G$sample_name))
results_file_G <- ggpicrust2(data = ko_abundance_G,
                                 metadata = metadata_G,
                                 group = "Group",
                                 pathway = "KO",
                                 daa_method = "ALDEx2",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")
