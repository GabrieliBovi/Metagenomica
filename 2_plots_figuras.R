library(tidyverse)
library(reshape2)
library(dplyr)
library(phyloseq)
library(ggsci)
library(ggpubr)

# Leitura de dados
ps <- readRDS("/Users/LS28_Seq/Documents/ref/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/ref/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")

# Seleciona amostras válidas
order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)

# Configurações
esconder_unclassified <- FALSE
categoria_others <- TRUE
abund_min <- 0.2

# Define o nível taxonômico desejado
taxrank <- "phylum"  # Alterar para "phylum", "class", "order", "family", "genus", ou "species"

# Aglomeração por nível taxonômico e cálculo de abundância
rank_seqtab <- tax_glom(seqtab, taxrank=taxrank)
rank_counts_tab <- otu_table(rank_seqtab)
rank_tax_vec <- as.vector(tax_table(rank_seqtab)[,taxrank])
colnames(rank_counts_tab) <- as.vector(rank_tax_vec)

rank_counts_tab2 <- otu_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))
rank_tax_vec2 <- as.vector(tax_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))[,taxrank])
colnames(rank_counts_tab2) <- as.vector(rank_tax_vec2)

# Cálculo de não classificados
unclassified_tax_counts <- rowSums(rank_counts_tab2) - rowSums(rank_counts_tab)
rank_and_unidentified_counts_tab <- cbind(rank_counts_tab, "Unclassified"=unclassified_tax_counts)

# Remover não classificados se necessário
if (esconder_unclassified){
  rank_and_unidentified_counts_tab <- rank_and_unidentified_counts_tab[
    ,!colnames(rank_and_unidentified_counts_tab) == "Unclassified"
  ]
}

# Filtrar por abundância mínima
asv_filt <- rank_and_unidentified_counts_tab[
  ,colSums(rank_and_unidentified_counts_tab) * 100 / sum(rank_and_unidentified_counts_tab) >= abund_min
]

# Adicionar categoria "Others" se necessário
if (categoria_others == TRUE & abund_min > 0){
  asv_filt <- cbind(
    asv_filt,
    data.frame(Others=rowSums(
      rank_and_unidentified_counts_tab[
        ,!(colnames(rank_and_unidentified_counts_tab) %in% colnames(asv_filt))
      ]
    )
    )
  )
}

# Calcular proporções
rank_taxa_proportions_tab <- apply(asv_filt, 1, function(x) x/sum(x)*100)
rank_taxa_proportions_tab_non0 <- rank_taxa_proportions_tab[
  rowSums(rank_taxa_proportions_tab) > 0,
]

# Transformar em formato "tidy" para ggplot
tidy_rank_taxa_proportions_tab <- reshape2::melt(rank_taxa_proportions_tab_non0,
                                                 value.name = "abundance",
                                                 varnames = c("taxa","sampleName")) %>%
  merge(data.frame(sample_data(seqtab)), by = "sampleName")

# Ajustar labels do gráfico
plot_labels <- sort(as.character(unique(tidy_rank_taxa_proportions_tab$taxa)))
if ("Others" %in% tidy_rank_taxa_proportions_tab$taxa){
  plot_labels <- c(plot_labels[plot_labels!="Others"],"Others")
}
if ("Unclassified" %in% plot_labels){
  plot_labels <- c(plot_labels[plot_labels!="Unclassified"],"Unclassified")
}

# Reorganizar fator de taxa
tidy_rank_taxa_proportions_tab$taxa <- factor(tidy_rank_taxa_proportions_tab$taxa, levels = plot_labels)
tidy_rank_taxa_proportions_tab <- tidy_rank_taxa_proportions_tab[rowSums(asv_filt) > 0, ]


## Criar o gráfico de barras phylum ------
plota_phylum <- ggplot(tidy_rank_taxa_proportions_tab) +
  geom_col(aes(x=sampleName2, y=abundance, fill=taxa), position = "stack") +
  theme_bw(base_size = 11) +
  labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(.8, "lines"),
        strip.background = element_rect(colour="white"),
        text=element_text(family="serif", color = "black", face = "plain")) +
  scale_fill_d3(palette = "category20c", alpha = 1)
print(plota_phylum)


## Criar o gráfico de heatmap para phylum ------
tidy_rank_taxa_proportions_tab_noUnclassOthers <- tidy_rank_taxa_proportions_tab[!tidy_rank_taxa_proportions_tab$taxa %in% c("Unclassified","Others"),]
tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa <- factor(tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa,levels=rev(plot_labels))

plota_phylum_h <- ggplot(tidy_rank_taxa_proportions_tab_noUnclassOthers) +
  geom_tile(aes(x = sampleName2, y = taxa, fill = abundance)) +
  scale_fill_gradient(high = "#7b4173f1", low = "gray95") +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 1.2,
    axis.title.x = element_blank(),           
    text = element_text(family = "serif", color = "black",  face = "plain"),
    strip.placement = "outside",
    panel.spacing = unit(.8, "lines"),
    legend.text = element_text(size = 11),    
    legend.title = element_text(size = 11),   
    axis.text.y = element_text(size = 11),    
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y = str_to_title(taxrank), fill = "Abundance (%)") +
  scale_y_discrete(
    position = "left",                         
    guide = guide_axis(position = "left")) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +  
  scale_y_discrete(expand = expansion(mult = c(0, 0)))

print(plota_phylum_h)
#ggsave('figuras/taxo_5/family_5.tiff', plot = plota_family_h, height = 1000, width = 1000, units = "px", dpi = 250)


# Class -------
# Leitura de dados
ps <- readRDS("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")

order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)
esconder_unclassified <- FALSE
categoria_others <- TRUE
abund_min <- 0.2

taxrank <- "class"  # Alterar para "phylum", "class", "order", "family", "genus", ou "species"

rank_seqtab <- tax_glom(seqtab, taxrank=taxrank)
rank_counts_tab <- otu_table(rank_seqtab)
rank_tax_vec <- as.vector(tax_table(rank_seqtab)[,taxrank])
colnames(rank_counts_tab) <- as.vector(rank_tax_vec)

rank_counts_tab2 <- otu_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))
rank_tax_vec2 <- as.vector(tax_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))[,taxrank])
colnames(rank_counts_tab2) <- as.vector(rank_tax_vec2)

unclassified_tax_counts <- rowSums(rank_counts_tab2) - rowSums(rank_counts_tab)
rank_and_unidentified_counts_tab <- cbind(rank_counts_tab, "Unclassified"=unclassified_tax_counts)

if (esconder_unclassified){
  rank_and_unidentified_counts_tab <- rank_and_unidentified_counts_tab[
    ,!colnames(rank_and_unidentified_counts_tab) == "Unclassified"
  ]
}

asv_filt <- rank_and_unidentified_counts_tab[
  ,colSums(rank_and_unidentified_counts_tab) * 100 / sum(rank_and_unidentified_counts_tab) >= abund_min
]

if (categoria_others == TRUE & abund_min > 0){
  asv_filt <- cbind(
    asv_filt,
    data.frame(Others=rowSums(
      rank_and_unidentified_counts_tab[
        ,!(colnames(rank_and_unidentified_counts_tab) %in% colnames(asv_filt))
      ]
    )
    )
  )
}

rank_taxa_proportions_tab <- apply(asv_filt, 1, function(x) x/sum(x)*100)
rank_taxa_proportions_tab_non0 <- rank_taxa_proportions_tab[
  rowSums(rank_taxa_proportions_tab) > 0,
]

tidy_rank_taxa_proportions_tab <- reshape2::melt(rank_taxa_proportions_tab_non0,
                                                 value.name = "abundance",
                                                 varnames = c("taxa","sampleName")) %>%
  merge(data.frame(sample_data(seqtab)), by = "sampleName")

plot_labels <- sort(as.character(unique(tidy_rank_taxa_proportions_tab$taxa)))
if ("Others" %in% tidy_rank_taxa_proportions_tab$taxa){
  plot_labels <- c(plot_labels[plot_labels!="Others"],"Others")
}
if ("Unclassified" %in% plot_labels){
  plot_labels <- c(plot_labels[plot_labels!="Unclassified"],"Unclassified")
}
tidy_rank_taxa_proportions_tab$taxa <- factor(tidy_rank_taxa_proportions_tab$taxa, levels = plot_labels)

## Criar o gráfico de barras 
plota_class <- ggplot(tidy_rank_taxa_proportions_tab) +
  geom_col(aes(x=sampleName, y=abundance, fill=taxa), position = "stack") +
  theme_bw(base_size = 11) +
  labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(.8, "lines"),
        strip.background = element_rect(colour="white"),
        text=element_text(family="serif", color = "black", face = "plain")) +
  scale_fill_d3(palette = "category20c", alpha = 1)
print(plota_class)


## Criar o gráfico de heatmap para class
tidy_rank_taxa_proportions_tab_noUnclassOthers <- tidy_rank_taxa_proportions_tab[!tidy_rank_taxa_proportions_tab$taxa %in% c("Unclassified","Others"),]
tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa <- factor(tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa,levels=rev(plot_labels))

plota_class_h <- ggplot(tidy_rank_taxa_proportions_tab_noUnclassOthers) +
  geom_tile(aes(x = sampleName2, y = taxa, fill = abundance)) +
  scale_fill_gradient(high = "#7b4173f1", low = "gray95") +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 1.2,
    axis.title.x = element_blank(),           
    text = element_text(family = "serif", color = "black",  face = "plain"),
    strip.placement = "outside",
    panel.spacing = unit(.8, "lines"),
    legend.text = element_text(size = 11),    
    legend.title = element_text(size = 11),   
    axis.text.y = element_text(size = 11),    
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y = str_to_title(taxrank), fill = "Abundance (%)") +
  scale_y_discrete(
    position = "left",                         
    guide = guide_axis(position = "left")) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +  
  scale_y_discrete(expand = expansion(mult = c(0, 0)))

print(plota_class_h)
#ggsave('figuras/taxo_5/class_5.tiff', plot = plota_class_h, height = 1000, width = 1000, units = "px", dpi = 250)

# Order -------
# Leitura de dados
ps <- readRDS("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")

order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)
esconder_unclassified <- FALSE
categoria_others <- TRUE
abund_min <- 0.2

taxrank <- "order"  # Alterar para "phylum", "class", "order", "family", "genus", ou "species"

rank_seqtab <- tax_glom(seqtab, taxrank=taxrank)
rank_counts_tab <- otu_table(rank_seqtab)
rank_tax_vec <- as.vector(tax_table(rank_seqtab)[,taxrank])
colnames(rank_counts_tab) <- as.vector(rank_tax_vec)

rank_counts_tab2 <- otu_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))
rank_tax_vec2 <- as.vector(tax_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))[,taxrank])
colnames(rank_counts_tab2) <- as.vector(rank_tax_vec2)

unclassified_tax_counts <- rowSums(rank_counts_tab2) - rowSums(rank_counts_tab)
rank_and_unidentified_counts_tab <- cbind(rank_counts_tab, "Unclassified"=unclassified_tax_counts)

if (esconder_unclassified){
  rank_and_unidentified_counts_tab <- rank_and_unidentified_counts_tab[
    ,!colnames(rank_and_unidentified_counts_tab) == "Unclassified"
  ]
}

asv_filt <- rank_and_unidentified_counts_tab[
  ,colSums(rank_and_unidentified_counts_tab) * 100 / sum(rank_and_unidentified_counts_tab) >= abund_min
]

if (categoria_others == TRUE & abund_min > 0){
  asv_filt <- cbind(
    asv_filt,
    data.frame(Others=rowSums(
      rank_and_unidentified_counts_tab[
        ,!(colnames(rank_and_unidentified_counts_tab) %in% colnames(asv_filt))
      ]
    )
    )
  )
}

rank_taxa_proportions_tab <- apply(asv_filt, 1, function(x) x/sum(x)*100)
rank_taxa_proportions_tab_non0 <- rank_taxa_proportions_tab[
  rowSums(rank_taxa_proportions_tab) > 0,
]

tidy_rank_taxa_proportions_tab <- reshape2::melt(rank_taxa_proportions_tab_non0,
                                                 value.name = "abundance",
                                                 varnames = c("taxa","sampleName")) %>%
  merge(data.frame(sample_data(seqtab)), by = "sampleName")

plot_labels <- sort(as.character(unique(tidy_rank_taxa_proportions_tab$taxa)))
if ("Others" %in% tidy_rank_taxa_proportions_tab$taxa){
  plot_labels <- c(plot_labels[plot_labels!="Others"],"Others")
}
if ("Unclassified" %in% plot_labels){
  plot_labels <- c(plot_labels[plot_labels!="Unclassified"],"Unclassified")
}
tidy_rank_taxa_proportions_tab$taxa <- factor(tidy_rank_taxa_proportions_tab$taxa, levels = plot_labels)

## Criar o gráfico de barras 
plota_order <- ggplot(tidy_rank_taxa_proportions_tab) +
  geom_col(aes(x=sampleName2, y=abundance, fill=taxa), position = "stack") +
  theme_bw(base_size = 11) +
  labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(.8, "lines"),
        strip.background = element_rect(colour="white"),
        text=element_text(family="serif", color = "black", face = "plain")) +
  scale_fill_d3(palette = "category20c", alpha = 1)

## Criar o gráfico de heatmap para class
tidy_rank_taxa_proportions_tab_noUnclassOthers <- tidy_rank_taxa_proportions_tab[!tidy_rank_taxa_proportions_tab$taxa %in% c("Unclassified","Others"),]
tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa <- factor(tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa,levels=rev(plot_labels))

plota_order_h <- ggplot(tidy_rank_taxa_proportions_tab_noUnclassOthers) +
  geom_tile(aes(x = sampleName2, y = taxa, fill = abundance)) +
  scale_fill_gradient(high = "#7b4173f1", low = "gray95") +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 1.2,
    axis.title.x = element_blank(),           
    text = element_text(family = "serif", color = "black",  face = "plain"),
    strip.placement = "outside",
    panel.spacing = unit(.8, "lines"),
    legend.text = element_text(size = 11),    
    legend.title = element_text(size = 11),   
    axis.text.y = element_text(size = 11),    
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y = str_to_title(taxrank), fill = "Abundance (%)") +
  scale_y_discrete(
    position = "left",                         
    guide = guide_axis(position = "left")) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +  
  scale_y_discrete(expand = expansion(mult = c(0, 0)))

print(plota_order_h)
#ggsave('figuras/taxo_5/order_5.tiff', plot = plota_order_h, height = 1000, width = 1000, units = "px", dpi = 250)


# Family -------
# Leitura de dados
ps <- readRDS("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")

order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)
esconder_unclassified <- FALSE
categoria_others <- TRUE
abund_min <- 0.2

taxrank <- "family"  # Alterar para "phylum", "class", "order", "family", "genus", ou "species"

rank_seqtab <- tax_glom(seqtab, taxrank=taxrank)
rank_counts_tab <- otu_table(rank_seqtab)
rank_tax_vec <- as.vector(tax_table(rank_seqtab)[,taxrank])
colnames(rank_counts_tab) <- as.vector(rank_tax_vec)

rank_counts_tab2 <- otu_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))
rank_tax_vec2 <- as.vector(tax_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))[,taxrank])
colnames(rank_counts_tab2) <- as.vector(rank_tax_vec2)

unclassified_tax_counts <- rowSums(rank_counts_tab2) - rowSums(rank_counts_tab)
rank_and_unidentified_counts_tab <- cbind(rank_counts_tab, "Unclassified"=unclassified_tax_counts)

if (esconder_unclassified){
  rank_and_unidentified_counts_tab <- rank_and_unidentified_counts_tab[
    ,!colnames(rank_and_unidentified_counts_tab) == "Unclassified"
  ]
}

asv_filt <- rank_and_unidentified_counts_tab[
  ,colSums(rank_and_unidentified_counts_tab) * 100 / sum(rank_and_unidentified_counts_tab) >= abund_min
]

if (categoria_others == TRUE & abund_min > 0){
  asv_filt <- cbind(
    asv_filt,
    data.frame(Others=rowSums(
      rank_and_unidentified_counts_tab[
        ,!(colnames(rank_and_unidentified_counts_tab) %in% colnames(asv_filt))
      ]
    )
    )
  )
}

rank_taxa_proportions_tab <- apply(asv_filt, 1, function(x) x/sum(x)*100)
rank_taxa_proportions_tab_non0 <- rank_taxa_proportions_tab[
  rowSums(rank_taxa_proportions_tab) > 0,
]

tidy_rank_taxa_proportions_tab <- reshape2::melt(rank_taxa_proportions_tab_non0,
                                                 value.name = "abundance",
                                                 varnames = c("taxa","sampleName")) %>%
  merge(data.frame(sample_data(seqtab)), by = "sampleName")

plot_labels <- sort(as.character(unique(tidy_rank_taxa_proportions_tab$taxa)))
if ("Others" %in% tidy_rank_taxa_proportions_tab$taxa){
  plot_labels <- c(plot_labels[plot_labels!="Others"],"Others")
}
if ("Unclassified" %in% plot_labels){
  plot_labels <- c(plot_labels[plot_labels!="Unclassified"],"Unclassified")
}
tidy_rank_taxa_proportions_tab$taxa <- factor(tidy_rank_taxa_proportions_tab$taxa, levels = plot_labels)

## Criar o gráfico de barras 
plota_family <- ggplot(tidy_rank_taxa_proportions_tab) +
  geom_col(aes(x=sampleName2, y=abundance, fill=taxa), position = "stack") +
  theme_bw(base_size = 11) +
  labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(.8, "lines"),
        strip.background = element_rect(colour="white"),
        text=element_text(family="serif", color = "black", face = "plain")) +
  scale_fill_d3(palette = "category20c", alpha = 1)

## Criar o gráfico de heatmap
tidy_rank_taxa_proportions_tab_noUnclassOthers <- tidy_rank_taxa_proportions_tab[!tidy_rank_taxa_proportions_tab$taxa %in% c("Unclassified","Others"),]
tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa <- factor(tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa,levels=rev(plot_labels))

plota_family_h <- ggplot(tidy_rank_taxa_proportions_tab_noUnclassOthers) +
  geom_tile(aes(x = sampleName2, y = taxa, fill = abundance)) +
  scale_fill_gradient(high = "#7b4173f1", low = "gray95") +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 1.2,
    axis.title.x = element_blank(),           
    text = element_text(family = "serif", color = "black",  face = "plain"),
    strip.placement = "outside",
    panel.spacing = unit(.8, "lines"),
    legend.text = element_text(size = 11),    
    legend.title = element_text(size = 11),   
    axis.text.y = element_text(size = 11),    
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y = str_to_title(taxrank), fill = "Abundance (%)") +
  scale_y_discrete(
    position = "left",                         
    guide = guide_axis(position = "left")) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +  
  scale_y_discrete(expand = expansion(mult = c(0, 0)))

print(plota_family_h)
#ggsave('figuras/taxo_5/family_5.tiff', plot = plota_family_h, height = 1000, width = 1000, units = "px", dpi = 250)


# Genus -------
# Leitura de dados
ps <- readRDS("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")

order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)
esconder_unclassified <- FALSE
categoria_others <- TRUE
abund_min <- 0.2

taxrank <- "genus"  # Alterar para "phylum", "class", "order", "family", "genus", ou "species"

rank_seqtab <- tax_glom(seqtab, taxrank=taxrank)
rank_counts_tab <- otu_table(rank_seqtab)
rank_tax_vec <- as.vector(tax_table(rank_seqtab)[,taxrank])
colnames(rank_counts_tab) <- as.vector(rank_tax_vec)

rank_counts_tab2 <- otu_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))
rank_tax_vec2 <- as.vector(tax_table(tax_glom(seqtab, taxrank=taxrank, NArm=FALSE))[,taxrank])
colnames(rank_counts_tab2) <- as.vector(rank_tax_vec2)

unclassified_tax_counts <- rowSums(rank_counts_tab2) - rowSums(rank_counts_tab)
rank_and_unidentified_counts_tab <- cbind(rank_counts_tab, "Unclassified"=unclassified_tax_counts)

if (esconder_unclassified){
  rank_and_unidentified_counts_tab <- rank_and_unidentified_counts_tab[
    ,!colnames(rank_and_unidentified_counts_tab) == "Unclassified"
  ]
}

asv_filt <- rank_and_unidentified_counts_tab[
  ,colSums(rank_and_unidentified_counts_tab) * 100 / sum(rank_and_unidentified_counts_tab) >= abund_min
]

if (categoria_others == TRUE & abund_min > 0){
  asv_filt <- cbind(
    asv_filt,
    data.frame(Others=rowSums(
      rank_and_unidentified_counts_tab[
        ,!(colnames(rank_and_unidentified_counts_tab) %in% colnames(asv_filt))
      ]
    )
    )
  )
}

rank_taxa_proportions_tab <- apply(asv_filt, 1, function(x) x/sum(x)*100)
rank_taxa_proportions_tab_non0 <- rank_taxa_proportions_tab[
  rowSums(rank_taxa_proportions_tab) > 0,
]

tidy_rank_taxa_proportions_tab <- reshape2::melt(rank_taxa_proportions_tab_non0,
                                                 value.name = "abundance",
                                                 varnames = c("taxa","sampleName")) %>%
  merge(data.frame(sample_data(seqtab)), by = "sampleName")

plot_labels <- sort(as.character(unique(tidy_rank_taxa_proportions_tab$taxa)))
if ("Others" %in% tidy_rank_taxa_proportions_tab$taxa){
  plot_labels <- c(plot_labels[plot_labels!="Others"],"Others")
}
if ("Unclassified" %in% plot_labels){
  plot_labels <- c(plot_labels[plot_labels!="Unclassified"],"Unclassified")
}
tidy_rank_taxa_proportions_tab$taxa <- factor(tidy_rank_taxa_proportions_tab$taxa, levels = plot_labels)

## Criar o gráfico de barras 
plota_genus <- ggplot(tidy_rank_taxa_proportions_tab) +
  geom_col(aes(x=sampleName2, y=abundance, fill=taxa), position = "stack") +
  theme_bw(base_size = 11) +
  labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(.8, "lines"),
        strip.background = element_rect(colour="white"),
        text=element_text(family="serif", color = "black", face = "plain")) +
  scale_fill_d3(palette = "category20c", alpha = 1)

## Criar o gráfico de heatmap
tidy_rank_taxa_proportions_tab_noUnclassOthers <- tidy_rank_taxa_proportions_tab[!tidy_rank_taxa_proportions_tab$taxa %in% c("Unclassified","Others"),]
tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa <- factor(tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa,levels=rev(plot_labels))

plota_genus_h <- ggplot(tidy_rank_taxa_proportions_tab_noUnclassOthers) +
  geom_tile(aes(x = sampleName2, y = taxa, fill = abundance)) +
  scale_fill_gradient(high = "#7b4173f1", low = "gray95") +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 1.2,
    axis.title.x = element_blank(),           
    text = element_text(family = "serif", color = "black", face = "plain"),
    strip.placement = "outside",
    panel.spacing = unit(.8, "lines"),
    legend.text = element_text(size = 11),    
    legend.title = element_text(size = 11),   
    axis.text.y = element_text(size = 11),    
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y = str_to_title(taxrank), fill = "Abundance (%)") +
  scale_y_discrete(
    position = "left",                         
    guide = guide_axis(position = "left")) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +  
  scale_y_discrete(expand = expansion(mult = c(0, 0)))

print(plota_genus_h)
#ggsave('figuras/taxo_5/genus_5.tiff', plot = plota_genus_h, height = 1000, width = 1000, units = "px", dpi = 250)


# Layout: Plotar todos os graficos de barra ------
# Organizar os gráficos em um layout de 2 colunas e 3 linhas
final_plot <- ggarrange(
  ggarrange(plota_phylum, plota_class, labels = c("A", "B"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_order, plota_family, labels = c("C", "D"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_genus, labels = c("E", "F"), ncol = 2, widths = c(1, 1)),
  nrow = 3, align = "hv") +bgcolor("white")
print(final_plot)
#ggsave('figuras/taxo_5_2/layout_barra_2.tiff', plot = final_plot, height = 2999, width = 2800, units = "px", dpi = 250)


# Layout: Plotar todos heatmaps -------
final_plot_h <- ggarrange(
  ggarrange(plota_phylum_h, plota_class_h, labels = c("A", "B"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_order_h, plota_family_h, labels = c("C", "D"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_genus_h, labels = c("E", "F"), ncol = 2, widths = c(1, 1)),
  nrow = 3, align = "hv") + bgcolor("white")
print(final_plot_h)
#ggsave('figuras/taxo_5_2/final_plot_h_2.tiff', plot = final_plot_h,height = 2800, width = 3500, units = "px", dpi = 290)



