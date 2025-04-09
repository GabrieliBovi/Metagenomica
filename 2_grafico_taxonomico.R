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
abund_min <- 5 # 5%

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

# Definir a ordem das fases - Por eemplo:
tidy_rank_taxa_proportions_tab$Fase <- factor(
  tidy_rank_taxa_proportions_tab$Fase, levels = c("Initial", "Final"))

tidy_rank_taxa_proportions_tab$sampleName2  <- factor(
  tidy_rank_taxa_proportions_tab$sampleName2, 
  levels = c("End. CT","NC CT","Blank", "EXP-BC-10", "EXP-BC-20","EXP-BC-30","EXP-CAG-10","EXP-CAG-20","EXP-CAG-30"))

## Criar o gráfico de barras phylum ------
barra_phylum <-  # barra_species  barra_genus  barra_family  barra_order
  ggplot(tidy_rank_taxa_proportions_tab) +
  geom_col(aes(x=sampleName2, y=abundance, fill=taxa), position = "stack") +
  facet_grid(~ Fase, scales = "free", space = "fixed") +
  theme_bw(base_size = 12) +
  labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
  theme(aspect.ratio = 1.5, panel.spacing = unit(-.1, "lines"), text=element_text(family="serif", color = "black", face = "plain"),
        strip.placement = "outside", strip.background = element_blank(), strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12, angle = 30, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 13, hjust = 1), 
        legend.text = element_text(size = 9), legend.title = element_text(size = 12), legend.key.size = unit(0.7, "lines")) + 
  scale_fill_manual(values = cores_taxa, guide = guide_legend(ncol = 1)) +
  scale_y_continuous(expand = expansion(mult = c(.01, .01)), labels = function(x) paste0(x, "%")) +
  scale_x_discrete(expand = expansion(mult = c(0.08, 0.08)))
#print(barra_phylum)
#ggsave('figuras/barra_phylum.tiff', plot = barra_phylum, height = 1500, width = 2800, units = "px", dpi = 300)

                     
## Criar o gráfico de heatmap para phylum ------
tidy_rank_taxa_proportions_tab_noUnclassOthers <- 
  tidy_rank_taxa_proportions_tab[!tidy_rank_taxa_proportions_tab$taxa %in% c("Unclassified","Others"),]

tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa <- 
  factor(tidy_rank_taxa_proportions_tab_noUnclassOthers$taxa,levels=rev(plot_labels))

plota_phylum <-
  ggplot(tidy_rank_taxa_proportions_tab_noUnclassOthers) +
  geom_tile(aes(x = sampleName2, y = taxa, fill = abundance), color = "gray98") +
  scale_fill_gradient(high = "#44b19f", low = "gray98", guide = guide_colorbar(barwidth = 1, barheight = 5)) +
  theme_bw(base_size = 12) +
  labs(y = str_to_title(taxrank), fill = "Abundance (%)") +
  theme(aspect.ratio = 1.5, text = element_text(family = "serif", color = "black",  face = "plain"),
    strip.placement = "outside", panel.spacing = unit(.3, "lines"),
    strip.background = element_blank(), strip.text = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10), legend.title = element_text(size = 11),
    axis.text.y = element_text(size = 11, hjust = 1),    
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1, color = "black"), axis.title.x = element_blank(),
    panel.grid.minor = element_line(color = "white", size = 0.1),
    panel.border = element_rect(color = "gray40", size =0.8),
    panel.grid.major = element_blank()) +
  scale_y_discrete(position = "left", guide = guide_axis(position = "left"), 
                   expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  facet_grid(~ Fase, scales = "free", space = "fixed")
#print(plota_phylum)
#ggsave('figuras/bacteria_arc/heat_phylum.tiff', plot = plota_phylum, height = 1500, width = 2800, units = "px", dpi = 300)

                     
# Layout: Plotar todos os graficos de barra ------
# Organizar os gráficos em um layout de 2 colunas e 3 linhas
final_plot <- ggarrange(
  ggarrange(barra_phylum, plota_class, labels = c("A", "B"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_order, plota_family, labels = c("C", "D"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_genus, labels = c("E", "F"), ncol = 2, widths = c(1, 1)),
  nrow = 3, align = "hv") +bgcolor("white")
print(final_plot)
#ggsave('figuras/taxo_5_2/layout_barra_2.tiff', plot = final_plot, height = 2999, width = 2800, units = "px", dpi = 250)

# Layout: Plotar todos heatmaps -------
final_plot_h <- ggarrange(
  ggarrange(plota_phylum, plota_class_h, labels = c("A", "B"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_order_h, plota_family_h, labels = c("C", "D"), ncol = 2, widths = c(1, 1)),
  ggarrange(plota_genus_h, labels = c("E", "F"), ncol = 2, widths = c(1, 1)),
  nrow = 3, align = "hv") + bgcolor("white")
print(final_plot_h)
#ggsave('figuras/taxo_5_2/final_plot_h_2.tiff', plot = final_plot_h,height = 2800, width = 3500, units = "px", dpi = 290)



