#Perfil funcional
# Análise Funcional usando PICRUSt2. Prever funções metabólicas com base nos dados taxonômicos e gerar gráficos e tabelas de abundância funcional.

library(tidyverse)
library(phyloseq)
library(readr)
library(ggsci)
library(writexl)
library(ggpubr)
#install.packages("writexl")
# Preparo de Dados para PICRUSt e Execução do PICRUSt2: As sequências ASV são convertidas em um arquivo .fasta.
# As contagens ASV são convertidas em formato .biom para uso no PICRUSt2.
# O comando picrust2_pipeline.py é utilizado para gerar previsões funcionais a partir dos arquivos .fasta e .biom.

ps <- readRDS("/Users/LS28_Seq/Documents/ref/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/ref/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")
order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)

# Run picrust
arquivo_fasta <- "ASVseqs.fasta"
ifelse(file.exists(arquivo_fasta),
       file.remove(arquivo_fasta),
       "Nada a remover")
# 
seqs_fasta <-
  seqtab %>%
  refseq() %>%
  data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename("ASV" = rowname, "FASTA" = ".")
# 
for (i in 1:nrow(seqs_fasta)) {
  cat(
    paste(">", seqs_fasta$ASV[i], sep = ""),
    file = arquivo_fasta,
    sep = "\n",
    append = TRUE
  )
  cat(seqs_fasta$FASTA[i],
      file = arquivo_fasta,
      sep = "\n",
      append = TRUE)
}
# 
asv_table <-
  seqtab %>%
  otu_table(taxa_are_rows = FALSE) %>%
  data.frame()
# 
biomformat::make_biom(t(asv_table)) %>%
  biomformat::write_biom(x = ., biom_file = "ASVcounts.biom")

# RUN PICRUST: Seleção de KOs (KEGG Orthology):
# Define conjuntos de genes de interesse para processos metabólicos, como nitrificação e síntese de glicogênio (ex.: AOB, NOB, etc.).

#conda activate picrust2
#picrust2_pipeline.py -s ASVseqs.fasta -i ASVcounts.biom -o output_picrust -p 6 --stratified

AOB <- c("K10944", "K10945", "K10946", "K10535", "K05601", "K15864") 	# ammonia monooxygenase;AMO EC 1.14.99.39, Hydroxylamine oxidoreductase, (HAO) EC 1.7.2.6 -> hydroxylamine oxidase
NOB <- c("K00370", "K00371") # Nitrite oxidoreductase (NOR or NXR)
ANAMMOX <- c("K20932", "K20933", "K20934", "K20935") # hydrazine hydrolase, hydrazine dehydrogenase
GAO <- c("K20812", "K00975", "K00688", "K02438") # glycogen synthase glgA, glucose-1-phosphate adenylyltransferas glgC, glycogen phosphorylase glgP, glycogen debranching enzyme glgX
PAO <- c("K00937", "K22468") # polyphosphate kinase ppk, polyphosphate kinase ppk2
DNB <- c("K00372", "K00360", "K00367", "K00370", "K00371", "K00373", "K00374", "K10534") # assimilatory nitrate reductase catalytic subunit nasA, assimilatory nitrate reductase electron transfer subunit nasB

ko <- c(AOB,NOB,ANAMMOX,GAO,PAO,DNB)
funcao <- c(rep("AOB",6),rep("NOB",2),rep("ANAMMOX",4),rep("GAO",4),rep("PAO",2),rep("DNB",8))
KOs <- data.frame(funcao,ko)

exclude_asvs <- c(
  as.vector(na.omit(rownames(tax_table(seqtab)[tax_table(seqtab)[,"order"] == "Chloroplast",]))),
  as.vector(na.omit(rownames(tax_table(seqtab)[tax_table(seqtab)[,"family"] == "Mitochondria",])))
)

taxons <- data.frame(tax_table(seqtab))
taxons$ASV <- rownames(taxons)

dados <- data.frame(sample_data(seqtab))

ko_table <- read.table("/Users/LS28_Seq/Documents/JM_03.10.24/artigo2/output_picrust/KO_metagenome_out/pred_metagenome_contrib.tsv.gz",fill=T,header=T)

merged <- ko_table %>%
  merge(KOs,by.x=c('function.'),by.y=c('ko')) %>%
  filter(!taxon %in% exclude_asvs) %>%
  merge(taxons,by.x=c('taxon'),by.y=c('ASV')) %>%
  merge(dados,by.x=c("sample"),by.y=c("sampleName")) %>%
  group_by(sampleName2,funcao) %>%
  summarize(abundance = sum(taxon_function_abun)) %>%
  drop_na()

merged$funcao <- factor(merged$funcao, levels = c("AOB", "NOB", "DNB", "GAO", "PAO"))

#Gráfico de Funções (barra): ggplot() é utilizado para criar gráficos de barras, 
#mostrando a abundância funcional agrupada por processo.

ggplot(merged, aes(x = sampleName2, y = abundance, fill = funcao)) +
  geom_bar(stat = "identity", width = 0.85, position=position_dodge()) +
  facet_grid(cols = vars(funcao), switch = "y") +
  theme_classic(base_size = 14) +
  theme(text = element_text(family = "serif", color = "black", face = "plain"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "gray90"),
        axis.title.x = element_blank()) +
  labs(y = "Functional abundance",
       fill = "Abundance") +
  scale_fill_d3(palette = "category20b")

#ggsave(filename="figuras/func.tiff", width = 1500, height = 900, units = "px", dpi = 200)

tabela_funcs <- 
  ko_table %>%
  merge(KOs,by.x=c('function.'),by.y=c('ko')) %>%
  filter(!taxon %in% exclude_asvs) %>%
  merge(taxons,by.x=c('taxon'),by.y=c('ASV')) %>%
  merge(dados,by.x=c("sample"),by.y=c("sampleName")) %>%
  group_by(sample,funcao,domain,phylum,class,order,family,genus) %>%
  summarize(abundance = sum(taxon_function_abun)) %>%
  select("Domínio"=domain,"Filo"=phylum,"Classe"=class,"Ordem"=order,"Família"=family,"Gênero"=genus,"Amostra"=sample,"Função"=funcao,"Abundância"=abundance)

write_xlsx(tabela_funcs,"tabelas/functions_table.xlsx")

library(scales) # Para scale_fill_d3

# Filtrar a tabela para remover filas com Filo NA
tabela_funcs_filtrada <- tabela_funcs %>%
  filter(!is.na(Filo)) %>%
  group_by(Filo, Amostra, Função) %>%
  summarise(Abundância = sum(Abundância), .groups = 'drop') %>%
  arrange(desc(Abundância)) %>%
  slice_head(n = 55)

ggplot(tabela_funcs_filtrada, aes(x = Amostra, y = Abundância, fill = Filo)) +
  geom_bar(stat = "identity", width = 0.85, position=position_dodge()) +
  facet_grid(cols = vars(Função), switch = "y") +
  theme_classic(base_size = 14) +
  theme(text = element_text(family = "serif", color = "black", face = "plain"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "gray98"),
        panel.grid.major.y = element_line(colour = "gray80"),
        axis.title.x = element_blank()) +
  labs(y = "Functional abundance",
       fill = "Phylum contribution") +
  scale_fill_d3(palette = "category20b") + 
  scale_x_discrete(labels = c("JM3I" = "I", 
                              "JM3R1" = "R1", 
                              "JM3R2" = "R2", "JM3R3" = "R3"))

ggsave(filename="figuras/func_filo_10.tiff", width = 1700, height = 900, units = "px", dpi = 200)


library(dplyr)
library(ggplot2)

# Filtrar a tabela para remover filas com Família NA
tabela_funcs_filtrada <- tabela_funcs %>%
  filter(!is.na(Família)) %>%
  group_by(Família, Amostra, Função) %>%
  summarise(Abundância = sum(Abundância), .groups = 'drop') %>%
  arrange(desc(Abundância)) %>%
  slice_head(n = 17)

ggplot(tabela_funcs_filtrada, aes(x = Amostra, y = Abundância, fill = Família)) +
  geom_bar(stat = "identity", width = 0.85, position = position_dodge()) +
  facet_grid(cols = vars(Função), switch = "y") +
  theme_classic(base_size = 14) +
  theme(text = element_text(family = "serif", color = "black", face = "plain"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "gray98"),
        panel.grid.major.y = element_line(colour = "gray80"),
        axis.title.x = element_blank()) +
  labs(y = "Functional abundance",
       fill = "Family contribution") +
  scale_fill_d3(palette = "category20b") + 
  scale_x_discrete(labels = c("JM3I" = "I", 
                              "JM3R1" = "R1", 
                              "JM3R2" = "R2", 
                              "JM3R3" = "R3"))

#ggsave(filename="figuras/func_family_10.tiff", width = 1700, height = 900, units = "px", dpi = 200)


