#1. Configuração do Ambiente
#Certifique-se de que o QIIME 2 esteja instalado e que as ferramentas CLI estejam funcionando 
#no seu sistema. No RStudio, usaremos o pacote system2 para chamar os comandos do QIIME 2.

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R", force = TRUE)

# Carregar pacotes necessários
.cran_packages <- c("tidyverse", "qiime2R", "ggplot2", "writexl", "reshape2", "dplyr", "vegan")
sapply(.cran_packages, require, character.only = TRUE)


#2. Análise de Qualidade Inicial
#conda activate qiime2-metagenome-2024.10
conda activate qiime2-amplicon-2024.10

qiime tools import \
   --type 'SampleData[PairedEndSequencesWithQuality]' \
   --input-path raw \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --output-path inter_files/demux-paired-end.qza

# 2.1. Executar análise de qualidade inicial
 qiime demux summarize --i-data inter_files/demux-paired-end.qza --o-visualization quality_reports/demux_summary.qzv

# https://view.qiime2.org - quality_reports/demux_summary.qzv

#3. Trimagem e Filtragem dos Dados
 conda deactivate
 conda activate qiime2-amplicon-2024.10
  
# Executar trimagem e filtragem com DADA2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs inter_files/demux-paired-end.qza \
--p-trim-left-f 5 \
--p-trim-left-r 5 \
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--p-trunc-q 15 \
--o-representative-sequences inter_files/asv-sequences-0.qza \
--o-table inter_files/feature-table-0.qza \
--o-denoising-stats inter_files/dada2-stats.qza
--verbose


#4. Análise de Qualidade Pós-Trimagem
# Gerar relatório de estatísticas de denoising
qiime metadata tabulate \
--m-input-file inter_files/dada2-stats.qza \
--o-visualization quality_reports/dada2-stats.qzv


#5. Tabela de Frequências de Sequências
# Gerar tabela de frequências. The feature table describes which amplicon sequence variants (ASVs) were observed in which samples, and how many times each ASV was observed in each sample. 
qiime metadata tabulate --m-input-file sample_metadata.txt --o-visualization sample-metadata.qzv # verificar se o metadado esta certo

qiime feature-table summarize \
--i-table inter_files/feature-table-0.qza \
--m-sample-metadata-file sample_metadata.txt \ 
--o-visualization quality_reports/feature-table-0-summ.qzv

qiime feature-table tabulate-seqs --i-data inter_files/asv-sequences-0.qza --o-visualization quality_reports/asv-sequences-0-summ.qzv


#6. Remoção de Quimeras
#No QIIME 2, isso já é realizado automaticamente no passo de trimagem com DADA2. Nenhum comando adicional é necessário para esta etapa.


#7. Atribuição Taxonômica 
#Opção: Re-treinar seu próprio classificador: Se precisar manter esse classificador, pode re-treiná-lo com a versão atual do scikit-learn. Esse processo é mais demorado e envolve baixar a base de dados Silva e rodar o comando:
#qiime feature-classifier fit-classifier-naive-bayes \
#--i-reference-reads silva-138-99-seqs.qza \
#--i-reference-taxonomy silva-138-99-tax.qza \
#--o-classifier silva-138-99-nb-classifier-new.qza

# Caminho do banco de dados SILVA   
qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-nb-classifier.qza \
--i-reads inter_files/asv-sequences-0.qza \
--o-classification inter_files/taxonomy.qza

# Visualizar taxonomia
qiime metadata tabulate \
--m-input-file inter_files/taxonomy.qza \
--o-visualization quality_reports/taxonomy.qzv

# In order to be able to download the sample OTU table need to do the taxonomy assignment and then make the taxa barplot. Then can download csv file with sequence number, samples and taxonomy.
qiime taxa barplot \
--i-table inter_files/feature-table-0.qza \
--i-taxonomy inter_files/taxonomy.qza \
--m-metadata-file sample_metadata.txt \
--o-visualization quality_reports/taxa-bar-plots.qzv

#Extra bit of code to genera te a taxonomy table table to tsv from the commandline
# Exportar a tabela de ASVs e taxonomia
sudo /home/ls28/anaconda3/envs/qiime2-amplicon-2024.10/bin/qiime tools export \
--input-path inter_files/taxonomy.qza \
--output-path inter_files/exports

sudo /home/ls28/anaconda3/envs/qiime2-amplicon-2024.10/bin/qiime tools export \
--input-path inter_files/feature-table-0.qza \
--output-path inter_files/exports

biom convert -i inter_files/feature-table.biom -o inter_files/exports/feature-table-0.tsv --to-tsv

#sudo /home/ls28/anaconda3/envs/qiime2-amplicon-2024.10/bin/qiime tools export \
#--input-path quality_reports/taxa-bar-plots.qzv \
#--output-path inter_files/exports


# ----- Analise do Silva ------- 

# Carregar a OTU table
otu_table <- read.table("inter_files/exports/feature-table-0.tsv", header = TRUE, sep = "\t", row.names = 1)

# Transpor a tabela de OTUs (amostras nas colunas, ASVs nas linhas)
otu_table_transposed <- t(otu_table)

# Carregar a tabela de taxonomia
tax_table <- read.table("inter_files/exports/taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1)

# Reorganizar tax_table para ter as mesmas IDs que otu_table
tax_table <- tax_table[rownames(otu_table), , drop = FALSE]

# Dividir a taxonomia nos diferentes ranks
tax_table_split <- strsplit(as.character(tax_table$Taxon), ";")

# Definir os ranks de interesse (domínio, filo, classe, ordem, família, gênero, espécie)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Criar uma nova tabela de taxonomia onde cada coluna representa um rank
tax_table_split_df <- do.call(rbind, lapply(tax_table_split, function(x) {
  rank_data <- rep(NA, length(ranks))
  rank_data[1:length(x)] <- gsub("^[a-zA-Z]+__", "", x)  # Remover o prefixo (por exemplo, "d__", "p__")
  return(rank_data)
}))

# Atribuir os nomes das linhas de tax_table_split_df para serem as mesmas de otu_table
rownames(tax_table_split_df) <- rownames(otu_table)

# Converter em um data.frame e atribuir os nomes dos ranks como colunas
tax_table_split_df <- as.data.frame(tax_table_split_df)
colnames(tax_table_split_df) <- ranks

# Transformar tax_table_split_df em uma matriz, como esperado para o phyloseq
tax_table_split_matrix <- as.matrix(tax_table_split_df)

# Criar o objeto phyloseq
library(phyloseq)
library(reshape2)
library(writexl)
library(dplyr)
library(vegan)
library(devtools)
library(QsRutils)

# Criar o objeto phyloseq com a tabela de OTUs transposta
ps_silva <- phyloseq(
  otu_table(otu_table_transposed, taxa_are_rows = FALSE),  # Tabela de OTUs transposta
  tax_table(tax_table_split_matrix)                      # Tabela de taxonomia
)

taxa_names(ps_silva) <- paste0("ASV", seq(ntaxa(ps_silva)))

# Salvar o objeto
saveRDS(ps_silva, "inter_files/ps_silva.rds")


ps <- readRDS("/Users/LS28_Seq/Documents/illumina_artigo/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/illumina_artigo/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")
order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)

asv_proportions_tab <- apply(otu_table(seqtab), 1, function(x) x/sum(x)*100)
ASVs <- asv_proportions_tab %>% as.data.frame()
taxa <- tax_table(seqtab) %>% as.data.frame()
ASVs$ASV <- rownames(ASVs)
tidy <- ASVs %>% melt()

taxa$ASV <- rownames(taxa)
merged <- merge(taxa,tidy,by='ASV')
colnames(merged) <- c("ASV", "Domínio", "Filo", "Classe", "Ordem", "Família", "Gênero", "Espécie", "Amostra", "Abundância relativa (%)")

#write_xlsx(merged,"tabelas/ASV_table_bruto.xlsx")

sample_names <- sample_data(seqtab)[[2]]

par(cex.axis =0.9, cex.lab =0.9)  # Ajuste conforme necessári

rarecurve_data <-rarecurve(otu_table(seqtab) %>% data.frame,
                           step = 500, cex=0.9,xlab = "Sample size",
                           ylab = "Species", label = "true", bty = "L", family = "serif")

for (i in seq_along(rarecurve_data)) { # Adicionar os rótulos manualmente ao final de cada curva
  x_values <- attr(rarecurve_data[[i]], "Subsample")
  y_values <- rarecurve_data[[i]]
  label_x <- x_values[length(x_values)]
  label_y <- y_values[length(y_values)]
  rect(xleft = label_x + - 2000, ybottom = label_y - 10, # Adicionar a borda preta (caixa) ao redor do texto
       xright = label_x - 0, ytop = label_y + 10, 
       border = "black", col = "white")  
  text(x = label_x, y = label_y,   # Adicionar o nome da amostra
       labels = sample_names[i], pos = 2, cex = 0.8,
       col = "black", family = "serif") 
}

a <- goods(otu_table(seqtab)) %>% rownames_to_column(var = "Samples")
b <- estimate_richness(seqtab, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) %>% select(-c("se.chao1", "se.ACE"))
indices <- bind_cols(a, b)
indices_tabela <- indices[,c("Samples", "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "no.seqs")]
write_xlsx(indices_tabela,"tabelas/indices_alpha_diversity_v2.xlsx")





