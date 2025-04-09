
library(stats)
.cran_packages <- c("tidyverse","knitr","BiocStyle","ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER")
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
library(Biostrings)
library(taxize)
library(ShortRead)
library(DECIPHER)


# Definição dos diretórios
filter_folder <- "/Users/LS28_Seq/Documents/raw_CB_nanopore/CB_nanopore/trimagem_q9_cut17_min100"
intermediate_folder <- "/Users/LS28_Seq/Documents/raw_CB_nanopore/CB_nanopore/inter_files"

dir.create(filter_folder, showWarnings = FALSE)
dir.create(intermediate_folder, showWarnings = FALSE)

# Lista de amostras e arquivos
sample_names <- c("CB1", "CB2") # "CB1", "CB2", "CB3", "CB4",  "CB7", "CB8", "CB9"

fnFiles <- file.path(filter_folder, paste0(sample_names, "_trimmed.fastq"))

# Carrega as sequências processadas
seqs <- getSequences(fnFiles)  

# Converter FASTQ para FASTA
convert_fastq_to_fasta <- function(fastq_file) {
  fq <- readFastq(fastq_file)  # Lê o FASTQ
  seqs <- sread(fq)  # Extrai as sequências
  
  # Define o nome do arquivo de saída (FASTA)
  fasta_file <- sub(".fastq", ".fasta", fastq_file)
  
  # Salva como FASTA
  writeXStringSet(seqs, fasta_file, format="fasta")
  return(fasta_file)
}

# Aplicar a conversão em todos os arquivos
fasta_files <- sapply(fnFiles, convert_fastq_to_fasta)

# Cria um conjunto de sequências DNAStringSet
seqs <- readDNAStringSet(fasta_files)

# Carrega o banco de dados treinado do SILVA
load("/Users/LS28_Seq/Documents/SILVA_SSU_r138_2019.RData")  

# Atribui taxonomia usando SILVA e Extrai a taxonomia formatada
library(parallel)

num_cores <- detectCores() - 1
ids <- IdTaxa(seqs, trainingSet, strand="top", processors=num_cores, verbose=TRUE)

ranks <- c("domain","phylum","class","order","family","genus","species")
taxid_silva <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA  # Remove classificações não resolvidas
  taxa
}))
colnames(taxid_silva) <- ranks
rownames(taxid_silva) <- getSequences(fnFiles)

# Salva os resultados
saveRDS(taxid_silva, file=paste0(intermediate_folder, "/ps_silva.rds"))

# Exibir uma prévia
head(taxid_silva)














# Rastreamento de progresso (contagem de reads por amostra)
track_df <- data.frame(
  sample = sample_names,
  input_reads = sapply(seq_list, length)  # Conta sequências por amostra
)

print(track_df)

# Carregando o banco de dados SILVA
silva_db <- "/Users/LS28_Seq/Documents/SILVA_SSU_r138_2019.RData"
load(silva_db)

# Atribuição taxonômica usando SILVA
dna <- DNAStringSet(getSequences(seqtabAll)) 
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=TRUE)

ranks <- c("domain","phylum","class","order","family","genus","species")

taxid_silva <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid_silva) <- ranks
rownames(taxid_silva) <- getSequences(seqtabAll)

# Salvando os resultados
saveRDS(taxid_silva, file=paste0(intermediate_folder,"/taxid_silva.rds"))

# Criando o objeto phyloseq
ps_silva <- phyloseq(
  otu_table(seqtabAll, taxa_are_rows=FALSE), 
  tax_table(taxid_silva)
)

# Adicionando nomes das ASVs
dna_silva <- Biostrings::DNAStringSet(taxa_names(ps_silva))
names(dna_silva) <- taxa_names(ps_silva)
ps_silva <- merge_phyloseq(ps_silva, dna_silva)
taxa_names(ps_silva) <- paste0("ASV", seq(ntaxa(ps_silva)))

saveRDS(ps_silva, paste0(intermediate_folder,"/ps_silva.rds"))
print(paste0("Done! Phyloseq object saved to >> ", intermediate_folder, "/ps_silva.rds <<"))





##Usando MiDAS: Atribui taxonomia usando o banco de dados MiDAS. - MiDAS 5.3.   

set.seed(100) # Initialize random number generator for reproducibility
taxid_midas <- assignTaxonomy(seqtabNoC, "/Users/LS28_Seq/Documents/DADA2_taxonomy.fa MiDAS 5.3.fa", multithread=TRUE)#.gz
colnames(taxid_midas) <- ranks

saveRDS(taxid_midas, file=paste0(intermediate_folder,"/taxid_midas.rds"))

###Criação de objetos Phyloseq - relatório - Cria um objeto phyloseq para análise de microbioma usando as tabelas de OTU ou ASV e taxonomia. MiDAS.

ps_midas <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
                     tax_table(taxid_midas))

dna_midas <- Biostrings::DNAStringSet(taxa_names(ps_midas))
names(dna_midas) <- taxa_names(ps_midas)
ps_midas <- merge_phyloseq(ps_midas, dna_midas)
taxa_names(ps_midas) <- paste0("ASV", seq(ntaxa(ps_midas)))

saveRDS(ps_midas, paste0(intermediate_folder,"/ps_midas.rds"))
print(paste0("Done! Phyloseq object saved to >> ",intermediate_folder,"/ps_midas.rds <<"))

#10. Finalização do Silva e Midas



# Tabela de ASVs/ taxons e rarecurve.
# Tabela taxons: Tabela de OTUs (ASVs) com suas taxonomias e abundâncias relativas

library(reshape2)
library(writexl)
library(dplyr)
library(vegan)
library(devtools)
library(QsRutils)

#1. Leitura dos dados e normalização:	otu_table(seqtab): Obtém as abundâncias relativas das ASVs.
# apply() calcula as proporções relativas (%).

ps <- readRDS("/Users/LS28_Seq/Documents/nano_24.10.24/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/nano_24.10.24/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")
order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)

asv_proportions_tab <- apply(otu_table(seqtab), 1, function(x) x/sum(x)*100)
ASVs <- asv_proportions_tab %>% as.data.frame()
taxa <- tax_table(seqtab) %>% as.data.frame()
ASVs$ASV <- rownames(ASVs)
tidy <- ASVs %>% melt()

#2. Criação da Tabela e Exportação com as abundâncias das ASVs -----
# merge() combina as informações taxonômicas (filo, família, etc.). write_xlsx() salva a tabela resultante em um arquivo Excel no diretório especificado.

taxa$ASV <- rownames(taxa)
merged <- merge(taxa,tidy,by='ASV')
colnames(merged) <- c("ASV", "Domínio", "Filo", "Classe", "Ordem", "Família", "Gênero", "Espécie", "Amostra", "Abundância relativa (%)")

write_xlsx(merged,"tabelas/ASV_table_bruto.xlsx")


#3. Grafico de rarecurve (5x4 pdf)

#tiff("figuras/rarecurve_plot.tiff", width = 1000, height =1000, res = 300)

sample_names <- sample_data(seqtab)[[2]]

par(cex.axis =0.9, cex.lab =0.9)  # Ajuste conforme necessári

rarecurve_data <-rarecurve(otu_table(seqtab) %>% data.frame,
                           step = 5, cex=0.6,xlab = "Sample size",
                           ylab = "Species", label = "true", bty = "L", family = "serif")

for (i in seq_along(rarecurve_data)) { # Adicionar os rótulos manualmente ao final de cada curva
  x_values <- attr(rarecurve_data[[i]], "Subsample")
  y_values <- rarecurve_data[[i]]
  label_x <- x_values[length(x_values)]
  label_y <- y_values[length(y_values)]
  rect(xleft = label_x + - 3000, ybottom = label_y - 2, # Adicionar a borda preta (caixa) ao redor do texto
       xright = label_x - 0, ytop = label_y + 2, 
       border = "black", col = "white")  
  text(x = label_x, y = label_y,   # Adicionar o nome da amostra
       labels = sample_names[i], pos = 2, cex = 0.5,
       col = "black", family = "serif") 
}
#dev.off()
#dev.new()
#dev.next()

#4. Apresentar tabela com índices de Chao 1, Shannon (Motteran et al., 2018)
#Calcular e gerar uma tabela com índices de diversidade alfa, como Chao1, Shannon, e Simpson.
#Calcular Índices: estimate_richness(): Calcula vários índices de diversidade (Chao1, Shannon, Simpson, etc.) para cada amostra.
#A tabela resultante é exportada com write_xlsx().

a <- goods(otu_table(seqtab)) %>% rownames_to_column(var = "Samples")
b <- estimate_richness(seqtab, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) %>% select(-c("se.chao1", "se.ACE"))
indices <- bind_cols(a, b)
indices_tabela <- indices[,c("Samples", "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "no.seqs")]
write_xlsx(indices_tabela,"tabelas/indices_alpha_diversity.xlsx")



