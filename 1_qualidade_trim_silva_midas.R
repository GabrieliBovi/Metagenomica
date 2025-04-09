#install.packages(c("tidyverse", "knitr", "BiocStyle", "ggplot2", "gridExtra"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("dada2", "phyloseq", "DECIPHER"))
#  if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install("BiocStyle")

#Dada2 para processamento de sequências; phyloseq para análise de dados de microbioma; DECIPHER para atribuição taxonômica ----
library(stats)
.cran_packages <- c("tidyverse","knitr","BiocStyle","ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER")
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
library(dada2)
library(stats)
library(parallel)

#1. Diretórios
##Define as pastas onde estão os arquivos FASTQ e onde os resultados intermediários e arquivos filtrados serão salvos.
fastq_folder <- "/Users/LS28_Seq/Documents/ref/raw"
filter_folder <- "/Users/LS28_Seq/Documents/ref/filtered_reads"
intermediate_folder <- "/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/inter_files"

filt_path <- file.path(filter_folder)
if(!file_test("-d", filt_path)) dir.create(filt_path)
intermediate_path <- file.path(intermediate_folder)
if(!file_test("-d", intermediate_path)) dir.create(intermediate_path)

#2. Lista os arquivos FASTQ (leitura forward (fnFs) e reverse(fnRs)) e extrai os nomes das amostras.
fnFs <- sort(list.files(fastq_folder, pattern="...R1_001.fastq.gz"))
fnRs <- sort(list.files(fastq_folder, pattern="...R2_001.fastq.gz"))
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(fastq_folder, fnFs)
fnRs <- file.path(fastq_folder, fnRs)

#3. Filtragem das amostras
##Função plotQualityProfile do pacote DADA2 para gerar gráficos da qualidade das leituras
plotQualityProfile(fnFs[1:5]) # 1000 x 600
plotQualityProfile(fnRs[1:5])

filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
names(filtFs) <- sampleNames
names(filtRs) <- sampleNames 

print("Filtering...")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0,                  ##remove leituras contendo Ns 
                     maxEE=c(2,2),              ##Define o número máximo permitido de Expected Errors
                     truncQ= c(15,10),                ##trunca as leituras ande a qualidade cai abaixo de 2
                     trimLeft = c(15, 20),        # Remove os ... primeiros pares de bases (F e R)
                     truncLen = c(260, 240), ##trunca as leituras ande a qualidade cai
                     rm.phix=TRUE,             ##remove contaminante
                     compress=TRUE,
                     multithread=TRUE,
                     verbose = TRUE)

saveRDS(out,paste0(intermediate_folder,"/out.rds"))

plotQualityProfile(filtFs[1:5])
plotQualityProfile(filtRs[1:5])

# 4. Aprendizado de erros e processamento de leituras: Usa learnErrors para aprender os padrões de erro nas leituras filtradas.
## Isso é parte crucial da correção de erros no DADA2.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

saveRDS(errF,paste0(intermediate_folder,"/errF.rds"))
saveRDS(errR,paste0(intermediate_folder,"/errR.rds"))

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

saveRDS(dadaFs,paste0(intermediate_folder,"/dadaFs.rds"))
saveRDS(dadaRs,paste0(intermediate_folder,"/dadaRs.rds"))

#5. Denoising e fusão de pares: 
##mergePairs combina leituras forward e reverse.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=4) # o ideal seria minimo de 12, mas como é esgoto podemos tentar com 10
seqtabAll <- makeSequenceTable(mergers)
rownames(seqtabAll) <- sampleNames
lengthsDist <- table(nchar(getSequences(seqtabAll)))

saveRDS(mergers,paste0(intermediate_folder,"/mergers.rds"))
saveRDS(lengthsDist,paste0(intermediate_folder,"/lengthsDist.rds"))

#6. Remoção de quimeras: 
##removeBimeraDenovo remove quimeras, que são artefatos comuns em sequenciamento de alto rendimento.
seqtabNoC <- removeBimeraDenovo(seqtabAll, method="consensus", verbose = TRUE, multithread=8)

#7. Rastreamento de progresso: 
##Cria uma tabela de resumo mostrando o número de leituras em cada etapa do processo.
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtabNoC))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
track_df <- as.data.frame(track)
track_df$nonchim_percent <- track_df$nonchim * 100 / track_df$input

save(track_df, file=paste0(intermediate_folder,"/track_df.Rda"))
save(seqtabNoC, file=paste0(intermediate_folder,"/seqtabNoC.rds"))

ranks <- c("domain","phylum","class","order","family","genus","species")

#8. Atribuição taxonômica: 
##Usando SILVA: Usa o banco de dados SILVA para atribuir as sequências a taxonomias.
dna <- DNAStringSet(getSequences(seqtabNoC)) # Create a DNAStringSet from the ASVs
load("/Users/LS28_Seq/Documents/SILVA_SSU_r138_2019.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=8, verbose=TRUE)

taxid_silva <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid_silva) <- ranks; rownames(taxid_silva) <- getSequences(seqtabNoC)

saveRDS(taxid_silva, file=paste0(intermediate_folder,"/taxid_silva.rds"))

## Criação de objetos Phyloseq - relatório. Cria um objeto phyloseq para análise de microbioma usando as tabelas de OTU ou ASV e taxonomia.
ps_silva <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
                     tax_table(taxid_silva))

dna_silva <- Biostrings::DNAStringSet(taxa_names(ps_silva))
names(dna_silva) <- taxa_names(ps_silva)
ps_silva <- merge_phyloseq(ps_silva, dna_silva)
taxa_names(ps_silva) <- paste0("ASV", seq(ntaxa(ps_silva)))

saveRDS(ps_silva, paste0(intermediate_folder,"/ps_silva.rds"))
print(paste0("Done! Phyloseq object saved to >> ",intermediate_folder,"/ps_silva.rds <<"))


##Usando MiDAS: Atribui taxonomia usando o banco de dados MiDAS.
##MiDAS 5.3 - demorou 7 horas meu deus, e não terminou  

set.seed(100) # Initialize random number generator for reproducibility
taxid_midas <- assignTaxonomy(seqtabNoC, "/Users/LS28_Seq/Documents/DADA2_taxonomy.fa MiDAS 5.3.fa", multithread=TRUE)
colnames(taxid_midas) <- ranks

saveRDS(taxid_midas, file=paste0(intermediate_folder,"/taxid_midas.rds"))

###Criação de objetos Phyloseq - relatório
###Cria um objeto phyloseq para análise de microbioma usando as tabelas de OTU ou ASV e taxonomia. MiDAS.
ps_midas <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
                     tax_table(taxid_midas))

dna_midas <- Biostrings::DNAStringSet(taxa_names(ps_midas))
names(dna_midas) <- taxa_names(ps_midas)
ps_midas <- merge_phyloseq(ps_midas, dna_midas)
taxa_names(ps_midas) <- paste0("ASV", seq(ntaxa(ps_midas)))

saveRDS(ps_midas, paste0(intermediate_folder,"/ps_midas.rds"))
print(paste0("Done! Phyloseq object saved to >> ",intermediate_folder,"/ps_midas.rds <<"))

#10. Finalização do Silva e Midas



## Tabela de ASVs, rarecurve e indices. ----
# Tabela taxons: Tabela de OTUs (ASVs) com suas taxonomias e abundâncias relativas
library(phyloseq)
library(reshape2)
library(writexl)
library(dplyr)
library(vegan)
library(tidyverse)
library(devtools)
library(QsRutils)

#1. Leitura dos dados e normalização:	otu_table(seqtab): Obtém as abundâncias relativas das ASVs.
# apply() calcula as proporções relativas (%).
ps <- readRDS("/Users/LS28_Seq/Documents/JM_03.10.24//artigo1/inter_files/ps_silva.rds")
samdf <- read.table("/Users/LS28_Seq/Documents/JM_03.10.24/artigo1/metadados.txt", header=TRUE, encoding="UTF-8", sep="\t")
order_samples <- paste0(samdf$sampleName)
seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)

asv_proportions_tab <- apply(otu_table(seqtab), 1, function(x) x/sum(x)*100)
ASVs <- asv_proportions_tab %>% as.data.frame()
taxa <- tax_table(seqtab) %>% as.data.frame()
ASVs$ASV <- rownames(ASVs)
tidy <- ASVs %>% melt()

# Criação da Tabela e Exportação com as abundâncias das ASVs. merge() combina as informações taxonômicas (filo, família, etc.).
taxa$ASV <- rownames(taxa)
merged <- merge(taxa,tidy,by='ASV')
colnames(merged) <- c("ASV", "Domínio", "Filo", "Classe", "Ordem", "Família", "Gênero", "Espécie", "Amostra", "Abundância relativa (%)")

write_xlsx(merged,"tabelas/ASV_table_bruto.xlsx")

#2. Grafico de rarecurve (5x4 pdf)
#tiff("figuras/rarecurve_plot.tiff", width = 1000, height =1000, res = 300)
sample_names <- sample_data(seqtab)[[2]]
par(cex.axis =0.6, cex.lab =0.6)  # Ajuste conforme necessári
rarecurve_data <-rarecurve(otu_table(seqtab) %>% data.frame,
                           step = 500, cex=0.6,xlab = "Sample size",
                           ylab = "Species", label = "false", bty = "L", family = "serif")

for (i in seq_along(rarecurve_data)) { # Adicionar os rótulos manualmente ao final de cada curva
  x_values <- attr(rarecurve_data[[i]], "Subsample")
  y_values <- rarecurve_data[[i]]
  label_x <- x_values[length(x_values)]
  label_y <- y_values[length(y_values)]
  rect(xleft = label_x + - 3000, ybottom = label_y - 6, # Adicionar a borda preta (caixa) ao redor do texto
       xright = label_x - 0, ytop = label_y + 6, 
       border = "black", col = "white")  
  text(x = label_x, y = label_y,   # Adicionar o nome da amostra
       labels = sample_names[i], pos = 2, cex = 0.5,
       col = "black", family = "serif") 
}
                             
#dev.off()

#3. Apresentar tabela com índices de Chao 1, Shannon (Motteran et al., 2018). Calcular e gerar uma tabela com índices de diversidade alfa, como Chao1, Shannon, e Simpson.
#Calcular Índices: estimate_richness(): Calcula vários índices de diversidade (Chao1, Shannon, Simpson, etc.) para cada amostra.
a <- goods(otu_table(seqtab)) %>% rownames_to_column(var = "Samples")
b <- estimate_richness(seqtab, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) %>% select(-c("se.chao1", "se.ACE"))
indices <- bind_cols(a, b)
indices_tabela <- indices[,c("Samples", "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "no.seqs")]
write_xlsx(indices_tabela,"tabelas/indices_alpha_diversity.xlsx")
