#!/bin/bash

# Script para análise de dados de amplicon usando QIIME2

# 2. Análise de Qualidade Inicial
echo "Ativando ambiente conda e importando dados..."
source ~/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

qiime tools import \
   --type 'SampleData[PairedEndSequencesWithQuality]' \
   --input-path raw \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt \
   --output-path inter_files/demux-paired-end.qza

# 2.1. Executar análise de qualidade inicial
echo "Gerando relatório de qualidade inicial..."
qiime demux summarize --i-data inter_files/demux-paired-end.qza --o-visualization quality_reports/demux_summary.qzv

# 3. Trimagem e Filtragem dos Dados
echo "Executando trimagem e filtragem..."
qiime dada2 denoise-paired \
--i-demultiplexed-seqs inter_files/demux-paired-end.qza \
--p-trim-left-f 2 \
--p-trim-left-r 2 \
--p-trunc-len-f 140 \
--p-trunc-len-r 150 \
--p-trunc-q 7 \
--o-representative-sequences inter_files/asv-sequences-0.qza \
--o-table inter_files/feature-table-0.qza \
--o-denoising-stats inter_files/dada2-stats.qza \
--verbose

# 4. Análise de Qualidade Pós-Trimagem
echo "Gerando relatório de estatísticas de denoising..."
qiime metadata tabulate \
--m-input-file inter_files/dada2-stats.qza \
--o-visualization quality_reports/dada2-stats.qzv

# 5. Tabela de Frequências de Sequências
echo "Gerando tabelas de frequência e sumários..."
qiime metadata tabulate --m-input-file metadados1.txt --o-visualization sample-metadata.qzv

qiime feature-table summarize \
--i-table inter_files/feature-table-0.qza \
--m-sample-metadata-file metadados1.txt \
--o-visualization quality_reports/feature-table-0-summ.qzv

qiime feature-table tabulate-seqs --i-data inter_files/asv-sequences-0.qza --o-visualization quality_reports/asv-sequences-0-summ.qzv

# Classificação taxonômica
echo "Classificando ASVs usando banco de dados SILVA..."
qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-nb-classifier.qza \
--i-reads inter_files/asv-sequences-0.qza \
--o-classification inter_files/taxonomy.qza

# Visualizar taxonomia
echo "Gerando visualização da taxonomia..."
qiime metadata tabulate \
--m-input-file inter_files/taxonomy.qza \
--o-visualization quality_reports/taxonomy.qzv

# Barplot taxonômico
echo "Criando barplot taxonômico..."
qiime taxa barplot \
--i-table inter_files/feature-table-0.qza \
--i-taxonomy inter_files/taxonomy.qza \
--m-metadata-file metadados1.txt \
--o-visualization quality_reports/taxa-bar-plots.qzv

# Exportação dos dados
echo "Exportando tabelas para formatos TSV..."
sudo /home/ls28/anaconda3/envs/qiime2-amplicon-2024.10/bin/qiime tools export \
--input-path inter_files/taxonomy.qza \
--output-path inter_files/exports

sudo /home/ls28/anaconda3/envs/qiime2-amplicon-2024.10/bin/qiime tools export \
--input-path inter_files/feature-table-0.qza \
--output-path inter_files

biom convert -i inter_files/feature-table.biom -o inter_files/feature-table-0.tsv --to-tsv

echo "Pipeline concluído com sucesso!"
