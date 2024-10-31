library(dplyr)
library(readxl)

# Carregar os dados do arquivo .xlsx
dados <- read_excel("/Users/LS28_Seq/Documents/ref/tabelas/ASV_table_bruto.xlsx")

# Exibir os primeiros registros para verificar se os dados foram carregados corretamente
head(dados)

# Filtrar ASVs que aparecem em todas as 
asvs_comuns <- dados %>% 
  group_by(ASV) %>% 
  summarise(total_amostras = sum(`Abundância relativa (%)` > 0, na.rm = TRUE)) %>% 
  filter(total_amostras == 4) # substitua 4 pelo número de suas amostras


# Filtrar ASVs que aparecem apenas na amostra JM3I
asvs_exclusivos_JM3I <- dados %>%
  group_by(ASV) %>%
  summarise(
    presente_v = sum(`Abundância relativa (%)`[Amostra == "JM3I"] > 0, na.rm = TRUE),
    presente_outros = sum(`Abundância relativa (%)`[Amostra != "JM3I"] > 0, na.rm = TRUE)
  ) %>%
  filter(presente_v > 0 & presente_outros == 0)


# Filtrar ASVs que aparecem apenas na amostra JM3R1
asvs_exclusivos_R1 <- dados %>%
  group_by(ASV) %>%
  summarise(
    presente_v1 = sum(`Abundância relativa (%)`[Amostra == "JM3R1"] > 0, na.rm = TRUE),
    presente_outros = sum(`Abundância relativa (%)`[Amostra != "JM3R1"] > 0, na.rm = TRUE)
  ) %>%
  filter(presente_v1 > 0 & presente_outros == 0)


# Filtrar ASVs que aparecem apenas na amostra JM3R2
asvs_exclusivos_R2 <- dados %>%
  group_by(ASV) %>%
  summarise(
    presente_v2 = sum(`Abundância relativa (%)`[Amostra == "JM3R2"] > 0, na.rm = TRUE),
    presente_outros = sum(`Abundância relativa (%)`[Amostra != "JM3R2"] > 0, na.rm = TRUE)
  ) %>%
  filter(presente_v2 > 0 & presente_outros == 0)


# Filtrar ASVs que aparecem apenas na amostra JM3R2
asvs_exclusivos_R3 <- dados %>%
  group_by(ASV) %>%
  summarise(
    presente_v3 = sum(`Abundância relativa (%)`[Amostra == "JM3R3"] > 0, na.rm = TRUE),
    presente_outros = sum(`Abundância relativa (%)`[Amostra != "JM3R3"] > 0, na.rm = TRUE)
  ) %>%
  filter(presente_v3 > 0 & presente_outros == 0)

# Exibir ASVs exclusivos da amostra
print(asvs_comuns)
print(asvs_exclusivos_JM3I)
print(asvs_exclusivos_R1)
print(asvs_exclusivos_R2)
print(asvs_exclusivos_R3)



#-----------
library(readxl)
library(dplyr)
library(VennDiagram)

# Carregar os dados do arquivo .xlsx
dados <- read_excel("/Users/LS28_Seq/Documents/JM_03.10.24/artigo2/tabelas/ASV_table_bruto.xlsx")

# Função para encontrar ASVs comuns entre um vetor de amostras
encontrar_asvs_comuns <- function(amostras) {
  dados %>%
    group_by(ASV) %>%
    summarise(presente = sum(`Abundância relativa (%)`[Amostra %in% amostras] > 0, na.rm = TRUE)) %>%
    filter(presente == length(amostras))
}

# Função para encontrar ASVs exclusivos para uma amostra
encontrar_asvs_exclusivos <- function(amostra, amostras_excluidas) {
  dados %>%
    group_by(ASV) %>%
    summarise(
      presente_v = sum(`Abundância relativa (%)`[Amostra == amostra] > 0, na.rm = TRUE),
      presente_outros = sum(`Abundância relativa (%)`[Amostra %in% amostras_excluidas] > 0, na.rm = TRUE)
    ) %>%
    filter(presente_v > 0 & presente_outros == 0)
}

# Encontrar ASVs comuns entre pares e grupos
asvs_comuns_2 <- list(
  JM3I_JM3R1 = encontrar_asvs_comuns(c("JM3I", "JM3R1")),
  JM3I_JM3R2 = encontrar_asvs_comuns(c("JM3I", "JM3R2")),
  JM3I_JM3R3 = encontrar_asvs_comuns(c("JM3I", "JM3R3")),
  JM3R1_JM3R2 = encontrar_asvs_comuns(c("JM3R1", "JM3R2")),
  JM3R1_JM3R3 = encontrar_asvs_comuns(c("JM3R1", "JM3R3")),
  JM3R2_JM3R3 = encontrar_asvs_comuns(c("JM3R2", "JM3R3"))
)

# Encontrar ASVs comuns entre todos
asvs_comuns_4 <- encontrar_asvs_comuns(c("JM3I", "JM3R1", "JM3R2", "JM3R3"))

# Filtrar ASVs exclusivos para cada amostra
asvs_exclusivos_JM3I <- encontrar_asvs_exclusivos("JM3I", c("JM3R1", "JM3R2", "JM3R3"))
asvs_exclusivos_JM3R1 <- encontrar_asvs_exclusivos("JM3R1", c("JM3I", "JM3R2", "JM3R3"))
asvs_exclusivos_JM3R2 <- encontrar_asvs_exclusivos("JM3R2", c("JM3I", "JM3R1", "JM3R3"))
asvs_exclusivos_JM3R3 <- encontrar_asvs_exclusivos("JM3R3", c("JM3I", "JM3R1", "JM3R2"))

# Criar uma lista para o gráfico de Venn
venn_list <- list(
  JM3I = unique(c(asvs_exclusivos_JM3I$ASV, asvs_comuns_2$JM3I_JM3R1$ASV, asvs_comuns_2$JM3I_JM3R2$ASV, asvs_comuns_2$JM3I_JM3R3$ASV, asvs_comuns_4$ASV)),
  JM3R1 = unique(c(asvs_exclusivos_JM3R1$ASV, asvs_comuns_2$JM3I_JM3R1$ASV, asvs_comuns_2$JM3R1_JM3R2$ASV, asvs_comuns_2$JM3R1_JM3R3$ASV, asvs_comuns_4$ASV)),
  JM3R2 = unique(c(asvs_exclusivos_JM3R2$ASV, asvs_comuns_2$JM3I_JM3R2$ASV, asvs_comuns_2$JM3R1_JM3R2$ASV, asvs_comuns_2$JM3R2_JM3R3$ASV, asvs_comuns_4$ASV)),
  JM3R3 = unique(c(asvs_exclusivos_JM3R3$ASV, asvs_comuns_2$JM3I_JM3R3$ASV, asvs_comuns_2$JM3R1_JM3R3$ASV, asvs_comuns_2$JM3R2_JM3R3$ASV, asvs_comuns_4$ASV))
)

# Plotar o gráfico de Venn (420x335)
venn.diagram(
  x = venn_list,
  category.names = c("I", "R1", "R2", "R3"),
  filename = 'figuras/venn.tiff',
  imagetype="tiff" ,
  fill = c("#393b79ff", "#BC3C29", "#637939FF", "#8C6D31FF"),
  fontfamily = "serif",
  cat.fontfamily = "serif",
  lwd = 0.25,
  resolution = 950,
  col = "black",
  cat.cex = 0.3,  # Tamanho das letras das categorias
  cex = 0.3,      # Tamanho das letras dos números
  height = 550,     # Altura da figura em px
  width = 500       # Largura da figura em px
)


