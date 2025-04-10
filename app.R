library(shiny)
library(tidyverse)
library(phyloseq)
library(reshape2)
library(ggsci)
library(ggpubr)
library(viridis)
library(RColorBrewer)

# Carrega dados
ps <- readRDS("ps_silva.rds")
samdf <- read.table("metadados.txt", header = TRUE, encoding = "UTF-8", sep = "\t")

# UI
ui <- fluidPage(
  titlePanel("Visualização Interativa de Metagenômica"),

  sidebarLayout(
    sidebarPanel(
      # Seletor para escolher o gráfico
      selectInput("grafico", "Escolha o Gráfico:",
                  choices = c("Gráfico de Barras", "Heatmap"),
                  selected = "Gráfico de Barras"),

      selectInput("taxrank", "Nível Taxonômico:",
                  choices = c("phylum", "class", "order", "family", "genus", "species"),
                  selected = "phylum"),
      sliderInput("abund_min", "Abundância Mínima (%):", min = 0, max = 20, value = 5, step = 1),

      selectInput("paleta", "Paleta para gráfico de barras:",
                  choices = c("ggsci::npg", "ggsci::nejm", "ggsci::lancet", "viridis", "brewer"),
                  selected = "ggsci::npg"),

      selectInput("heat_color", "Paleta para heatmap:",
                  choices = c("Verde-Azul" = "green", "Roxo-Rosa" = "magma", "Viridis" = "viridis", "Plasma" = "plasma"),
                  selected = "green"),

      downloadButton("downloadBar", "Baixar gráfico de barras (.tiff)"),
      downloadButton("downloadHeat", "Baixar heatmap (.tiff)")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Gráfico de Barras", plotOutput("barPlot", height = "800px")),
        tabPanel("Heatmap", plotOutput("heatPlot", height = "800px"))
      )
    )
  )
)

# Funções auxiliares
get_palette <- function(paleta, n) {
  if (paleta == "ggsci::npg") return(ggsci::pal_npg("npg")(n))
  if (paleta == "ggsci::nejm") return(ggsci::pal_nejm("default")(n))
  if (paleta == "ggsci::lancet") return(ggsci::pal_lancet("lanonc")(n))
  if (paleta == "viridis") return(viridis::viridis(n))
  if (paleta == "brewer") return(RColorBrewer::brewer.pal(min(n, 8), "Set2"))
  return(hcl.colors(n, "Dark 3"))
}

get_heatmap_colors <- function(option) {
  switch(option,
         "green" = scale_fill_gradient(high = "#44b19f", low = "gray98",
                                       guide = guide_colorbar(barwidth = 1, barheight = 5)),
         "magma" = scale_fill_viridis_c(option = "magma", guide = guide_colorbar(barwidth = 1, barheight = 5)),
         "viridis" = scale_fill_viridis_c(option = "viridis", guide = guide_colorbar(barwidth = 1, barheight = 5)),
         "plasma" = scale_fill_viridis_c(option = "plasma", guide = guide_colorbar(barwidth = 1, barheight = 5))
  )
}

# Server
server <- function(input, output) {

  dados_processados <- reactive({
    order_samples <- paste0(samdf$sampleName)
    seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
    rownames(samdf) <- samdf$sampleName
    sample_data(seqtab) <- as.data.frame(samdf)

    taxrank <- input$taxrank
    abund_min <- input$abund_min

    rank_seqtab <- tax_glom(seqtab, taxrank = taxrank)
    rank_counts_tab <- otu_table(rank_seqtab)
    rank_tax_vec <- as.vector(tax_table(rank_seqtab)[, taxrank])
    colnames(rank_counts_tab) <- as.vector(rank_tax_vec)

    rank_counts_tab2 <- otu_table(tax_glom(seqtab, taxrank = taxrank, NArm = FALSE))
    rank_tax_vec2 <- as.vector(tax_table(tax_glom(seqtab, taxrank = taxrank, NArm = FALSE))[, taxrank])
    colnames(rank_counts_tab2) <- as.vector(rank_tax_vec2)

    unclassified_tax_counts <- rowSums(rank_counts_tab2) - rowSums(rank_counts_tab)
    rank_and_unidentified_counts_tab <- cbind(rank_counts_tab, "Unclassified" = unclassified_tax_counts)

    asv_filt <- rank_and_unidentified_counts_tab[
      , colSums(rank_and_unidentified_counts_tab) * 100 / sum(rank_and_unidentified_counts_tab) >= abund_min
    ]

    if (abund_min > 0) {
      asv_filt <- cbind(
        asv_filt,
        data.frame(Others = rowSums(
          rank_and_unidentified_counts_tab[
            , !(colnames(rank_and_unidentified_counts_tab) %in% colnames(asv_filt))
          ]
        ))
      )
    }

    rank_taxa_proportions_tab <- apply(asv_filt, 1, function(x) x / sum(x) * 100)
    rank_taxa_proportions_tab_non0 <- rank_taxa_proportions_tab[
      rowSums(rank_taxa_proportions_tab) > 0,
    ]

    tidy_rank_taxa_proportions_tab <- reshape2::melt(rank_taxa_proportions_tab_non0,
                                                     value.name = "abundance",
                                                     varnames = c("taxa", "sampleName")) %>%
      merge(data.frame(sample_data(seqtab)), by = "sampleName")

    tidy_rank_taxa_proportions_tab$Fase <- factor(tidy_rank_taxa_proportions_tab$Fase,
                                                  levels = c("Initial", "Final"))

    tidy_rank_taxa_proportions_tab$sampleName2 <- factor(
      tidy_rank_taxa_proportions_tab$sampleName2,
      levels = c("End. CT", "NC CT", "Blank",
                 "EXP-BC-10", "EXP-BC-20", "EXP-BC-30",
                 "EXP-CAG-10", "EXP-CAG-20", "EXP-CAG-30")
    )

    list(tidy = tidy_rank_taxa_proportions_tab, taxrank = taxrank)
  })

  output$barPlot <- renderPlot({
    if(input$grafico == "Gráfico de Barras") {
      df <- dados_processados()
      tidy <- df$tidy
      taxrank <- df$taxrank

      plot_labels <- sort(as.character(unique(tidy$taxa)))
      if ("Others" %in% tidy$taxa) plot_labels <- c(setdiff(plot_labels, "Others"), "Others")
      if ("Unclassified" %in% tidy$taxa) plot_labels <- c(setdiff(plot_labels, "Unclassified"), "Unclassified")

      tidy$taxa <- factor(tidy$taxa, levels = plot_labels)

      cores_taxa <- get_palette(input$paleta, length(levels(tidy$taxa)))

      ggplot(tidy) +
        geom_col(aes(x = sampleName2, y = abundance, fill = taxa), position = "stack") +
        facet_grid(~Fase, scales = "free", space = "fixed") +
        theme_bw(base_size = 12) +
        labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
        theme(
          aspect.ratio = 1.5,
          panel.spacing = unit(-.1, "lines"),
          text = element_text(family = "serif", color = "black", face = "plain"),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle = 30, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 13, hjust = 1),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 12),
          legend.key.size = unit(0.7, "lines")
        ) +
        scale_fill_manual(values = cores_taxa, guide = guide_legend(ncol = 1)) +
        scale_y_continuous(expand = expansion(mult = c(.01, .01)), labels = function(x) paste0(x, "%")) +
        scale_x_discrete(expand = expansion(mult = c(0.08, 0.08)))
    }
  })

  output$heatPlot <- renderPlot({
    if(input$grafico == "Heatmap") {
      df <- dados_processados()
      tidy <- df$tidy
      taxrank <- df$taxrank

      tidy <- tidy[!tidy$taxa %in% c("Unclassified", "Others"), ]
      tidy$taxa <- factor(tidy$taxa, levels = rev(sort(unique(tidy$taxa))))

      ggplot(tidy) +
        geom_tile(aes(x = sampleName2, y = taxa, fill = abundance), color = "gray98") +
        get_heatmap_colors(input$heat_color) +
        theme_bw(base_size = 12) +
        labs(y = str_to_title(taxrank), fill = "Abundance (%)") +
        theme(
          aspect.ratio = 1.5,
          text = element_text(family = "serif", color = "black", face = "plain"),
          strip.placement = "outside",
          panel.spacing = unit(.3, "lines"),
          strip.background = element_blank(),
          strip.text = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11),
          axis.text.y = element_text(size = 11, hjust = 1),
          axis.text.x = element_text(size = 12, angle = 30, hjust = 1, color = "black"),
          axis.title.x = element_blank(),
          panel.grid.minor = element_line(color = "white", size = 0.1),
          panel.border = element_rect(color = "gray40", size = 0.8),
          panel.grid.major = element_blank()
        ) +
        scale_y_discrete(position = "left", expand = expansion(mult = c(0, 0))) +
        scale_x_discrete(expand = expansion(mult = c(0, 0))) +
        facet_grid(~Fase, scales = "free", space = "fixed")
    }
  })

  output$downloadBar <- downloadHandler(
    filename = function() {"grafico_barras.tiff"},
    content = function(file) {
      ggsave(file, plot = output$barPlot(), width = 12, height = 8, dpi = 300, units = "in", device = "tiff")
    }
  )

  output$downloadHeat <- downloadHandler(
    filename = function() {"heatmap.tiff"},
    content = function(file) {
      ggsave(file, plot = output$heatPlot(), width = 12, height = 8, dpi = 300, units = "in", device = "tiff")
    }
  )
}

shinyApp(ui = ui, server = server)
