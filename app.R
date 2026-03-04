
# -------------------------
# Pakete laden
# -------------------------
library(shiny)
library(golubEsets)
library(RColorBrewer)
library(pheatmap)
library(heatmaply)

# -------------------------
# Golub-Daten laden
# -------------------------
data(Golub_Train)
x <- exprs(Golub_Train)

# Spaltennamen ALL/AML anhängen
colnames(x) <- paste(
  pData(Golub_Train)$Samples,
  pData(Golub_Train)$ALL.AML,
  sep = "_"
)

# Werte < 1 ersetzen (wegen log2)
x[x < 1] <- 1

# Log2
x_log <- log2(x)

# -------------------------
# Patientenzahlen berechnen
# -------------------------
num_ALL <- sum(pData(Golub_Train)$ALL.AML == "ALL")
num_AML <- sum(pData(Golub_Train)$ALL.AML == "AML")

# -------------------------
# UI
# -------------------------
ui <- fluidPage(

  titlePanel("Golub Heatmap – pheatmap & heatmaply"),

  sidebarLayout(

    sidebarPanel(

      h3("Patientenübersicht"),
      tags$div(
        style="background:#e8f1ff;padding:8px;border-radius:5px;margin-bottom:6px;",
        strong("ALL: "), num_ALL
      ),
      tags$div(
        style="background:#ffe8e8;padding:8px;border-radius:5px;margin-bottom:12px;",
        strong("AML: "), num_AML
      ),

      hr(),

      sliderInput("numberOfGenes",
                  "Anzahl Gene (höchste Varianz)",
                  min = 10, max = 200, value = 50),

      selectInput("distMea",
                  "Distanzmaß",
                  choices = c("euclidean", "maximum",
                              "manhattan", "canberra",
                              "binary", "minkowski")),

      selectInput("clustMeth",
                  "Clustering-Methode",
                  choices = c("ward.D", "ward.D2",
                              "single", "complete",
                              "average", "mcquitty",
                              "median", "centroid"))
    ),

    mainPanel(
      h3("Statische Heatmap (pheatmap)"),
      plotOutput("pheat", height = 600),

      hr(),

      h3("Interaktive Heatmap (heatmaply)"),
      heatmaplyOutput("interactive")
    )
  )
)

# -------------------------
# SERVER
# -------------------------
server <- function(input, output) {

  # Reaktives Objekt: Gene mit höchster Varianz
  selectedGenes <- reactive({
    vars <- apply(x_log, 1, var)
    names(sort(vars, decreasing = TRUE)[1:input$numberOfGenes])
  })

  # Reaktive Matrix (Z-Score Normalisierung)
  selectedMatrix <- reactive({
    m <- x_log[selectedGenes(), ]
    t(scale(t(m)))   # Z-Score pro Gen
  })

  # Annotation: ALL / AML
  annot <- reactive({
    data.frame(
      Typ = pData(Golub_Train)$ALL.AML,
      row.names = colnames(x_log)
    )
  })

  # -------------------------
  # pheatmap (statisch)
  # -------------------------
  output$pheat <- renderPlot({

    pheatmap(
      selectedMatrix(),
      annotation_col = annot(),
      clustering_distance_rows = input$distMea,
      clustering_distance_cols = input$distMea,
      clustering_method = input$clustMeth,
      scale = "none",
      color = colorRampPalette(brewer.pal(9, "Blues"))(200),
      show_rownames = FALSE,
      fontsize_col = 8
    )

  })

  # -------------------------
  # Interaktive Heatmap
  # -------------------------
  output$interactive <- renderHeatmaply({

    heatmaply(
      selectedMatrix(),
      colors = colorRampPalette(brewer.pal(9, "Blues"))(200),
      k_col = 2,
      k_row = 2,
      Rowv = TRUE,
      Colv = TRUE,
      xlab = "Gene",
      ylab = "Patienten",
      main = "Interaktive Heatmap"
    )

  })
}

# -------------------------
# App starten
# -------------------------
shinyApp(ui, server)
