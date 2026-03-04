library(shiny)
library(golubEsets)
library(RColorBrewer)
library(pheatmap)

# -------------------------
# Daten laden und vorbereiten
# -------------------------
data(Golub_Train)

x <- exprs(Golub_Train)

# Sample-Namen erstellen
sample_labels <- paste(
  pData(Golub_Train)$Samples,
  pData(Golub_Train)$ALL.AML,
  sep = "_"
)
colnames(x) <- sample_labels

# Log2-Transformation
xLog <- log2(pmax(x, 1))

# Patientenzahlen
num_ALL <- sum(pData(Golub_Train)$ALL.AML == "ALL")
num_AML <- sum(pData(Golub_Train)$ALL.AML == "AML")

# Varianzberechnung
geneVariance <- apply(xLog, 1, var)
sortedGenes <- names(sort(geneVariance, decreasing = TRUE))

# Annotation (für Spalten!)
annotation <- data.frame(
  Leukemia = pData(Golub_Train)$ALL.AML
)
rownames(annotation) <- colnames(xLog)

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  
  titlePanel("Heatmap of Patients and Genes"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Anzahl Patienten"),
      tags$ul(
        tags$li(paste("ALL:", num_ALL)),
        tags$li(paste("AML:", num_AML))
      ),
      
      sliderInput("numberOfGenes",
                  "Number of Genes",
                  min = 10,
                  max = 100,
                  value = 50),
      
      selectInput("distMea",
                  "Distance Measure",
                  choices = c("euclidean", "maximum",
                              "manhattan", "canberra",
                              "binary", "minkowski")),
      
      selectInput("clustMeth",
                  "Clustering Method",
                  choices = c("ward.D", "ward.D2",
                              "single", "complete",
                              "average", "mcquitty",
                              "median", "centroid"))
    ),
    
    mainPanel(
      plotOutput("heatmap", height = "900px")
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output) {
  
  selectedGenes <- reactive({
    xLog[sortedGenes[1:input$numberOfGenes], , drop = FALSE]
  })
  
  output$heatmap <- renderPlot({
    
    mat <- selectedGenes()
    
    pheatmap(
      mat,   # << NICHT transponiert!
      clustering_distance_rows = input$distMea,
      clustering_distance_cols = input$distMea,
      clustering_method = input$clustMeth,
      annotation_col = annotation,
      color = colorRampPalette(brewer.pal(8, "Blues"))(25),
      fontsize = 9,
      main = "Heatmap der Top-variablen Gene"
    )
  })
}

shinyApp(ui, server)