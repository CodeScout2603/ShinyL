
# benötigte Pakete
library(shiny)
library(golubEsets)
library(RColorBrewer)
library(pheatmap)

# -------------------------
# Daten laden und vorbereiten
# -------------------------
data(Golub_Train)

x <- exprs(Golub_Train)

# Sample-Namen schöner gestalten
colnames(x) <- paste(
  pData(Golub_Train)$Samples, 
  pData(Golub_Train)$ALL.AML, 
  sep = "_"
)

# Werte < 1 auf 1 setzen und log2 transformieren
xLog <- log2(pmax(x, 1))

# Patientenzahlen bestimmen
num_ALL <- sum(pData(Golub_Train)$ALL.AML == "ALL")
num_AML <- sum(pData(Golub_Train)$ALL.AML == "AML")

# Varianz der Gene einmalig berechnen (Performance!)
geneVariance <- apply(xLog, 1, var)
sortedGenes <- names(sort(geneVariance, decreasing = TRUE))

# Annotation für Heatmap
annotation <- data.frame(
  Type = pData(Golub_Train)$ALL.AML
)
rownames(annotation) <- colnames(xLog)

# -------------------------
# UI Definition
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
      plotOutput("heatmap", height = reactive({ input$numberOfGenes * 18 }))
    )
  )
)

# -------------------------
# Server Logik
# -------------------------
server <- function(input, output) {
  
  # reaktive Auswahl der variabelsten Gene
  selectedGenes <- reactive({
    xLog[sortedGenes[1:input$numberOfGenes], ]
  })
  
  output$heatmap <- renderPlot({
    
    tx <- t(selectedGenes())
    
    # Validierung: binary Distanz benötigt binäre Matrix
    validate(
      need(!(input$distMea == "binary" && !all(tx %in% c(0, 1))),
           "Binary distance erfordert binäre Daten (0/1).")
    )
    
    # Heatmap zeichnen
    pheatmap(
      tx,
      clustering_distance_rows = input$distMea,
      clustering_distance_cols = input$distMea,
      clustering_method = input$clustMeth,
      annotation_row = annotation,
      color = colorRampPalette(brewer.pal(8, "Blues"))(25),
      fontsize = 9,
      main = "Heatmap der Top-variablen Gene"
    )
  })
}

# Shiny App starten
shinyApp(ui, server)
