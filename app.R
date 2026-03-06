layout(matrix(c(1,2), nrow=2), heights = c(8,1))
par(mar = c(2, 8, 4, 2))
library(shiny)
library(golubEsets)
library(RColorBrewer)
library(DT)

# ---------------------
# Datenvorbereitung
# ---------------------

data(Golub_Train)

x <- exprs(Golub_Train)
colnames(x) <- paste(pData(Golub_Train)$Samples, pData(Golub_Train)$ALL.AML, sep = "_")

x[x < 1] <- 1
xLogarithmised <- log2(x)


cleanGenes <- sub("_.*$", "", sort(rownames(x)))

genes <- data.frame(
  Gene = paste0(
    '<a href="https://www.ncbi.nlm.nih.gov/gene/?term=',
    cleanGenes,
    '" target="_blank">',
    cleanGenes,
    '</a>'
  )
)


num_ALL <- sum(pData(Golub_Train)$ALL.AML == "ALL")
num_AML <- sum(pData(Golub_Train)$ALL.AML == "AML")


# ---------------------
# Benutzeroberfläche
# ---------------------

ui <- fluidPage(

  # kleines globales CSS
  tags$head(tags$style(HTML("
    body { font-family: Arial; }
    .sidebarPanel { font-size: 14px; }

    .heatmap-desc-box {
      background: #f8f9fa;
      border: 1px solid #d3d3d3;
      border-radius: 8px;
      padding: 12px 16px;
      margin: 10px 0 20px 0;
      font-size: 14px;
      color: #333;
    }
  "))),

  titlePanel("Heatmap of Patients and Genes"),

  sidebarLayout(
    sidebarPanel(
      h4("📊 Anzahl Patienten"),
      p(paste("ALL:", num_ALL)),
      p(paste("AML:", num_AML)),
      br(),

      h4("⚙️ Einstellungen"),
      sliderInput("numberOfGenes",
                  "Number of Genes",
                  min = 10, max = 100, value = 50),

      selectInput("distMea", "Distance Measure",
                  choices = c("euclidean", "maximum", "manhattan",
                              "canberra", "binary", "minkowski")),

      
selectInput("clustMeth", "Clustering Method",
            choices = c("ward.D", "ward.D2", "single", "complete",
                        "average", "mcquitty", "median", "centroid")),

# >>> Punkt B (hier einfügen) <<<
        h4("📝 Heatmap-Beschreibung"),
        checkboxInput("showDesc", "Beschreibung anzeigen", value = TRUE),
        radioButtons(
        "descPos", "Position der Beschreibung",
        choices = c("Über der Heatmap" = "above", "Unter der Heatmap" = "below"),
        inline = TRUE
    ),
# <<< Ende Punkt B >>>

       
       
       br(),
      

      h4("🔗 GitHub"),
      tags$a(href = "https://github.com/CodeScout2603/ShinyL-",
             target = "_blank", "Mein Code auf GitHub"),
      br(), br(),

      h4("🧬 Genliste"),
      div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 6px;",
          DTOutput("geneTable"))
    ),

    mainPanel(
      uiOutput("descAbove"),
      plotOutput("heatmap", height = 900),
      uiOutput("descBelow")

    )
  )
)


# ---------------------
# Server-Logik
# ---------------------

server <- function(input, output) {

  # reaktive Auswahl der Top-Gene (performanter!)
  selectedGenes <- reactive({
    sorted <- sort(apply(xLogarithmised, 1, var), decreasing = TRUE)
    topGenes <- names(sorted)[1:input$numberOfGenes]
    xLogarithmised[topGenes, ]
  })

    


desc_ui <- reactive({
  req(input$showDesc)
  HTML(sprintf(
    '
    <div class="heatmap-desc-box">
      <b>Heatmap-Beschreibung</b><br>
      Diese Heatmap zeigt die %d Gene mit der höchsten Varianz.<br>
      Abstand: <b>%s</b>, Clustering: <b>%s</b>.
    </div>
    ',
    input$numberOfGenes, input$distMea, input$clustMeth
  ))
})

output$descAbove <- renderUI({
    if (isTRUE(input$showDesc) && input$descPos == "above") desc_ui()
})

    output$descBelow <- renderUI({
    if (isTRUE(input$showDesc) && input$descPos == "below") desc_ui()
})
  
 
output$heatmap <- renderPlot({
  tx <- t(selectedGenes())

  # Layout: Panel 1 = Heatmap, Panel 2 = X-Achsen-Titel
  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(8, 1))

  # Panel 1: Heatmap + linker Rand für Y-Titel
  par(mar = c(2, 8, 4, 2))
  heatmap(
    tx,
    distfun = function(c) dist(c, method = input$distMea),
    hclustfun = function(c) hclust(c, method = input$clustMeth),
    col = colorRampPalette(brewer.pal(8, "Blues"))(25),
    labRow = FALSE,   # << verhindert doppelte/überlappende Zeilenlabels
    labCol = FALSE    # << verhindert doppelte/überlappende Spaltenlabels
  )

  # (1) Genau EIN Y-Achsentitel (links, nah an der Heatmap)
  mtext("Gene (Top-Varianz-Auswahl)", side = 2, line = 4, cex = 1.1)

  # Panel 2: Nur für den X-Achsentitel (sauber unter der Heatmap)
  par(mar = c(5, 2, 1, 2))
  plot.new()

  # (2) Genau EIN X-Achsentitel (unten, zentriert)
  mtext("Patienten (Samples mit ALL/AML-Diagnose)", side = 1, line = 2, cex = 1.1)
  })

output$geneTable <- renderDT({
    datatable(
      genes,
      escape = FALSE,   # << NICHT VERGESSEN! HTML erlauben
      options = list(
        pageLength = 8,
        autoWidth = TRUE,
        searching = TRUE,
        ordering = TRUE
      ),
      rownames = FALSE
    )
  })
}

# ---------------------
# App starten
# ---------------------
shinyApp(ui, server)