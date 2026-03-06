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
    "
    <div>
      <h4>Heatmap Erklärung</h4>
      Diese Heatmap zeigt die %d Gene mit der höchsten Varianz.
      Farbwerte entsprechen log2-transformierten Expressionswerten.
      Abstand: <b>%s</b>, Clustering-Methode: <b>%s</b>.
    </div>
    ",
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

    heatmap(
      tx,
      distfun = function(c) dist(c, method = input$distMea),
      hclustfun = function(c) hclust(c, method = input$clustMeth),
      col = colorRampPalette(brewer.pal(8, "Blues"))(25)
    )
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