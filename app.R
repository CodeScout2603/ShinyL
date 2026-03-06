
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

# Werte < 1 auf 1 setzen, dann log2
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
    .heatmap-desc { 
      background: #f8f9fb; 
      border: 1px solid #e5e7eb; 
      border-radius: 6px; 
      padding: 12px; 
      margin: 10px 0 18px 0;
    }
    .heatmap-desc h4 { margin-top: 0; }
    .muted { color: #6b7280; font-size: 12px; }
  "))),

  titlePanel("Heatmap of Patients and Genes"),

  sidebarLayout(
    sidebarPanel(
      h4("📊 Anzahl Patienten"),
      p(paste("ALL:", num_ALL)),
      p(paste("AML:", num_AML)),
      br(),

      h4("⚙️ Einstellungen"),
      sliderInput(
        "numberOfGenes",
        "Number of Genes",
        min = 10, max = 100, value = 50
      ),

      selectInput(
        "distMea", "Distance Measure",
        choices = c("euclidean", "maximum", "manhattan",
                    "canberra", "binary", "minkowski"),
        selected = "euclidean"
      ),

      selectInput(
        "clustMeth", "Clustering Method",
        choices = c("ward.D", "ward.D2", "single", "complete",
                    "average", "mcquitty", "median", "centroid"),
        selected = "complete"
      ),
      br(),

      h4("📝 Heatmap-Beschreibung"),
      checkboxInput("showDesc", "Beschreibung anzeigen", value = TRUE),
      radioButtons(
        "descPos", "Position der Beschreibung",
        choices = c("Über der Heatmap" = "above", "Unter der Heatmap" = "below"),
        selected = "above", inline = TRUE
      ),
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
      # Beschreibung oben (konditional)
      uiOutput("descAbove"),
      plotOutput("heatmap", height = 900),
      # Beschreibung unten (konditional)
      uiOutput("descBelow")
    )
  )
)

# ---------------------
# Server-Logik
# ---------------------

server <- function(input, output) {

  # reaktive Auswahl der Top-Gene
  selectedGenes <- reactive({
    sorted <- sort(apply(xLogarithmised, 1, var), decreasing = TRUE)
    topGenes <- names(sorted)[1:input$numberOfGenes]
    xLogarithmised[topGenes, ]
  })

  # Dynamische Beschreibung als HTML
  desc_ui <- reactive({
    req(input$showDesc)

    HTML(sprintf(
      '
      <div class="heatmap-desc">
        <h4>🗺️ Was zeigt diese Heatmap?</h4>
        <p>
          Die Heatmap visualisiert die <b>Expressionswerte</b> der %d Gene mit der
          höchsten Varianz (log2-transformiert; Werte &lt; 1 wurden vorab auf 1 gesetzt),
          gemessen in den <b>Golub</b>-Trainingsdaten (ALL/AML).<br>
          Jede <b>Spalte</b> entspricht einem Patienten (Spaltennamen enthalten die Diagnose),
          jede <b>Zeile</b> einem Gen. Die Farbskala kodiert relative Expressionsniveaus
          (hell = niedriger, dunkel = höher).
        </p>
        <p>
          Die hierarchische Clusteranalyse ordnet <b>ähnliche Patienten</b> und
          <b>ähnliche Gene</b> zu Gruppen. Nähe/Distanz wird mit
          <code>%s</code> berechnet; Dendrogramme entstehen mit
          <code>%s</code>. Cluster können Hinweise auf biologische Muster
          und eine Trennung zwischen <b>ALL</b> und <b>AML</b> geben.
        </p>
        <p class="muted">
          Hinweis: Die Auswahl der Top-Gene nach Varianz erhöht den Kontrast der Muster,
          ist aber keine Feature-Selection für Vorhersagen.
        </p>
      </div>
      ',
      input$numberOfGenes, input$distMea, input$clustMeth
    ))
  })

  output$descAbove <- renderUI({
    if (isTRUE(input$showDesc) && identical(input$descPos, "above")) desc_ui()
  })

  output$descBelow <- renderUI({
    if (isTRUE(input$showDesc) && identical(input$descPos, "below")) desc_ui()
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
      escape = FALSE,   # HTML erlauben
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
