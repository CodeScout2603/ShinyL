library(golubEsets)
# collect the data
data(Golub_Train)
# get the expression data
x = exprs(Golub_Train)
# indicate for each patient ALL or AML
colnames(x) <- paste(pData(Golub_Train)$Samples, pData(Golub_Train)$ALL.AML, sep="_")
# set all values to at least 1 to avoid NaNs
xWithoutLT1 = replace(x, x<1,1)
# logarithmize x 
xLogarithmised = log2(xWithoutLT1)
# sort the genenames for the table on the more info tab
genes = data.frame(sort(rownames(x)))

# Anzahl ALL und AML bestimmen
num_ALL <- sum(pData(Golub_Train)$ALL.AML == "ALL")
num_AML <- sum(pData(Golub_Train)$ALL.AML == "AML")
# ------------------------

# RColorBrewer for better color of the heatmap
library("RColorBrewer")
# DT for inserting a table in tab: More Info
library(DT)

# user interface object
ui <- fluidPage(
  
  titlePanel("Heatmap of Patients and Genes"),
  
  sidebarLayout(
    sidebarPanel(
     
	 
 # >>> HIER: https://github.dev/CodeScout2603/ShinyL-/blob/main/app.R
      tags$a(
        href = "https://github.com/CodeScout2603/ShinyL-",
        target = "_blank",
        "🔗 Mein Code auf GitHub"
      ),
      br(), br(),

	 h4("Anzahl Patienten"),
      p(paste("ALL:", num_ALL)),
      p(paste("AML:", num_AML)),
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
      plotOutput("heatmap", height = 900)
    )
  )
)

# server logic unit
server <- function(input, output) {
  # rendering the heatmap plot
  output$heatmap <- renderPlot({
    # first the user chosen number of genes with the highest expression are selected
    xHighestEX = xLogarithmised[names(sort(apply(xLogarithmised,1,var), decreasing=TRUE)[1:input$numberOfGenes]),]
    # this matrix has now to be tranposed for better understanding of the heatmap
    tx = t(xHighestEX)
    # the heatmap is printet with the matrix the chosen dist and clust function and the blue color of RColorBrewer
    heatmap(tx,distfun=function(c){dist(c,method=input$distMea)}, hclustfun=function(c){hclust(c,method=input$clustMeth)}, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
  })
  #Render table
  output$geneTable <- renderDT(genes)
}

# Generate the app
shinyApp(ui, server)