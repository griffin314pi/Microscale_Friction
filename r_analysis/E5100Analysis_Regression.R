# Include Libraries

library(shiny)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(readr)
library(gridExtra)
library(plotly)


# Define constants and parameters
z.quartz <- 8800000
area.effective <- 0.0025^2*pi
amplitude.factor <- 1.4 #pm/V
d26 <- 3.1*0.000000000001

U.star = 1
F.elas.star = 1
F.k = 1
K = 1
G.star = 14
a = K/(8*G.star)

# Define UI
ui <- fluidPage(
  titlePanel(
    h1("E5100 Data Analysis", align = "left")
  ),
  # Page will be in row format
  fluidRow(column(width=4,offset=0,
                  div(style = "height:50px;", "by: Eskil Irgens"))),
  # File input
  fluidRow(column(width=4,
                  fileInput("upload", NULL, buttonLabel = "Choose File:", multiple = FALSE))),
  # Zoom function
  fluidRow(
    column(width = 6, offset = 3,
           sliderInput(width = 600, inputId = "ZoomSlider", label = h3("Zoom in on graph"), 
                       min = 0, max = 100, value = c(0, 100)))
  ),
  # Row 1 - Frequency and Bandwith
  fluidRow(column(width=6,plotOutput("Frequency")),
           column(width=6,plotOutput("Bandwith"))),
  fluidRow(
    column(width = 6, offset = 3,
           sliderInput(width = 600, inputId = "SelectionSlider", label = h3("Area selected"),
                       min = 0, max = 100, value = c(0,100)))
  ),
  fluidRow(
    column(width = 6,plotOutput("RegressionPlotFrequency")),
    column(width = 6,plotOutput("RegressionPlotBandwith"))
  ),
  fluidRow(column(width = 6,tableOutput("RegressionResultsFrequency")),
           column(width = 6, tableOutput("RegressionResultsBandwith")))
)

# Server Logic to make plots
server <- function(input, output) {
  # Code for file input
  test.tbl <- reactive({
    req(input$upload)
    inFile <- input$upload
    tbl <- read.delim(inFile$datapath, header=FALSE)
    return(tbl)
  })
  
  # Extract segments when analysis button pressed
  
  
  # Extract the testing conditions and units of measurement from input data
  testing.conditions <- reactive ({ test.tbl()[1,c(2:11)] })
  units.of.measurement <- reactive ({ test.tbl()[c(2,3),c(3:20)] })
  # print(testing.conditions)
  # print(units.of.measurement)
  normal.load <- reactive ({ as.numeric(gsub(" uN", "", testing.conditions()[7])) })
  output$Load <- reactive ({ paste("Normal Load:",normal.load(),"uN") })
  
  # Code for tidying the input dataset
  
  # Remove empty columns and remove testing condition rows
  tidy.tbl <- reactive({ new.tbl <- test.tbl()[-c(1,3),c(3:20)]
  # Set column names to first row
  colnames(new.tbl) <- new.tbl[1,] 
  # Remove first row (which is now the column names)
  new.tbl[-1,]
  })
  
  # Add gamma and RMS delta V columns to full dataset
  
  # remove and parse time column for later
  Time.vec <- reactive({ parse_time(tidy.tbl()$Time) })
  segType <- reactive ({ tidy.tbl()$Seg })
  segID <- reactive ({ tidy.tbl()$`ID Tag` })
  
  tidy.tbl2 <- reactive ({ tidy.tbl()[,c(1,3,5:18)]%>%
      # Change data type of all columns from character to numeric
      mutate_if(is.character,as.numeric)%>%
      # Compute gamma
      mutate(Gamma = R1/(4*pi*L1))%>%
      mutate(RMS.deltaV = `Ch1 (RMS)` - `Ch2 (RMS)`)%>%
      mutate(Time = Time.vec())%>%
      mutate(segType=segType())%>%
      mutate(segID=segID())%>%
      mutate(Time0 = as.numeric(Time) - as.numeric(Time[1]))
  })
  # Creating range for zoom
  observeEvent(!is.null(tidy.tbl2()),{
    maxpoint <- rev(tidy.tbl2()$Time0)[1]%>%as.numeric
    updateSliderInput(inputId = "ZoomSlider", max = maxpoint,value = c(0,maxpoint))
  })
  # Extracting zoomed in data
  zoom.data <- reactive({ #reactive object updates every time the user changes the slider
    filter(tidy.tbl2(), Time0>input$ZoomSlider[1] & Time0<input$ZoomSlider[2])
  })
  observeEvent(input$ZoomSlider, {
    updateSliderInput(inputId = "SelectionSlider", min = input$ZoomSlider[1], max = input$ZoomSlider[2])
  })
  regression.data <- reactive({ #reactive object updates every time the user changes the slider
    filter(tidy.tbl2(), Time0>input$SelectionSlider[1] & Time0<input$SelectionSlider[2])
  })
  # Graph Frequency data
  output$Frequency <- renderPlot({
    ggplot(zoom.data())+
      geom_point(aes(x = Time0, y = fs), color = "blue")+
      geom_point(data = regression.data(), aes(x = Time0, y = fs), color = "red")+
      labs(title = "Frequency of test",
           x = "Time",
           y = "Frequency")+
      # edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
    
  })
  # Graph Bandwidth data
  output$Bandwith <- renderPlot({
    ggplot(zoom.data())+
      geom_point(aes(x = Time0, y = Gamma), color = "blue")+
      geom_point(data = regression.data(), aes(x = Time0, y = Gamma), color = "red")+
      labs(title = "Bandwith of test",
           x = "Time",
           y = "Bandwith")+
      # edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  output$RegressionPlotFrequency <- renderPlot({
    ggplot(regression.data())+
      geom_point(aes(x = Time0, y = fs), color = "blue")+
      #geom_smooth(aes(x = Time0, y = fs), method = "lm", se=TRUE, color="black", formula =
      #y ~ log(x))+
      labs(title = "Frequency of test",
           x = "Time",
           y = "Frequency")+
      # edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  output$RegressionPlotBandwith <- renderPlot({
    ggplot(regression.data())+
      geom_point(aes(x = Time0, y = Gamma), color = "blue")+
      #geom_smooth(aes(x = Time0, y = Gamma), method = "lm", se=TRUE, color="black", formula =
      #y ~ log(x))+
      labs(title = "Bandwith of test",
           x = "Time",
           y = "Bandwith")+
      # edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  regression.data.filtered <- reactive({
    regression.data()%>%
      select(Time0,fs,Gamma)%>%
      mutate(dTime = 0.5)%>%
      mutate(dfs = 0.1)%>%
      mutate(dGamma = 0.1)
  })
  output$RegressionResultsFrequency <- renderTable(regression.data.filtered()%>%select(Time0,dTime,fs,dfs))
  output$RegressionResultsBandwith <- renderTable(regression.data.filtered()%>%select(Time0,dTime,Gamma,dGamma))
}


# Run the application 
shinyApp(ui = ui, server = server)
