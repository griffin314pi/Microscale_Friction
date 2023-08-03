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
amplitude.factornoamplitude <- 0.012247

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
  # Helper text
  fluidRow(column(width=6,
                  h6("Use the slider below to select the baseline range")),
           column(width = 6,
                  h6("Use the slider below to select the test range"))
  ),
  # Row 2 - Slider bars for selecting baseline and test range
  fluidRow(
    column(width = 6,
           sliderInput(width = 600, inputId = "BaselineSlider", label = h3("Index Range Baseline (orange)"), 
                       min = 0, max = 100, value = c(40, 50))),
    column(width = 6,
           sliderInput(width = 600, inputId = "TestSlider", label = h3("Index Range Test (green)"),
                       min = 0, max = 100, value = c(50, 60)))
  ),
  fluidRow(
    column(width = 6, offset = 0,
           tableOutput(outputId = "dataTable")),
    column(width = 6,
           plotOutput("Power"))
  ),
  #File name input for results file
  fluidRow(
    column(width=4,
           textInput("filename", "Results file name:"))),
  fluidRow(
    column(width=6,
           downloadButton("downloadData", "Download results file")))
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
  maxIndex <- reactive ({ as.numeric(gsub("\\D", "", testing.conditions()[10]))})
  
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
      mutate(PoweruW=Power*10^6)
  })
  observeEvent(!is.null(tidy.tbl2()),{
    maxpoint <- rev(tidy.tbl2()$Number)[1]%>%as.numeric
    updateSliderInput(inputId = "ZoomSlider", max = maxpoint,value = c(0,maxpoint))
  })
  observeEvent(input$ZoomSlider,{
    zoomMin <- input$ZoomSlider[1]
    zoomMax <- input$ZoomSlider[2]
    updateSliderInput(inputId = "BaselineSlider", min = zoomMin, max = zoomMax)
    updateSliderInput(inputId = "TestSlider", min = zoomMin, max = zoomMax)
  })
  # extract slider min and max for baseline
  BaselineMin <- reactive ({ input$BaselineSlider[1] })
  BaselineMax <- reactive ({ input$BaselineSlider[2] })
  #extract slider min and max for test
  TestMin <- reactive ({ input$TestSlider[1] })
  TestMax <- reactive ({ input$TestSlider[2] })
  
  # Extracting zoomed in data
  zoom.data <- reactive({ #reactive object updates every time the user changes the slider
    req(input$BaselineSlider)
    filter(tidy.tbl2(), Number>input$ZoomSlider[1] & Number<input$ZoomSlider[2])
  })
  # Subsetting baseline data
  filteredBaseline.data <- reactive({ #reactive object updates every time the user changes the slider
    req(input$BaselineSlider)
    filter(tidy.tbl2(), Number>input$BaselineSlider[1] & Number<input$BaselineSlider[2])
  })
  # Subsetting test data
  filteredTest.data <- reactive({ #reactive object updates every time the user changes the slider
    req(input$TestSlider)
    filter(tidy.tbl2(), Number>input$TestSlider[1] & Number<input$TestSlider[2])
  })
  # Filtering first datapoints in each segment
  firstPoint.data <- reactive({
    filter(zoom.data(), Number %% (maxIndex()+1) == 0)
  })
  
  # Graph Frequency data
  output$Frequency <- renderPlot({
    ggplot(zoom.data())+
      geom_point(aes(x = Number, y = fs), color = "blue")+
      geom_point(data = filteredBaseline.data(), aes(x = Number, y = fs), color = "orange")+
      geom_point(data = filteredTest.data(), aes(x = Number, y = fs), color = "green")+
      geom_point(data = firstPoint.data(), aes(x = Number, y = fs), color =  "red")+
      labs(title = "Frequency of test",
           x = "Index",
           y = "Frequency")+
      # edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
    
  })
  output$Bandwith <- renderPlot({
    ggplot(zoom.data())+
      geom_point(aes(x = Number, y = Gamma), color = "blue")+
      geom_point(data = filteredBaseline.data(), aes(x = Number, y = Gamma), color = "orange")+
      geom_point(data = filteredTest.data(), aes(x = Number, y = Gamma), color = "green")+
      geom_point(data = firstPoint.data(), aes(x = Number, y = Gamma), color =  "red")+
      labs(title = "Bandwith of test",
           x = "Index",
           y = "Bandwith")+
      # edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  output$Power <- renderPlot({
    ggplot(zoom.data())+
      geom_point(aes(x = Number, y = PoweruW), color = "blue")+
      geom_point(data = filteredBaseline.data(), aes(x = Number, y = PoweruW), color = "orange")+
      geom_point(data = filteredTest.data(), aes(x = Number, y = PoweruW), color = "green")+
      geom_point(data = firstPoint.data(), aes(x = Number, y = PoweruW), color =  "red")+
      labs(title = "Power of test",
           x = "Index",
           y = "Power uW")+
      # edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  # Calculating friction forces from test points with st. dev.
  results <- reactive({
    # Finding averages for Baseline and Test
    meanBaselineFreq <- mean(filteredBaseline.data()$fs)
    meanBaselineGamma <- mean(filteredBaseline.data()$Gamma)
    meanBaselinePower <- mean(filteredBaseline.data()$Power)
    meanBaselinePoweruW <- mean(filteredBaseline.data()$PoweruW)
    meanTestFreq <- mean(filteredTest.data()$fs)
    meanTestGamma <- mean(filteredTest.data()$Gamma)
    meanTestL1 <- mean(filteredTest.data()$L1)
    meanTestR1 <- mean(filteredTest.data()$R1)
    meanTestRMS.deltaV <- mean(filteredTest.data()$RMS.deltaV)
    meanTestPower <- mean(filteredTest.data()$Power)
    meanTestPoweruW <- mean(filteredTest.data()$PoweruW)
    # Finding standard deviation for Frequency and Gamma
    sdBaselineFreq <- sd(filteredBaseline.data()$fs)
    sdBaselineGamma <- sd(filteredBaseline.data()$Gamma)
    sdBaselinePower <- sd(filteredBaseline.data()$Power)
    sdBaselinePoweruW <- sd(filteredBaseline.data()$PoweruW)
    sdTestFreq <- sd(filteredTest.data()$fs)
    sdTestGamma <- sd(filteredTest.data()$Gamma)
    sdTestL1 <- sd(filteredTest.data()$L1)
    sdTestR1 <- sd(filteredTest.data()$R1)
    sdTestRMS.deltaV <- sd(filteredTest.data()$RMS.deltaV)
    sdTestPower <- sd(filteredTest.data()$Power)
    sdTestPoweruW <- sd(filteredTest.data()$PoweruW)
    
    
    #Finding frequency shifts
    deltaFs <- meanTestFreq - meanBaselineFreq
    sddeltaFs <- sdTestFreq + sdBaselineFreq
    deltaGamma <- meanTestGamma - meanBaselineGamma
    sddeltaGamma <- sdTestGamma + sdBaselineGamma
    
    # Quality factor
    Q <- 2*pi*meanTestFreq*(meanTestL1/meanTestR1)
    sdQ <- Q*sqrt((sdTestFreq/meanTestFreq)^2+(sdTestL1/meanTestL1)^2+(sdTestR1/meanTestR1)^2)
    
    # Amplitude
    if (meanTestRMS.deltaV == 0){
      Amp <- amplitude.factornoamplitude*Q*sqrt(meanTestPower)
      sdAmp <- Amp*sqrt((sdQ/Q^2)^2+(sdTestPower/(2*sqrt(meanTestPower)))^2)
    }
    else {
      Amp <- amplitude.factor*Q*sqrt(2)*meanTestRMS.deltaV/1000
      sdAmp <- Amp*sqrt((sdQ/Q)^2+(sdTestRMS.deltaV/meanTestRMS.deltaV)^2)
    }
    # Ki elastic
    Ki.elas <- 2*(pi^2)*z.quartz*area.effective*deltaFs/1000
    sdKi.elas <- 2*(pi^2)*z.quartz*area.effective*sddeltaFs/1000
    # 2pi fb
    TwoPi.fb <- 2*(pi^2)*z.quartz*area.effective*deltaGamma/1000
    sdTwoPi.fb <- 2*(pi^2)*z.quartz*area.effective*sddeltaGamma/1000
    # Elastic Force
    F.elas <- Ki.elas*Amp
    sdF.elas <- F.elas*sqrt((sdKi.elas/Ki.elas)^2+(sdAmp/Amp)^2)
    # Damping Force
    F.damp <- TwoPi.fb*Amp
    sdF.damp <- F.damp*sqrt((sdTwoPi.fb/TwoPi.fb)^2+(sdAmp/Amp)^2)
    # deltaE Elastic
    deltaE.elas <- F.elas*Amp/2000
    sddeltaE.elas <- deltaE.elas*sqrt((sdF.elas/F.elas)^2+(sdAmp/Amp)^2)
    # deltaE Damping
    deltaE.damp <- 2*(pi^3)*z.quartz*area.effective*deltaGamma*(Amp*0.000000001)^2*1000000000000
    sddeltaE.damp <- deltaE.damp*sqrt((sddeltaGamma/deltaGamma)^2+(2*sdAmp/Amp)^2)
    # deltaGamma/deltaFs
    delGamma.delFs <- deltaGamma/deltaFs
    sddelGamma.delFs <- delGamma.delFs*sqrt((sddeltaGamma/deltaGamma)^2+(sddeltaFs/deltaFs)^2)
    # damping force/load (ratio)
    Fdamp.load <- F.damp/normal.load()
    sdFdamp.load <- Fdamp.load*sdF.damp/F.damp
    resultstbl <- data_frame(Data = c("Baseline Frequency", "Baseline Gamma", "Baseline Power (uW)", "Test Frequency", "Test Gamma", "Test Power (uW)", "Test L1", "Test R1", "Test RMS.deltaV", "Delta Frequency", "Delta Gamma", "Quality Factor", "Amplitude", "Ki Elastic", "2 pi Fb", "Elastic Force", "Damping Force", "Delta E Elastic", "Delta E Damping", "Delta Gamma/Delta Frequency", "Ratio damping and load"))
    resultstbl <- resultstbl%>%
      mutate(Mean = c(meanBaselineFreq, meanBaselineGamma, meanBaselinePoweruW, meanTestFreq, meanTestGamma, meanTestPoweruW, meanTestL1, meanTestR1, meanTestRMS.deltaV, deltaFs, deltaGamma, Q, Amp, Ki.elas, TwoPi.fb, F.elas, F.damp, deltaE.elas, deltaE.damp, delGamma.delFs, Fdamp.load))%>%
      mutate(StDev = c(sdBaselineFreq, sdBaselineGamma, sdBaselinePoweruW, sdTestFreq, sdTestGamma, sdTestPoweruW, sdTestL1, sdTestR1, sdTestRMS.deltaV, sddeltaFs, sddeltaGamma, sdQ, sdAmp, sdKi.elas, sdTwoPi.fb, sdF.elas, sdF.damp, sddeltaE.elas, sddeltaE.damp, sddelGamma.delFs, sdFdamp.load))
    
  })
  output$dataTable <- renderTable(results())
  
  #Code for downloading results table
  output$downloadData <- downloadHandler(
    filename = function() {paste(input$filename, ".csv", sep = "")
    },
    content = function(file){
      write.csv(results(), file, row.names = FALSE)
      #capture.output(results.data(), file = "my_list.txt")
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
