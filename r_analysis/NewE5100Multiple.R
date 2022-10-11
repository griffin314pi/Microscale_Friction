#Include Libraries

library(shiny)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(readr)
library(gridExtra)
library(plotly)


#Define constants and parameters
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
    h1("E5100 Data Analysis (Multiple Datasets)", align = "left")
  ),
  fluidRow(column(width=4,offset=0,
                  div(style = "height:50px;", "by: Lucas Kramarczuk"))),
  
  fluidRow(column(width=4,
                  fileInput("files", NULL, buttonLabel = "Choose Files:", multiple = TRUE))),
  
  fluidRow(
    column(width = 12,
           plotOutput("stiffness", click = "plot_click2"))),
  fluidRow(
    column(width = 4, offset=4,
           verbatimTextOutput("info2"))),
  fluidRow(column(width = 12,
           plotOutput("deltaQinv"))),
  
  
  #Row 2
  fluidRow(
    column(width = 12,
           plotOutput("damping"))),
  fluidRow(column(width = 12,
           plotOutput("elastic", click = "plot_click"))),
  fluidRow(
    column(width = 4, offset = 4,
           verbatimTextOutput("info"))),
  
  fluidRow(
    column(width = 12,
           plotOutput("ELossGlobal"))),
  fluidRow(column(width = 12,
           plotOutput("ELossLinear"))),
  
  #Row 2 contains the COF and rsq of Energy loss linear regime
  fluidRow(
    column(width = 6, offset = 6,
           verbatimTextOutput("normal.load.vec")),
    column(width = 6, 
           verbatimTextOutput("KinForce"))),
  fluidRow(
    column(width = 6, offset = 6,
           verbatimTextOutput("COF")),
    column(width = 6,
           verbatimTextOutput("rsq"))),
  
  
  #Row 3 contains the slider for inputting the linear range of energy loss
  fluidRow(
    column(width = 12,
           sliderInput(width = 1200, inputId = "sliderRange", label = h3("Amplitude Range"), 
                       min = 0, max = 40, value = c(10, 20)))),
  
  #Row 4 shows the slider range values
  # fluidRow(
  #   column(width = 6, offset = 6, verbatimTextOutput("range"))),
  
  fluidRow(                
    column(width = 12,
           plotOutput("eStored"))),
  fluidRow(
    column(width = 12,
           plotOutput("lossTangent"))),
  
  # fluidRow(
  #   column(width=4,
  #          textInput("filename", "Results file name:"))),
  
  fluidRow(
    column(width=6,
           # downloadButton("downloadData", "Download results file")
           verbatimTextOutput("Results files coming soon...")),
#     column(width=6,
#            downloadButton("reorderData", "Download reordered data")
))
# )

server <- function(input, output) {
  
  #Make a list of all input datasets
  full.data.list <- reactive({
    req(input$files)
    upload <- list()
    
    for(i in 1:length(input$files[,1])){
      upload[[i]] <- read.delim(input$files[[i, 'datapath']], header=FALSE)
    }
    
    return(upload)
  })
  
  #Table of filenames and file IDs for baseline comparison between tests
  filename.tbl <- reactive({
    req(input$files)
    names <- c()

    for(i in 1:length(input$files[,1])){
      names[i] <- input$files[[i, 'name']]
    }
    df <- data.frame(ID=seq(1,length(names)),Test.ID=names)
    df
  })

  normal.load.vec <- reactive({
    vec <- c()
    for(i in 1:length(full.data.list())){
      testing.conditions <- full.data.list()[[i]][1,c(2:11)]
      vec[i] <- as.numeric(gsub(" uN", "", testing.conditions[7]))
    }
    vec
  })
  
  #Format the list of input datasets for plotting
  tidy.data.long <- reactive({
    #Initialize an empty list to store each dataset
    data.lst <- list()
    for(i in 1:length(full.data.list())){
      
      #Remove empty columns and remove testing condition rows
      tidy.tbl <- full.data.list()[[i]][-c(1,3),c(3:20)]
      #Set column names to first row
      colnames(tidy.tbl) <- tidy.tbl[1,]
      #Remove first row (which is now the column names)
      tidy.tbl <- tidy.tbl[-1,]
      
      #remove and parse time column for later
      segType <- tidy.tbl$Seg
      segID <- tidy.tbl$`ID Tag`
      
      tidy.tbl2 <- tidy.tbl[,c(3,5:18)]%>%
          #Change data type of all columns from character to numeric
          mutate_if(is.character,as.numeric)%>%
          #Compute gamma
          mutate(Gamma = R1/(4*pi*L1))%>%
          mutate(RMS.deltaV = `Ch1 (RMS)` - `Ch2 (RMS)`)%>%
          mutate_if(is.character,as.numeric)%>%
          mutate(segType=segType)%>%
          mutate(segID=segID)
      
      seg1 <- tidy.tbl2%>%
          filter(segID=="s1")%>%
          select(segType, Number, Power, C0, C1, L1, R1, fs, Gamma)%>%
          mutate(ID=i)%>%
          arrange(Power)
      
      seg2 <- tidy.tbl2%>%
          filter(segID=="s2 ")%>%
          select(segType, Number, Power, C0, C1, L1, R1, fs, Gamma, RMS.deltaV)%>%
          mutate(ID=i)%>%
          arrange(Power)
      
      seg3 <- tidy.tbl2%>%
          filter(segID=="s3")%>%
          select(segType, Number, Power, C0, C1, L1, R1, fs, Gamma)%>%
          mutate(ID=i)%>%
          arrange(Power)
      
      #Create final summary dataset with calculated metrics
      final.data <- data.frame(RMS.deltaV = seg2$RMS.deltaV,ID=seg2$ID)
      #delta fs and gamma
      #will use segType=B for caluclating these metrics
      final.data <-
        if(seg1$segType[1] == "B" & seg3$segType[1] != "B"){
          mutate(final.data,
                 deltaFs = seg2$fs - seg1$fs,
                 deltaGamma = seg2$Gamma - seg1$Gamma)
        } else if(seg3$segType[1] == "B" & seg1$segType[1] != "B"){
          mutate(final.data,
                 deltaFs = seg2$fs - seg3$fs,
                 deltaGamma = seg2$Gamma - seg3$Gamma)
        } else if(seg1$segType[1] == "B" & seg3$segType[1] == "B"){
          mutate(final.data,
                 deltaFs = seg2$fs - (seg1$fs+seg3$fs)/2,
                 deltaGamma = seg2$Gamma - (seg1$Gamma+seg3$Gamma)/2)
        }
      final.data <- final.data%>%
        #quality factor
        mutate(Q = 2*pi*seg2$fs*(seg2$L1/seg2$R1))%>%
        
        #absolutes
        
        #amplitude
        mutate(Amp = amplitude.factor*Q*sqrt(2)*RMS.deltaV/1000)%>%
        #Ki elastic
        mutate(Ki.elas = 2*(pi^2)*z.quartz*area.effective*deltaFs/1000)%>%
        #2pi fb
        mutate(TwoPi.fb = 2*(pi^2)*z.quartz*area.effective*deltaGamma/1000)%>%
        #Elastic Force
        mutate(F.elas = Ki.elas*Amp)%>%
        #Damping Force
        mutate(F.damp = TwoPi.fb*Amp)%>%
        #deltaE Elastic
        mutate(deltaE.elas = F.elas*Amp/2000)%>%
        #deltaE Damping
        mutate(deltaE.damp = 2*(pi^3)*z.quartz*area.effective*deltaGamma*(Amp*0.000000001)^2*1000000000000)%>%
        #deltaGamma/deltaFs
        mutate(delGamma.delFs = deltaGamma/deltaFs)%>%
        mutate(ID=i)
      
      data.lst[[i]] <- final.data
    }
    full.data <- bind_rows(data.lst)
    final.data <- full_join(filename.tbl(),full.data,by="ID")
  })
  
  #Clicked input for stiffness plot
  x2 <- reactiveVal(0)
  y2 <- reactiveVal(0)
  observeEvent(input$plot_click2,
               {
                 near <- nearPoints(tidy.data.long(), input$plot_click2)[1,]
                 newx2 <- near$Amp
                 x2(newx2)
               })
  observeEvent(input$plot_click2,
               {
                 near <- nearPoints(tidy.data.long(), input$plot_click2)[1,]
                 newy2 <- near$Ki.elas
                 y2(newy2)
               })
  observeEvent(input$plot_click2,
               {output$info2 <- renderText({
                 paste0("Amplitude: ", round(x2(),3), "\nStiffness: ", round(y2(),3))
               })})
  clicked2 <- reactive  ({ data.frame(x2=x2(),y2=y2()) })
  
  #Code for stiffness plot
  output$stiffness <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x = Amp, y = Ki.elas, color=Test.ID))+
      geom_point(data=clicked2(), aes(x2,y2),color="yellow")+
      labs(title = "Stiffness",
           x = "Amplitude (nm)",
           y = "K (μN/nm)")+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for delta Q inverse plot
  output$deltaQinv <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x = Amp, y = TwoPi.fb, color=Test.ID))+
      labs(title = expression(paste(Delta, Q^-1)),
           x = "Amplitude (nm)",
           y = expression(paste(Delta, Q^-1)))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for outputting elastic force maximum
  x <- reactiveVal(0)
  y <- reactiveVal(0)
  observeEvent(input$plot_click,
               {
                 near <- nearPoints(tidy.data.long(), input$plot_click)[1,]
                 newx <- near$Amp
                 x(newx)
               })
  observeEvent(input$plot_click,
               {
                 near <- nearPoints(tidy.data.long(), input$plot_click)[1,]
                 newy <- near$F.elas
                 y(newy)
               })
  observeEvent(input$plot_click,
               {output$info <- renderText({
                 paste0("Amplitude: ", round(x(),3), "\nForce:     ", round(y(),3))
               })})
  clicked <- reactive  ({ data.frame(x=x(),y=y()) })
  
  #Code for average elastic force plot
  output$elastic <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x = Amp,y = F.elas,color=Test.ID))+
      geom_point(data=clicked(), aes(x,y),color="Yellow")+
      labs(title = "Average Elastic Force",
           x = "Amplitude (nm)",
           y = "F (μN)")+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for average damping force plot
  output$damping <- renderPlot({
    ggplot()+
      geom_point(data=tidy.data.long(),aes(x = Amp, y = F.damp,color=Test.ID))+
      labs(title = "Average Damping Force",
           x = "Amplitude (nm)",
           y = "F (μN)")+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for energy stored vs. amplitude plot
  output$eStored <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x = Amp, y = deltaE.elas,color=Test.ID))+
      labs(title = "Energy Stored",
           x = "Amplitude (nm)",
           y = expression(paste(Delta,E," (pJ)")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for loss tangent plot
  output$lossTangent <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x = Amp, y = delGamma.delFs,color=Test.ID))+
      labs(title = "Loss Tangent",
           x = "Amplitude (nm)",
           y = "")+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for global energy loss vs. amplitude plot
  output$ELossGlobal <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x = Amp, y = deltaE.damp,color=Test.ID))+
      labs(title = "Energy Lost per Cycle (Global)",
           x = "Amplitude (nm)",
           y = expression(paste(Delta,E," (pJ)")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=21,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #extract slider min and max for output file
  sliderMin <- reactive ({ input$sliderRange[1] })
  sliderMax <- reactive ({ input$sliderRange[2] })
  
  #Subsetting data for energy loss linear regime
  filtered.data <- reactive({ #reactive object updates every time the user changes the slider
    req(input$sliderRange)
    filter(tidy.data.long(), Amp>input$sliderRange[1] & Amp<input$sliderRange[2])
  })

maxID <- reactive({ 
  obj <- tidy.data.long()$ID
  max(obj)
  })  
  
#vectors of lm coefficients
b.vec <- reactive({ 
  lm.list <- list()
  b.vec <- c()
  for (i in 1:maxID()){
  #Make linear model for energy lost vs. amplitude based on selected amplitude range
  filter.long <- filter(tidy.data.long(),ID==i)
  lm.list[[i]] <- lm(deltaE.damp ~ Amp, filter.long)
  #Extract intercept
  b.vec[i] <- round(unname(lm.list[[i]]$coefficients[1]),3) 
  }
  b.vec
    })
  
m.vec <- reactive({ 
  lm.list <- list()
  m.vec <- c()
  for (i in 1:maxID()){
    #Make linear model for energy lost vs. amplitude based on selected amplitude range
    filter.long <- filter(tidy.data.long(),ID==i)
    lm.list[[i]] <- lm(deltaE.damp ~ Amp, filter.long)
    #Extract intercept
    m.vec[i] <- round(unname(lm.list[[i]]$coefficients[2]),3) 
  }
  m.vec
  })

kinetic.friction.force <- reactive ({ 
  lm.list <- list()
  force.vec <- c()
  for (i in 1:maxID()){
    #Make linear model for energy lost vs. amplitude based on selected amplitude range
    filter.long <- filter(tidy.data.long(),ID==i)
    lm.list[[i]] <- lm(deltaE.damp ~ Amp, filter.long)
    #Extract intercept
    force.vec[i] <- round(1000*unname(lm.list[[i]]$coefficients[2])/4,3) 
  }
  force.vec
  })

cof.vec <- reactive({ 
  lm.list <- list()
  cof.vec <- c()
  for (i in 1:maxID()){
    #Make linear model for energy lost vs. amplitude based on selected amplitude range
    filter.long <- filter(tidy.data.long(),ID==i)
    lm.list[[i]] <- lm(deltaE.damp ~ Amp, filter.long)
    #Extract intercept
    cof.vec[i] <- round(1000*unname(lm.list[[i]]$coefficients[2])/4/normal.load.vec()[i],3) 
  }
  cof.vec
  })

rsq.vec <- reactive ({ 
  lm.list <- list()
  rsq.vec <- c()
  for (i in 1:maxID()){
    #Make linear model for energy lost vs. amplitude based on selected amplitude range
    filter.long <- filter(tidy.data.long(),ID==i)
    lm.list[[i]] <- lm(deltaE.damp ~ Amp, filter.long)
    #Extract intercept
    rsq.vec[i] <- round(summary(lm.list[[i]])$r.squared,3)
  }
  rsq.vec
  })

#Test the output
output$rsq <- reactive ({ rsq.vec() })
output$COF <- reactive ({ cof.vec() })

  #Code for energy loss linear regime plot
  output$ELossLinear <- renderPlot ({
    ggplot(filtered.data())+ #use the reactive filtered data for plotting
      geom_point(aes(x = Amp, y = deltaE.damp,color=Test.ID))+
      geom_smooth(aes(x = Amp, y = deltaE.damp,color=Test.ID),method = 'lm',
                  se=FALSE, formula = y ~ x)+
      labs(title = "Energy Lost per Cycle Linear Regime",
           x = "Amplitude",
           y = expression(paste(Delta,E," (pJ)")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=21,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  
}

shinyApp(ui, server)

