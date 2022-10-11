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
    h1("E5100 Data Analysis", align = "left")
  ),
  #Page will be in row format
  fluidRow(column(width=4,offset=0,
                  div(style = "height:50px;", "by: Lucas Kramarczuk"))),
  #File input
  fluidRow(column(width=4,
                  fileInput("upload", NULL, buttonLabel = "Choose File:", multiple = FALSE))),
  #Helper text
  fluidRow(column(width=6,offset=0,
                  h6("Click on the y intercept of stiffness vs. amplitude")
  )),
  #Row 1 - Stiffness and deltaQ inverse
  fluidRow(
    column(width = 6,
           plotOutput("stiffness", click = "plot_click2")),
    column(width = 6,
           plotOutput("deltaQinv"))),
  #Text under plots displaying clicked data and stiffness calculations
  fluidRow(
    column(width = 3,
           verbatimTextOutput("stiffVec")),
    column(width = 3,
           verbatimTextOutput("info2"))),
  fluidRow(
    column(width = 3,
           verbatimTextOutput("stiffVecSeg1"))),
  
  
  #Helper text
  fluidRow(column(width=6,offset=6,
                  h6("Click and drag your cursor to estimate the maximum elastic force ")
  )),
  #Row 2 - Damping force and Elastic force
  fluidRow(
    column(width = 6,
           plotOutput("damping")),
    column(width = 6,
           plotOutput("elastic", brush = "plot_brush"))),
  #Info from brushed points
  fluidRow(
    column(width = 4, offset = 7,
           verbatimTextOutput("info"))),
  
  #Helper text
  fluidRow(column(width=6,offset=6,
                  h6("Use the slider below to select the linear range of Energy Lost per Cycle")
  )),
  
  #Row 3 - energy loss global and local
  fluidRow(
    column(width = 6,
           plotOutput("ELossGlobal")),
    column(width = 6,
           plotOutput("ELossLinear"))),
  
  #Output normal load and COF linear model details
  fluidRow(
    column(width = 3, offset = 6,
           verbatimTextOutput("Load")),
    column(width = 3, 
           verbatimTextOutput("KinForce"))),
  fluidRow(
    column(width = 3, offset = 6,
           verbatimTextOutput("COF")),
    column(width = 3,
           verbatimTextOutput("rsq"))),
  
  #Slider range UI
  fluidRow(
    column(width = 6, offset = 6,
           sliderInput(width = 600, inputId = "sliderRange", label = h3("Amplitude Range"), 
                       min = 0, max = 50, value = c(10, 20)))),
  #Row 4 - energy stored and loss tangent
  fluidRow(                
    column(width = 6,
           plotOutput("eStored")),
    column(width = 6,
           plotOutput("lossTangent"))),
  #File name input for results file
  fluidRow(
    column(width=4,
           textInput("filename", "Results file name:"))),
  #Helper text
  fluidRow(actionButton("show", "Ready to Download? \nClick This Button")),
  #Download buttons at bottom of page
  fluidRow(
    column(width=6,
           downloadButton("downloadData", "Download results file")),
    column(width=6,
           downloadButton("reorderData", "Download reordered data")))
)

# Define server logic required to make the plots
server <- function(input, output) {
  
  #Warning message to make sure users did everything correctly
  observeEvent(input$show, {
    showModal(modalDialog(
      title = "Important message",
      HTML("Before downloading the results, ensure that:<br> 
      1) The correct stiffness intercept has been selected<br>
      2) The correct points have been highlighted for the maximum elastic force<br>
      3) The linear range of energy lost per cycle is accurate")
    ))
  })
  
  #Code for file input
  test.tbl <- reactive({
    req(input$upload)
    inFile <- input$upload
    tbl <- read.delim(inFile$datapath, header=FALSE)
    return(tbl)
  })
  
  #Extract the testing conditions and units of measurement from input data
  testing.conditions <- reactive ({ test.tbl()[1,c(2:11)] })
  units.of.measurement <- reactive ({ test.tbl()[c(2,3),c(3:20)] })
  #print(testing.conditions)
  #print(units.of.measurement)
  normal.load <- reactive ({ as.numeric(gsub(" uN", "", testing.conditions()[7])) })
  output$Load <- reactive ({ paste("Normal Load:",normal.load(),"uN") })
  
  #Code for tidying the input dataset
  
  #Remove empty columns and remove testing condition rows
  tidy.tbl <- reactive({ new.tbl <- test.tbl()[-c(1,3),c(3:20)]
  #Set column names to first row
  colnames(new.tbl) <- new.tbl[1,] 
  #Remove first row (which is now the column names)
  new.tbl[-1,]
  })
  
  #Add gamma and RMS delta V columns to full dataset
  
  #remove and parse time column for later
  Time.vec <- reactive({ parse_time(tidy.tbl()$Time) })
  segType <- reactive ({ tidy.tbl()$Seg })
  segID <- reactive ({ tidy.tbl()$`ID Tag` })
  
  tidy.tbl2 <- reactive ({ tidy.tbl()[,c(3,5:18)]%>%
      #Change data type of all columns from character to numeric
      mutate_if(is.character,as.numeric)%>%
      #Compute gamma
      mutate(Gamma = R1/(4*pi*L1))%>%
      mutate(RMS.deltaV = `Ch1 (RMS)` - `Ch2 (RMS)`)%>%
      mutate(Time = Time.vec())%>%
      mutate(segType=segType())%>%
      mutate(segID=segID())
  })
  
  #Create separate datasets for segments 1,2,3
  seg1 <- reactive({ tidy.tbl2()%>%
      filter(segID=="s1")%>%
      select(segType, Number, Power, C0, C1, L1, R1, fs, Gamma)%>%
      arrange(Power)
  })
  
  seg2 <- reactive({ tidy.tbl2()%>%
      filter(segID=="s2 ")%>%
      select(segType, Number, Power, C0, C1, L1, R1, fs, Gamma, RMS.deltaV)%>%
      arrange(Power)
  })
  
  seg3 <- reactive({ tidy.tbl2()%>%
      filter(segID=="s3")%>%
      select(segType, Number, Power, C0, C1, L1, R1, fs, Gamma)%>%
      arrange(Power)
  })
  
  #Create final summary dataset with calculated metrics
  final.data <- reactive({ data.tbl <- data.frame(RMS.deltaV = seg2()$RMS.deltaV)
  #delta fs and gamma
  #this part of the code detects which segment is the baseline 
  #and calculates gamma and frequency shifts accordingly
  data.tbl <-
    if(seg1()$segType[1] == "B" & seg3()$segType[1] != "B"){
      mutate(data.tbl,
             deltaFs = seg2()$fs - seg1()$fs,
             deltaGamma = seg2()$Gamma - seg1()$Gamma)
    } else if(seg3()$segType[1] == "B" & seg1()$segType[1] != "B"){
      mutate(data.tbl,
             deltaFs = seg2()$fs - seg3()$fs,
             deltaGamma = seg2()$Gamma - seg3()$Gamma)
    } else if(seg1()$segType[1] == "B" & seg3()$segType[1] == "B"){
      mutate(data.tbl,
             deltaFs = seg2()$fs - (seg1()$fs+seg3()$fs)/2,
             deltaGamma = seg2()$Gamma - (seg1()$Gamma+seg3()$Gamma)/2)
    }
  data.tbl <- data.tbl%>%
    #quality factor
    mutate(Q = 2*pi*seg2()$fs*(seg2()$L1/seg2()$R1))%>%
    
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
    mutate(deltaE.damp = 
             2*(pi^3)*z.quartz*area.effective*deltaGamma*(Amp*0.000000001)^2*1000000000000)%>%
    #deltaGamma/deltaFs
    mutate(delGamma.delFs = deltaGamma/deltaFs)%>%
    arrange(Amp)
  data.tbl
  })
  
  #Stiffness plot click
  observeEvent(input$plot_click2,
               {output$info2 <- renderText({
                 paste0("Clicked Point:","\nAmplitude: ", round(x2(),3), "\nStiffness: ", round(y2(),3))
               })})
  
  #Clicked input for stiffness plot
  x2 <- reactiveVal(0)
  y2 <- reactiveVal(0)
  
  #Set the x and y values closest to the click
  observeEvent(input$plot_click2,
               {
                 near <- nearPoints(final.data(), input$plot_click2)[1,]
                 newx2 <- near$Amp
                 x2(newx2)
               })
  observeEvent(input$plot_click2,
               {
                 near <- nearPoints(final.data(), input$plot_click2)[1,]
                 newy2 <- near$Ki.elas
                 y2(newy2)
               })
  #Store x and y values as a 1x2 data frame
  clicked2 <- reactive  ({ data.frame(x2=x2(),y2=y2()) })
  
  #First 5 stiffness data points (test segment)
  first5 <- reactive ({ 
    data.frame(Amp=final.data()$Amp[c(1:5)],Ki.elas=final.data()$Ki.elas[c(1:5)])
  })
  
  #First 5 stiffness data points of seg1
  first5seg1 <- reactive ({ 
    deltaFs <- seg1()$fs - seg3()$fs
    data.frame(Amp=final.data()$Amp[c(1:5)],
               Ki.elas=2*(pi^2)*z.quartz*area.effective*deltaFs[c(1:5)]/1000)
  })
  
  #Code for stiffness plot
  output$stiffness <- renderPlot({
    ggplot(final.data())+
      geom_point(aes(x = Amp, y = Ki.elas), color = "blue")+
      geom_point(data=first5(), aes(Amp,Ki.elas),color="green")+
      geom_point(data=first5seg1(), aes(Amp,Ki.elas),color="orange")+
      geom_point(data=clicked2(), aes(x2,y2),color="red")+
      labs(title = "Stiffness",
           x = "Amplitude (nm)",
           y = "K (μN/nm)")+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Average and standard deviation of first 5 stiffness points
  stiffness.vec <- reactive({
    data <- final.data()$Ki.elas[c(1:5)]
    mean <- round(mean(data),3)
    sd <- round(sd(data),3)
    paste("Green Data (seg2):","\nMean: ",mean,"\nStDev:",sd)
  })
  
  stiff.mean <- reactive({
    data <- final.data()$Ki.elas[c(1:5)]
    round(mean(data),3)
  })
  
  stiff.sd <- reactive({
    data <- final.data()$Ki.elas[c(1:5)]
    round(sd(data),3)
  })
  
  output$stiffVec <- reactive ({ 
    stiffness.vec()
  })
  
  #Segment 1 stiffness  
  stiffness.seg1.vec <- reactive({
    data <- first5seg1()$Ki.elas
    mean <- round(mean(data),3)
    sd <- round(sd(data),3)
    paste("Orange Data (seg1):","\nMean: ",mean,"\nStDev:",sd)
  })
  
  output$stiffVecSeg1 <- reactive ({ 
    stiffness.seg1.vec()
  })
  
  #Code for delta Q inverse plot
  output$deltaQinv <- renderPlot({
    ggplot(final.data())+
      geom_point(aes(x = Amp, y = TwoPi.fb), color = "blue")+
      labs(title = expression(paste(Delta, Q^-1)),
           x = "Amplitude (nm)",
           y = expression(paste(Delta, Q^-1)))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for outputting elastic force maximum
  observeEvent(input$plot_brush,
               {output$info <- renderText({
                 paste0("Mean Force:     ", round(static.friction.force(),3),"\nMean Amplitude: ", 
                        round(x(),3),
                        "\nStatic COF:     ",static.COF())
               })})
  
  #Create data frame of highlighted points for elastic force
  brushed.tbl <- reactive ({ 
    brushed <- brushedPoints(final.data(), input$plot_brush)
    x <- brushed$Amp
    y <- brushed$F.elas
    data_frame(Amp=x,F.elas=y)
  })
  
  #Code for highlighting the brushed points on elastic force plot
  x <- reactive({ mean(brushed.tbl()$Amp )})
  static.friction.force <- reactive({ mean(brushed.tbl()$F.elas )})
  
  #StdDev of highlighted points for elastic force
  xsd <- reactive({ sd(brushed.tbl()$Amp )})
  ysd <- reactive({ sd(brushed.tbl()$F.elas )})
  
  #Code for average elastic force plot
  output$elastic <- renderPlot({
    ggplot()+
      geom_point(data=final.data(),aes(x = Amp, y = F.elas),color="Blue")+
      geom_point(data=brushed.tbl(), aes(Amp,F.elas),color="Red")+
      labs(title = "Average Elastic Force",
           x = "Amplitude (nm)",
           y = "F (μN)")+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Code for extracting static COF from elastic force brushed points
  static.COF <- reactive({
    val<-round(static.friction.force()/normal.load(),3)
  })
  
  #Code for average damping force plot
  output$damping <- renderPlot({
    ggplot()+
      geom_point(data=final.data(),aes(x = Amp, y = F.damp), color = "blue")+
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
    ggplot(final.data())+
      geom_point(aes(x = Amp, y = deltaE.elas), color = "blue")+
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
    ggplot(final.data())+
      geom_point(aes(x = Amp, y = delGamma.delFs), color = "blue")+
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
    ggplot(final.data())+
      geom_point(aes(x = Amp, y = deltaE.damp), color = "blue")+
      labs(title = "Energy Lost per Cycle (Global)",
           x = "Amplitude (nm)",
           y = expression(paste(Delta,E," (pJ)")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=21,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #extract slider min and max for output file linear range
  sliderMin <- reactive ({ input$sliderRange[1] })
  sliderMax <- reactive ({ input$sliderRange[2] })
  
  #Subsetting data for energy loss linear regime
  filtered.data <- reactive({ #reactive object updates every time the user changes the slider
    req(input$sliderRange)
    filter(final.data(), Amp>input$sliderRange[1] & Amp<input$sliderRange[2])
  })
  
  #Make linear model for energy lost vs. amplitude based on selected amplitude range
  lin.model <- reactive({ lm(deltaE.damp ~ Amp, filtered.data()) })
  
  #Extract slope and intercept
  b <- reactive({ unname(lin.model()$coefficients[1]) })
  m <- reactive({ unname(lin.model()$coefficients[2]) })
  
  #Calculate COF from linear model coefficient
  kinetic.friction.force<-reactive({round(1000*m()/4,3)})
  output$KinForce <- reactive ({ paste("Fric. Force:",round(kinetic.friction.force(),3),"uN") })
  cof <- reactive({ round(kinetic.friction.force()/normal.load(),3) })
  output$COF <- reactive({ 
    as.character(sprintf("Kinetic COF = %s", 
                         cof()))
  })
  
  #Extract r squared from the linear model
  rsquared <- reactive ({ round(summary(lin.model())$r.squared,3) })
  output$rsq <- reactive({ as.character(sprintf("Rsq = %s",
                                                round(summary(lin.model())$r.squared,3))) })
  
  #Code for energy loss linear regime plot
  output$ELossLinear <- renderPlot ({
    ggplot(filtered.data())+ #use the reactive filtered data for plotting
      geom_point(aes(x = Amp, y = deltaE.damp), color = "blue")+
      geom_smooth(aes(x = Amp, y = deltaE.damp), method = "lm", se=FALSE, color="black", formula =
                    y ~ x)+
      labs(title = "Energy Lost per Cycle Linear Regime",
           x = "Amplitude",
           y = expression(paste(Delta,E," (pJ)")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=21,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Show slider range
  output$range <- renderPrint({input$sliderRange})
  
  #Code for making reordered dataset
  reordered.data <- reactive({
    cbind(seg1(),seg2(),seg3())
  })
  
  #Code for dowmloading reordered data
  output$reorderData <- downloadHandler(
    filename = function() {paste("Reordered ", input$upload$name, ".csv", sep = "")
    },
    content = function(file){
      write.csv(reordered.data(), file, row.names = FALSE)
    }
  )
  
  #Code for dowmloading results table
  output$downloadData <- downloadHandler(
    filename = function() {paste(input$filename, ".csv", sep = "")
    },
    content = function(file){
      write.csv(results.data(), file, row.names = FALSE)
      #capture.output(results.data(), file = "my_list.txt")
    }
  )
  
  #Code for creating the results table
  #If you are reading this good luck
  results.data <- reactive ({
    names <- c("Filename","Normal.Load (uN)","Kinetic.Friction.Force (uN)",
               "COF.Kinetic","Static.Friction.Force (uN)","COF.Static",
               "max.elas.force.sd","m","b","rsq","lin.amp.min",
               "lin.amp.max","Stiffness.First5.Avg","Stiffness.First5.SD","Stiffness.k0")
    vals <- reactive ({
      c(input$upload$name,normal.load(),kinetic.friction.force(),cof(),
        x(),static.COF(),xsd(),m(),b(),rsquared(),sliderMin(),
        sliderMax(),stiff.mean(),stiff.sd(),y2())
    })
    #Assign static objects to reactives for changing lengths
    vals2 <- reactive({
      filtered.data()$Amp
    })
    vals3 <- reactive({
      filtered.data()$deltaE.damp
    })
    vals4 <-   
      reactive({
        final.data()$Amp
      })
    vals5 <- reactive({
      final.data()$F.elas
    })
    
    #Make vectors same length by adding NAs
    length(names) <- max(length(names),length(vals()),length(vals2()),length(vals3())
                         ,length(vals4()),length(vals5()))
    valu <- vals()
    length(valu) <- max(length(names),length(vals()),length(vals2()),length(vals3())
                        ,length(vals4()),length(vals5()))
    valu2 <- vals2()
    length(valu2) <- max(length(names),length(vals()),length(vals2()),length(vals3())
                         ,length(vals4()),length(vals5()))
    valu3 <- vals3()
    length(valu3) <- max(length(names),length(vals()),length(vals2()),length(vals3())
                         ,length(vals4()),length(vals5()))
    valu4 <- vals4()
    length(valu4) <- max(length(names),length(vals()),length(vals2()),length(vals3())
                         ,length(vals4()),length(vals5()))
    valu5 <- vals5()
    length(valu5) <- max(length(names),length(vals()),length(vals2()),length(vals3())
                         ,length(vals4()),length(vals5()))
    #Combine vectors into data frame
    df <- data.frame(Variables=names,Results=valu,Amp.Linear=valu2,E.Loss.Linear=valu3,
                     Amp.Global=valu4,Elas.Force=valu5)
    #Set NAs to blanks
    df[is.na(df)] <- ""
    #Return the data frame
    df
  })
  
  #Phew that's it
  
}

# Run the application 
shinyApp(ui = ui, server = server)


