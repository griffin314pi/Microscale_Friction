#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

#Define constants and parameters
z.quartz <- 8800000 
area.effective <- 0.0025^2*pi
amplitude.factor <- 1.4 #pm/V
d26 <- 3.1*0.000000000001

U.star <- 1
F.elas.star = 1
F.k = 1
K = 1
G.star = 14
a = K/(8*G.star)

ui <- fluidPage(
  fluidRow(column(width=4,
                  fileInput("files", NULL, buttonLabel = "Choose Files:", multiple = TRUE)
  )),
  fluidRow(
    tableOutput("contents")
  ),
  
  fluidRow(
    plotOutput("freq")
  ),
  fluidRow(
    plotOutput("dFs")
  ),
  fluidRow(
    plotOutput("resistance")
  ),
  fluidRow(
    plotOutput("dR")
  ),
  fluidRow(
    tableOutput("table")
  ),
  
)

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
      names[i] <- substr(names[i],1,25)
    }
    df <- data.frame(ID=seq(1,length(names)),Test.ID=names)
    df
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
        mutate_if(is.character,as.numeric)%>%
        #Compute gamma
        mutate(Gamma = R1/(4*pi*L1))%>%
        mutate(segType=segType)%>%
        mutate(segID=segID)
      
      baseline.tbl <- tidy.tbl2%>%
        select(segType, segID, Number, Power, R1, L1, fs,`Ch1 (RMS)`,`Ch2 (RMS)`,Gamma)%>%
        mutate(ID=i)
        # mutate(L1 = as.numeric(L1))%>%
        # mutate(R1 = as.numeric(R1))%>%
        # mutate(`Ch1 (RMS)` = as.numeric(`Ch1 (RMS)`))%>%
        # mutate(`Ch2 (RMS)` = as.numeric(`Ch2 (RMS)`))%>%
        # mutate(Gamma = R1/(4*pi*L1))%>%
        # mutate(RMS.deltaV = `Ch1 (RMS)` - `Ch2 (RMS)`)%>%
        # mutate(Power = as.numeric(Power))%>%
        # mutate(Number = as.numeric(Number))%>%
        # mutate(R1 = as.numeric(R1))%>%
        # mutate(fs = as.numeric(fs))
      
      s1 <- dplyr::filter(baseline.tbl,segID=="s1")
      s2 <- dplyr::filter(baseline.tbl,segID=="s2 ")
      s3 <- dplyr::filter(baseline.tbl,segID=="s3")
      
      Q <- 2*pi*s2$fs*(s2$L1/s2$R1)
      RMS.deltaV <- s2$`Ch1 (RMS)` - s2$`Ch2 (RMS)`
      Amp.vec <- amplitude.factor*Q*sqrt(2)*RMS.deltaV/1000
      
      baseline.tbl <- baseline.tbl%>%
        mutate(dR = ifelse(segID=="s1",s2$R1-s1$R1,ifelse(segID=="s3",s2$R1-s3$R1,0)))%>%
        mutate(dFs = ifelse(segID=="s1",s2$fs-s1$fs,ifelse(segID=="s3",s2$fs-s3$fs,0)))%>%
        filter(segType == "B")%>%
        mutate(Amp = ifelse(segID=="s1",Amp.vec,Amp.vec))%>%
        mutate(Baseline.Type = ifelse(segID=="s1","Baseline Before","Baseline After"))%>%
        mutate(ID=i)%>%
        mutate(rootPower = sqrt(Power))
      
      data.lst[[i]] <- baseline.tbl
    }
    full.data <- bind_rows(data.lst)
    final.data <- full_join(filename.tbl(),full.data,by="ID")
  })
  
  #Frequency baseline comparison vs. sqrt power plot 
  output$freq <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x=rootPower,y=fs,color=Test.ID,shape=Baseline.Type))+
      labs(title = "Frequency Baseline Comparison",
           x = expression(paste(sqrt(Power), " (",mu,"W)")),
           y = "f (Hz)")+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #Resistance baseline comparison vs. sqrt power plot
  output$resistance <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x=rootPower,y=R1,color=Test.ID,shape=Baseline.Type))+
      labs(title = "Resistance Baseline Comparison",
           x = expression(paste(sqrt(Power), " (",mu,"W)")),
           y = expression(paste("R (",Omega,")")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #dR vs. amplitude plot
  output$dR <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x=Amp,y=dR,color=Test.ID,shape=Baseline.Type))+
      labs(title = "dR Comparison",
           x = paste("Amplitude", "(nm)"),
           y = expression(paste("dR (",Omega,")")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  #dFs vs. amplitude plot
  output$dFs <- renderPlot({
    ggplot(tidy.data.long())+
      geom_point(aes(x=Amp,y=dFs,color=Test.ID,shape=Baseline.Type))+
      labs(title = "dF Comparison",
           x = paste("Amplitude", "(nm)"),
           y = expression(paste("dF"," (Hz)")))+
      #edit text sizes for title and axes labels
      theme(plot.title = element_text(size=24,face="bold"),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
  })
  
  # output$table <- renderTable(
  #   tidy.data.long()
  # )
  
}

shinyApp(ui, server)

