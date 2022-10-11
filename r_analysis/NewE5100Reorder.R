#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(stringr)
library(tidyverse)
library(readr)

# Define UI for application 
ui <- fluidPage(
  #Choose file button
  fluidRow(
    fileInput("upload", NULL, buttonLabel = "Choose File:", multiple = FALSE)
  ),
  #Shows preview of reordered data
  fluidRow(
    column(width=12,
           tableOutput("ShowData"))
  ),
  #Download data button
  fluidRow(
    downloadButton("reorderData", "Download reordered data")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #Code for reading single file input
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
  
  #Create separate datasets for baseline before, during test, and baseline after 
  #And sort by power using arrange()
  #Select only the necessary variables using select()
  
  #filter segment 1 data based on segID
  seg1 <- reactive({ tidy.tbl2()%>%
      #This filter will eventually be automated based on segType == R,T,B from labview
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
  
  #Use cbind() to merge the 3 data segments side-by-side
  reordered.data <- reactive({
    cbind(seg1(),seg2(),seg3())
  })
  
  #Show a preview of the reordered data
  output$ShowData <- renderTable({
    head(reordered.data(), 5)
  })
  
  #Code for dowmloading reordered data
  output$reorderData <- downloadHandler(
    filename = function() {paste("Reordered ", input$upload$name, ".csv", sep = "")
    },
    content = function(file){
      write.csv(reordered.data(), file, row.names = FALSE)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)


