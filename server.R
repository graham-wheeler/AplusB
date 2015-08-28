# Here is server.R

library(shiny)
library(shinyjs)
library(knitr)
library(bcrm)
library(ggplot2)
library(mvtnorm)
library(grid)
library(plyr)
library(boot)
library(xtable)
library(parallel)
source("helpers.R")

appCSS <- "
#loading-content {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}
"

shinyServer(function(input, output) {
  
  output$activeTab <- reactive({
    return(input$tab)
  })
  outputOptions(output, 'activeTab', suspendWhenHidden=FALSE)
  
  output$inumdoses <- renderUI({
    sliders<- list()
    sliders[[1]]<-sliderInput(paste0('slider', 1), paste('True probability of toxicity at dose level ', 1), min=0, max=1, value=0.1, step=0.01)
    for(i in 2:input$numdoses){
      sliders[[i]]<-sliderInput(paste0('slider', i), paste('True probability of toxicity at dose level ', i), min=0, max=1, value=min(0.99,i*0.1), step=0.01)
    }
    paste_all = function(...) paste(..., collapse = '\n')
    HTML(do.call('paste_all', sliders))
  })
  
  output$C <- renderUI({
    x<-input$A
    sliderInput("C", "C (Minimum number of DLTs needed out of 'A' patients in order to assign 'B' more)", min=0, max=x, value=1, step=0)
    })
  
  output$D <- renderUI({
    y<-input$A
    z<-input$C
  sliderInput("D", "D (Maximum number of DLTs needed out of 'A' patients in order to assign 'B' more - otherwise, trial stops/de-escalates)", min=z, max=y, value=z, step=1)
  })
  
  output$E <- renderUI({
    y<-input$D
    z1<-input$A
    z2<-input$B
    sliderInput("E", "E (Maximum number of DLTs permitted out of 'A+B' patients so trial may proceed/not de-escalate)", min=y, max=z1+z2, value=y, step=1)
  })
  
  
  dataInput <- reactive({
    beginning <- Sys.time()
    if(input$checkbox==FALSE){foo<-threep3ABCDE(sapply(1:input$numdoses, function(j) eval(parse(text=paste("input$slider",j,sep="")))), A=input$A, B=input$B, C=input$C, D=input$D, E=input$E)
    }else{
      foo<-threep3ABCDE.desc(sapply(1:input$numdoses, function(j) eval(parse(text=paste("input$slider",j,sep="")))), A=input$A, B=input$B, C=input$C, D=input$D, E=input$E)
    }
    end <- Sys.time()
    print(end - beginning)
    foo
  })
  
  data<-eventReactive(input$actionButton, {
    dataInput()
  })
  
  output$threep3Plot <- renderPlot({
    plot.threep3(data())
  })
  
  output$text1 <- renderText({
    if(input$actionButton >0){
      "The plots displayed here are specific to the choice of design parameters as well as scenario specific parameters (number of dose levels and true dose-toxicity probabilities)."
      }else{"The plots displayed here are specific to the choice of design parameters as well as scenario specific parameters (number of dose levels and true dose-toxicity probabilities). To generate these, click ''Get Design Properties!''."}
  })
  
  output$text2 <- renderText({
    if(input$actionButton>0){"The operating characteristics displayed here are specific to the choice of design parameters as well as scenario specific parameters (number of dose levels and true dose-toxicity probabilities)."}else{"The operating characteristics displayed here are specific to the choice of design parameters as well as scenario specific parameters (number of dose levels and true dose-toxicity probabilities). To generate these, click ''Get Design Properties!''."}
  })
  
  output$plottext <- renderText({
       if(input$actionButton >0){
    "Figure of operating characteristics averaged over all possible trials: distribution of sample size [top left]; probability of experimentation at each dose level [top right]; probability of recommending each dose level as the MTD [bottom left]; distribution of DLT rates [bottom right]."}
  })
  
  output$tabtext1 <- renderText({
    if(input$actionButton>0){
      "Table 1: Sample size distribution summary."
    }
  })
  output$tabtext2 <- renderText({
    if(input$actionButton>0){
      "Table 2: Experimentation and MTD recommendation proportions per dose level."
    }
  })
  output$tabtext3 <- renderText({
    if(input$actionButton>0){
      "Table 3: Experimentation and MTD recommendation proportions relating to true probabilities of DLT."
    }
  })
  output$tabtext4 <- renderText({
    if(input$actionButton>0){
      "Table 4: Expected Toxicity Level (ETL - expected DLT rate at the MTD averaged over all trials), mean number of DLTs per trial and Expected Overall Toxicity Rate (EOTR - mean number of DLTs divided by the mean sample size)."
    }
  })
  output$tabtext5 <- renderText({
    if(input$actionButton>0){
      paste("Table 5: ", input$w, "% Clopper-Pearson confidence intervals for all possible data collections at the chosen MTD.", sep="")
    }
  })
  output$tabtext6 <- renderText({
    if(input$actionButton>0){
      "This is the true probability of DLT at which, given the values of A, B, C, D and E chosen, the chance of escalating to the next dose is equal to either stopping the trial or de-escalating."
    }
  })
  
  output$tabtext <- renderText({
    if(input$actionButton>0){
      "* Of all trials that recommend an MTD"
    }
  })

  x1<-reactive(print.threep3(data())$tab0)
  x2<-reactive(print.threep3(data())$tab)
  x3<-reactive(print.threep3(data())$tab2)
  x4<-reactive(stats.exp.fn(data())$tab2)
  
  output$threep3Tab0 <-renderTable({
    x1()
  })
  
  output$threep3Tab1 <-renderTable({
   x2()  
  })
  
  output$threep3Tab2 <-renderTable({
    x3()
  })
  
  output$ETLEOTR <-renderTable({
    x4()
  })
  
  output$ClopperPearson <- renderTable({
    if(input$actionButton >0){
    clopperPearsonAplusB(input$w,input$A,input$B,input$C,input$E,deesc=input$checkbox)}
  }, include.rownames=FALSE)
  
  output$TippingPoint <- renderText({
    if(input$actionButton >0){
      x<-round(foofind(input$A,input$B,input$C,input$D,input$E),3)
      paste("Tipping Point = ",x,".", sep="")}
  })
  
output$report = downloadHandler(
filename = function() { paste(input$A,input$B,input$C,input$D,input$E,input$checkbox,input$numdoses,Sys.time(), '.pdf', sep='') },
content = function(file) {
out = knit2pdf('input.Rnw', clean = TRUE)
file.copy(out, file) # move pdf to file for downloading
  },
  contentType = 'application/pdf'
  )

})
