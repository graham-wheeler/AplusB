# Here is ui.R

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
library(binom)
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
  color: #36CE85;

}
"

shinyUI(fluidPage(

  div(
    id = "app-content",

  # Application title
  titlePanel("AplusB: A + B design investigator for phase I dose-escalation studies"),
  h4("Graham Wheeler"),
  h5("Cancer Research UK & UCL Cancer Trials Centre, University College London, UK;", a("graham.wheeler@ucl.ac.uk", href="mailto:graham.wheeler@ucl.ac.uk")),
  h5("MRC Biostatistics Unit Hub for Trials Methodology Research, Cambridge, UK;", a("graham.wheeler@mrc-bsu.cam.ac.uk", href="mailto:graham.wheeler@mrc-bsu.cam.ac.uk")),
  
  
  # Sidebar with a slider input for the number of bins

  sidebarLayout(
    sidebarPanel(h4("Scenario parameters"),
                 tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #36CE85;
               z-index: 105;
             }
          ")),
      sliderInput("numdoses",
                  "Number of dose levels:",
                  min = 1,
                  max = 10,
                  value = 4),
      uiOutput("inumdoses"),
      br(),
      h4("Design parameters"),
      sliderInput("A",
                  "A (Number of patients in first cohort that receive a certain dose)",
                  min = 1,
                  max = 6,
                  value = 3),
      sliderInput("B",
                  "B (Number of patients in second cohort (if required) that receive a certain dose)",
                  min = 1,
                  max = 6,
                  value = 3),
      uiOutput("C"),
      uiOutput("D"),
      uiOutput("E"),
      column(12, checkboxInput("checkbox", label = "Allow dose de-escalation?", value = FALSE)),
      sliderInput("w", "Confidence level for confidence intervals around identified MTD (%)", min=50, max=99, value=95),
      actionButton("actionButton", "Get design properties"),
      downloadButton('report', 'Download report (PDF)'),
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("Loading...",id="loadmessage"))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      strong("This application provides exact operating characteristics for a phase I dose-escalation study conducted using the rule-based A+B design [1]."),
      p("The operating characteristics presented are determined using the probability of every possible A+B trial occurring so that approximations are avoided [2, 3]."),
      p("Click",a("here",target="_blank",href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0159026"),"to view the paper for this application [2], which includes details on the formulae used to generate the ouputs below (see Supporting Information)."),
      p("Click",a("here", target="_blank", href="https://github.com/graham-wheeler/AplusB"),"to download the AplusB application to your computer from GitHub (set-up instructions provided)."),
      strong("References"),
      p("[1] Lin, Y. and Shih, W.J. (2001). Statistical properties of the traditional algorithm-based designs for phase I cancer clinical trials, ", em("Biostatistics"), strong("2"),"2:203-215."),
      p("[2] Wheeler, G.M., Sweeting, M.J. and Mander, A.P. (2016). AplusB: a web application for investigating A + B designs for phase I cancer clinical trials, ", em("PLoS ONE"), strong("11"), "7:e0159026. doi:10.1371/journal.pone.0159026"),
      p("[3] Wheeler, G.M. (2014). Adaptive designs for phase I dose-escalation studies, ", em("PhD thesis, University of Cambridge.")),
      tabsetPanel(
        tabPanel("Scenario plots",
                 tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #36CE85;
               z-index: 105;
             }
          ")),
      textOutput("text1"),
      plotOutput("threep3Plot"),
      textOutput("plottext")
    ),
    tabPanel("Scenario operating characteristics",
             tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #36CE85;
               z-index: 105;
             }
          ")),
            textOutput("text2"),
      tableOutput("threep3Tab0"),
      textOutput("tabtext1"),
      br(),
      tableOutput("threep3Tab1"),
      textOutput("tabtext2"),
      br(),
      tableOutput("threep3Tab2"),
      textOutput("tabtext"),
      textOutput("tabtext3"),
      br(),
      tableOutput("ETLEOTR"),
      textOutput("tabtext4"),
      br()
    ),
    tabPanel("Design operating characteristics",
             tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #36CE85;
               z-index: 105;
             }
          ")),
        p("The operating characteristics displayed here depend on the design only, i.e. values of A, B, C, D, E, the option to include dose de-escalation and the confidence level for the confidence interval calculations. These values will not change when the number of doses or the dose-toxicity probabilities are changed."),
      tableOutput("exact"),
      textOutput("tabtext5a"),
      br(),
      tableOutput("wilson"),
      textOutput("tabtext5b"),
      br(),
      textOutput("TippingPoint"),
      textOutput("tabtext6"),
      br(),
      br()
    ),id="tab"
    )
    )
  )  
  )
))
