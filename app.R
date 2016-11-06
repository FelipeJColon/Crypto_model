## ############################################################################
##
## DISCLAIMER: 
##
## This script has been developed for illustrative purposes only. 
## The script is provided without any warranty of any kind, either express or 
## implied. The entire risk arising out of the use or performance of the sample 
## script and documentation remains with you. In no event shall its
## author, or anyone else involved in the creation, production, or delivery of 
## the script be liable for any damages whatsoever (including, without 
## limitation, damages for loss of business profits, business interruption, 
## loss of business information, or other pecuniary loss) arising out of the use 
## of or inability to use the sample scripts or documentation, even if the 
## author has been advised of the possibility of such damages. 
##
## ############################################################################
##
## DESCRIPTION
## Shiny app to simulate epidemic curves for point source Cryptosporidium 
## outbreaks
##
## Version 1: Initially created on 2 Nov 2016
##
## Dependencies: None
##
## Written by: Felipe J Colón-González
## For any problems with this code, please contact f.colon@uea.ac.uk
## 
## ############################################################################

## -------------------------
## Source packages
## -------------------------

require(shiny)
require(data.table)
require(ggplot2)
require(deSolve)
require(tidyr)
require(plyr)
require(DT)
require(BH)


## -------------------------
## Define UI options
## -------------------------

ui <- shinyUI(fluidPage(
       
       #  Application title
       headerPanel(
              h1(em("Cryptosporidium"), "Outbreak Simulator")),
       
       ##--- Initial conditions
       
       # Sidebar with available options for initial conditions
       tabsetPanel(type="tabs",
                   
                   # Tab 1
                   tabPanel("Documentation", br(), 
                            h3("Welcome to the online", em("Cryptosporidium"),
                               "outbreak simulator!"),
                            
                            h4("What is this App about?"),
                            
                            p("This Shiny App has been developed", strong("to simulate
                              outbreaks"), "of a disease called cryptosporidiosis
                              which is caused by a microscopic  parasite called",
                              em("Cryptosporidium spp."), "Cryptosporidiosis is
                              the leading cause of human waterborne disease in 
                              many developed countries including the UK and the
                              USA. The most common 
                              transmission pathway is through contaminated water 
                              (e.g. dringking or recreational). 
                              The disease may affect all people, however, people with
                              compromised immune systems (e.g. those with AIDS), are
                              more at risk of developing severe symptoms. The most 
                              frequent symptom is watery diarrhoea but may also include
                              nausea, vomiting, abdominal cramps, and fever. Symptoms
                              may last between 1 and 2 weeks.",
                              
                              
                              br(),br(), # breaks
                              strong("This App allows you to predict the number of
                              people that could become, infected and
                              symptomatic"), "if oocysts (i.e. the immature forms
                              of the parasite that infect people) are deliberately 
                              or accidentally released into a water source.",
                              strong("The main goal of the App is to serve as a 
                                    support tool for public health decision-making.")),
                            
                            br(),
                            h4("Who should use the App?"),
                            
                            p("The App is intended for public health professionals
                              either in the decision-making or academic fields. 
                              However, the App could be easily used by anyone."),
                            
                            br(),
                            h4("How can I use the App?"),
                            
                            p("The App has of three main components. First, 
                              you will find a", strong("SIDE PANEL"), "where you 
                              will be able 
                              to tell the App the characteristics of the population
                              where you would like to model the outbreak.
                              You can easily do that by changing the values in 
                              the boxes. For 
                              example, you may indicate the number of people that
                              were infected with the disease before the outbreak
                              started and the amount of un-boiled water people
                              consume, on average, every day. Second, you have 
                              a", strong("LOWER PANEL"), "where you may tell the 
                              App the 
                              characteristics of the outbreak you would like to
                              model using the slider. For example, you may indicate 
                              the proportion
                              of the population that has access to contaminated 
                              water, or how many days does it take for the symptoms
                              to appear. Finally, you have a", strong("MAIN PANEL"),
                              "where 
                              the output of the computations are shown in a 
                              plot. The output plot shows the predicted number of
                              people in each category (i.e. symptomatic, asymptomatic,
                              recovered...) on a given day for a period that you can
                              define (in days). You may notice that you can download 
                              the output by clicking on the button below the graph."),
                            
                            br(),
                            h4("Before you use the App"),
                            
                            p("Make sure you have the following packages installed:
                              data.table, ggplot2, deSolve, tidyr, plyr, DT, and BH.
                              If you are not familiar with compartmental modelling
                              for infectious diseases I would strongly encourage
                              you to get familiar with some terminology before you
                              use this App.",
                              
                              br(), br(),
                              strong("- Compartmental model:"), " A series of differential
                            equations that are used to simulate the behaviour, in this
                            case, of a disease across time.",
                              
                              br(),
                              strong("- Initial infected and latent:"), "The number of
                            infected people that have not yet started to show 
                            symptoms.",
                              
                              br(),
                              strong("- Initial infectious and symptomatic:"), "The 
                            number of infected and infectious people with symptoms.",
                              
                              br(),
                              strong("- Initial infectious and asymptomatic:"), "The 
                            number of infected and infectious people that do not 
                            develop any symptoms.",
                              
                              br(),
                              strong("- Initial recovered:"), "The 
                            number of infected people that no longer suffer the
                            disease.",
                              
                              br(),
                              strong("- Days before public health action:"), "The 
                            number of days before the health authorities implement
                            a mechanism (intervention) to prevent further 
                            transmission.",
                              
                              br(),
                              strong("- Incubation period:"), "The time between the
                            acquisition of the parasite and the onset of symptoms.",
                              
                              br(),
                              strong("- Infectious period:"), "The time over which
                            a person can transmit the disease to a susceptible
                            host."),
                            
                            br(),
                            h5("Created by:"),
                            tags$a("Felipe J Colon-Gonzalez", 
                                   href="https://github.com/FelipeJColon")),
                   
                   # Tab 2
                   tabPanel("Model simulations", br(),
                            fluidRow(
                                   column(3, 
                                          wellPanel(
                                                 tags$h4("Initial conditions"),
                                                 
                                                 numericInput(inputId="N", label="Population:", value=1e6, 
                                                              step=1e4),
                                                 numericInput(inputId="L", 
                                                              label="Initial infected and latent:", 
                                                              value=0),
                                                 numericInput(inputId="I",
                                                              label="Initial infectious and symptomatic:", 
                                                              value=0),
                                                 numericInput(inputId="A", 
                                                              label="Initial infectious and asymptomatic:", 
                                                              value=0),
                                                 numericInput(inputId="R", 
                                                              label="Initial recovered:", value=0),
                                                 numericInput(inputId="intake", 
                                                              label="Average water intake (L/day):", 
                                                              value=1.8, min=1, max=2.4, step=0.1),
                                                 numericInput(inputId="ndays", label="Time span:", 
                                                              value=50),
                                                 numericInput(inputId="remed", 
                                                              label="Days before public health action:", 
                                                              value=3)
                                          )
                                   ),
                                   
                                   ##--- Output options 
                                   column(9,
                                          # Check box for plotting options
                                          checkboxGroupInput("Indicators", "",
                                                             c("Latent", 
                                                               "Symptomatic",
                                                               "Asymptomatic",
                                                               "Recovered"),
                                                             selected=c(
                                                                    "Symptomatic", 
                                                                    "Asymptomatic"),
                                                             inline=TRUE),
                                          # Output plot
                                          plotOutput("plot1"),
                                          tableOutput("datatable"),
                                          downloadButton('downloadData', 'Download data'))
                            ),
                            
                            ##--- Model parameters
                            
                            # Define model parameters as slidebars
                            fluidRow(
                                   column(3, 
                                          br(),
                                          h5("Created by:"),
                                          tags$a("Felipe J Colon-Gonzalez", 
                                                 href="https://github.com/FelipeJColon")),
                                   
                                   wellPanel(
                                          tags$h4("Model parameters"),
                                          
                                          column(3,
                                                 sliderInput(inputId="sigma", 
                                                             label="Incubation period (days):", 
                                                             min=1, max=21, value=5, step=0.1), 
                                                 sliderInput(inputId="gamma",
                                                             label="Infectious period (days):", 
                                                             min=3, max=50, value=7, step=0.1)), 
                                          column(3, 
                                                 sliderInput(inputId="p",  
                                                             label=paste("Fraction of Latent that become",
                                                                         "infectious and symptomatic:"), 
                                                             min=0.2, max=0.9, value=0.75, step=0.01),
                                                 sliderInput(inputId="oocysts",
                                                             label=paste("Number of oocysts released into",
                                                                         "water system (per litre)"),
                                                             min=10, max=1e6, value=100, step=5)),
                                          column(3, 
                                                 sliderInput(inputId="prop",
                                                             label=paste("Proportion of the population in",
                                                                         "contact with infected water",
                                                                         "(per day)"),
                                                             min=1e-6, max=1e-1, value=1e-2, step=1e-4)
                                                 
                                          ) ))
                   )
       )
))

## -------------------------
## Compartmental model 
## -------------------------

# Define the model structure
liarModel <- function (t, x, param) {
       
       # Set initial conditions and parameters as list
       with(as.list(c(x, param)),{    
              
              # Model differential equations
              dL <- -sigma * L 
              dI <- p * sigma * L - (gamma * I)
              dA <- (1 - p) * sigma * L - (gamma * A)
              dR <- gamma * (I + A)
              
              # Output
              list(c(dL, dI, dA, dR))
       }
       
       )
       
}

## -------------------------
## Define Server options
## -------------------------

# Simulation of Pandemic Influenza (H1N1) dynamics 
server <- shinyServer(function(input, output) {
       
       ##--- Reactive dataset
       mydata <- reactive({
              
              ##--- Model Parameters:
              
              # Define time span (from 0 to ndays with intervals of 1 day)
              times <- seq(0, input$ndays) 
              
              # Initial conditions
              init <- c(L=input$L, I=input$I, A=input$A, R=input$R)
              
              # Parameters for the hyper-geometric (i.e. beta-Poisson) 
              # dose-response function based on Teunis et al, 2002, 
              # Risk Anal 22(1): 175-183
              ialpha <- 0.115
              ibeta  <- 0.176
              
              # Per capita daily dose of oocysts based on an average daily consumption of 
              # un-boiled water
              dose <- input$oocysts * input$intake
              
              # Probability of infection (based on Teunis et al 2002)
              pinf <- round(1 - (1+(dose / ibeta))^ -ialpha, digits=1)
              
              # Define compartmental model parameters
              param <- c(p=input$p, sigma=1/input$sigma, gamma=1/input$gamma)
              
              # Define number of imported infectious cases per day
              prim_cases <- data.frame(var=c("L"), time=seq(0, input$remed),
                                       value=rpois(length(seq(0, input$remed)),
                                                   lambda=input$N*pinf*input$prop),
                                       method=c("add"))
              
              # Run the model and store output on data.table
              wide <- data.table(ode(y=init,times=times,
                                     fun=liarModel,parms=param,
                                     events=list(data=prim_cases)))
              
              # Re-format the output in long format
              long <- gather(wide, Indicator, Individuals, -time)
              
              # Rename indicators
              long$Indicator <- factor(long$Indicator)
              long$Indicator <- mapvalues(long$Indicator, from=c("L", "I", 
                                                                 "A", "R"), 
                                          to=c("Latent", "Symptomatic", 
                                               "Asymptomatic", "Recovered"))
              
              
              list(long=long, wide=wide)
       })
       
       ##--- Render plot
       output$plot1 <- renderPlot({
              
              # Get reactive dataset
              dataLong <- mydata()[["long"]]
              
              # Define plot structure
              myPlot <- ggplot(dataLong[dataLong$Indicator %in% input$Indicators,], 
                               aes(x=time, y=Individuals, group=Indicator)) + 
                     geom_line(aes(colour=Indicator), size=1.5, alpha=.75) + 
                     ggtitle("Population totals") +
                     theme_bw() + 
                     scale_x_continuous(name="Days") +
                     scale_y_continuous(name="Individuals") + 
                     theme(plot.title=element_text(size=22)) + 
                     theme(axis.title.y=element_text(size=rel(1.4),angle=90))+
                     theme(axis.title.x=element_text(size=rel(1.4),angle=00))+
                     theme(axis.text.x=element_text(angle=00,hjust=0.75,
                                                    size=rel(1.5),
                                                    color="black"))+
                     theme(axis.text.y=element_text(angle=00,hjust=0.75,
                                                    size=rel(1.2),
                                                    color="black")) +
                     theme(legend.position="right") +
                     theme(legend.title=element_blank())+
                     scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9",
                                                 "#009E73", "#F0E442", "#0072B2", 
                                                 "#D55E00", "#CC79A7"))
              
              # Print plot to screen
              print(myPlot)
       })
       
       output$downloadData <- downloadHandler(
              filename <- function() { paste0("simulated_data", ".csv") },
              content <- function(file) {
                     write.csv(mydata()[["wide"]], file, row.names = FALSE)
              }
       )
       
       #      output$datatable <- 
       #           renderTable({
       #                Tdata <- mydata()[["wide"]]
       #                Tdata <- cbind(day=1:nrow(Tdata), Tdata)
       #                Tdata[seq(1, nrow(Tdata), length.out=50),]
       #           })
       
       
})


shinyApp(ui=ui, server=server)