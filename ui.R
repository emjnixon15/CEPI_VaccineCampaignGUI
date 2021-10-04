#library(shinymanager)

##======================================##
## Author: Anita Lerch &  Alex Perkins
##======================================##

homeParamTab <- function(){
  tabPanel(title="Background", value="backTab",
           br(),
           p("This project was implemented by", strong("Perkins Lab at University of Notre Dame"), "and funded by the", strong("Coalition for Epidemic Preparedness Innovations (CEPI).")),
           p("CEPI is currently investing in the development of a number of vaccine candidates and platforms, as well as undertaking planning exercises around the manufacturing and future deployment of vaccines once they are approved for use. This modeling project was launched with the aim of informing CEPI’s", strong("Sustainable Manufacturing Strategy by"), "providing estimates of vaccine demand for stockpiling and outbreak response."),
           p(strong("Aim:"), "To inform CEPI’s Sustainable Manufacturing strategy by undertaking modeling exercises to understand the magnitude of vaccine manufacturing and stockpiling needed for CEPI target pathogen vaccines (Lassa, MERS, Nipah, Rift Valley Fever, chikungunya)."),
           p(strong("Objectives:"),
             tags$ul(
               tags$li("To understand estimates for vaccine manufacturing for a) stockpile needs and b) outbreak response for CEPI target pathogens"),
               tags$li("To estimate vaccine demand curve needs across a number of different outbreak scenarios and different vaccination responses"),
               tags$li("To develop a modeling approach that can be modified for use in estimating vaccine demand for other pathogens")
             )
           ),
           # [Add text on process, modeling approach, methodology, other background]
           hr(),
           # [Add text on Notre Dame, research group and/or individual researchers?]
           p("CEPI is an innovative global partnership between public, private, philanthropic, and civil society organisations launched in Davos in 2017 to develop vaccines to stop future epidemics. CEPI’s mission is to accelerate the development of vaccines against emerging infectious diseases and enable equitable access to these vaccines for people during outbreaks.")
           
  )
}

epiParamTab <- function(){
  tabPanel(title="Epidemiology *", value="epiTab", 
           br(),
           p(strong("Instructions:"),"Adjust the epidemiological parameters. The black solid line in the figures shows the distribution for the default values and the red solid line shows the distribution for the modified parameter values. The corresponding mean values are represented by dashed vertical lines."),
           hr(),
           h4("Yearly Spillover Size:"),
           fluidRow(
             column(4,
                    numericInput(inputId="spo_mean", "Mean", min=0.0, max=5000, step=10, value=NA),
                    numericInput(inputId="spo_sd", label="SD", min=0, step=0.1, value=NA),
                    numericInput(inputId="spo_z", label="Zero probability (z)", min=0, max=1, step=0.01, value=0)
                    ),
             column(4, uiOutput("SpillOverSize")),
             column(4, p("Distribution of the number of spillovers per year globally. ‘Mean’ represents the average number of spillovers observed per year. ‘SD’ represents the standard deviation of the numbers of spillovers between different years. ‘z’ represents the probability of years with no observed and unobserved spillovers. Distribution of spillover per year is modeled using a zero-inflated negative binomial distribution with this mean, standard deviation and probability of zero spillover years."))
             ),
           hr(),
           h4("Spillover Probability of Catchment Areas:"),
           fluidRow(
             column(4,
                    numericInput(inputId="clu_mu", "Mean", step = 0.01, value=NA),
                    numericInput(inputId="clu_sd", "SD", step = 0.01, value=NA)
             ),
             column(4, uiOutput("SpillOverProb")),
             column(4, p("Distribution of spillover probability across catchment areas. ‘Mean’ represents the average spillover probability over all catchment areas. ‘SD’ represents the standard deviation of the spillover probability across catchment areas. Distribution of spillover probability is modeled using a beta distribution with this mean and standard deviation."))
             ),
           hr(),
           h4("Seasonal Timing of Spillovers:"),
           fluidRow(
             column(4,
                    numericInput(inputId="sea_mu", "Mean", step = 1, min=1, max=365, value=NA),
                    numericInput(inputId="sea_sd", "SD", step = 1, min=1, max=365, value=NA)
             ),
             column(4, uiOutput("Timing")),
             column(4, p("Seasonal distribution of spillovers. ‘Mean’ represents the average week on which spillovers are observed. ‘SD’ represents the standard deviation of the weeks across which spillovers are observed. Seasonal distribution of spillover is modeled using a beta distribution with this mean and standard deviation."))
           ),
           hr(),
           h4("Basic Reproduction Number:"),
           fluidRow(
             column(4,
                    numericInput(inputId="R0_mean", "Mean", min=0.0, max=0.99, step = 0.01, value=NA),
                    checkboxInput(inputId="R0_nbDist", label="Negative Binomial distribution?", value=FALSE),
                    conditionalPanel(condition="input.R0_nbDist",
                                     numericInput(inputId="R0_sd", label="SD", min=0, step=0.01, value=NA)
                                     )
                    ),
             column(4, uiOutput("R0")),
             column(4, p("Distribution of the number of secondary cases per primary case. ‘Mean’ represents the average number of secondary cases per primary case. ‘SD’ represents the standard deviation of the number of secondary cases per primary case. Distribution of secondary cases per primary case is modeled using a Poisson or negative binomial distribution with this mean and standard deviation."))
             ),
           hr(),
           h4("Incubation Period:"),
           fluidRow(
             column(4,
                    numericInput(inputId="inc_mu", "Mean", step = 0.1, value=NA),
                    numericInput(inputId="inc_sd", "SD", step = 0.1, value=NA)
                    ),
             column(4, uiOutput("IncubationTime")),
             column(4, p("Distribution of number of days until a person becomes infectious after being infected. ‘Mean’ represents the average number of days until a person becomes infectious. ‘SD’ represents the variation in number of days a person becomes infectious. Distribution of incubation period is modeled using a gamma distribution."))
           ),
           hr(),
           h4("Infectious Period:"),
           fluidRow(
             column(4,
                    numericInput(inputId="rec_mu", "Mean", step = 0.1, value=NA),
                    numericInput(inputId="rec_sd", "SD", step = 0.1, value=NA)
             ),
             column(4, uiOutput("RecoveryTime")),
             column(4, p("Distribution of the number of days a person is infectious. ‘Mean’ represents the average number of days a person is infectious. ‘SD’ represents the standard deviation of the numbers of days a person is infectious. Distribution of infectious period is modeled using a gamma distribution with this mean and standard deviation."))
           ),
           hr()
  )
}

demoParamTab <- function(){
  tabPanel("Demographics *", value="demoTab", #h4("Catchment Area Parametes")
           br(),
           p(strong("Instructions:"),"Adjust the demographic parameters. The black solid line in the figure shows the distribution for the default values and the red solid line shows the distribution for the modified parameter values. The corresponding mean values are represented by dashed vertical lines."),
           hr(),
           h4("Total population:"),
           p("The estimated total population across all catchment areas with the current parameters is: ", textOutput("globalPopEst", inline=T)),
           hr(),
           h4("Type of Catchment Areas"),
           p("The type of catchment areas is implicitly represented by the number of catchment areas and population size. Simulations of admin 1 are represented by a smaller number of catchment areas but with larger average population size. Hospital, health care facilities and admin 2 are represented by a higher number of catchment areas but with smaller average population size.  Overall the estimated total population across all catchment areas should sum up to the same value irrespective of the type of catchment areas."),
           hr(),
           h4("Number of Catchment Areas:"),
           fluidRow(
             column(4,
                    numericInput(inputId="numCluster", label=NULL, min=1, step=100, value=100)
             ),
             column(4, p("Default:"),
                    p(textOutput("defaultTypeCluster")),
                    p(textOutput("defaultNumCluster"))),
             column(4, p("The 'Number of catchment areas' defines the number of catchment areas used to distribute the yearly spillovers. These are considered to be catchment areas globally where this is a risk of spillover of the pathogen."))
           ),
           hr(),
           h4("Population Size in Catchment Area:"),
           fluidRow(
             column(4,
                    numericInput(inputId="cluster_mean", "Mean", min=0.0, step = 10, value=NA),
                    numericInput(inputId="cluster_sd", label="SD", min=0, step=0.1, value=NA)
             ),
             column(4, uiOutput("Population")),
             column(4, p("The 'Clusters Population' distribution is used to estimate the number of regimens needed based on the population size in the catchment areas. ‘Mean’ represents the average population size in a catchment area. ‘SD’ represents the standard deviation of the population sizes across catchment areas. Distribution of population size is modeled using a beta distribution with this mean and standard deviation."))
           )
  )
}

otherParamTab <- function(){
  tabPanel("Simulation *", value="simTab",
           h4('Number of Replicates'),
           fluidRow(
             column(4, numericInput(inputId="replicats", label=NULL, min=1, max=1000, step=10, value=100)),
             column(8, p("The number of replicates used to run the simulation. Fewer replicates run faster but with more variability in the results. More replicates run slower but with more precise estimation of the number of regimens required. Use a smaller number of replicates for testing different parameters and larger numbers of replicates – e.g. 500 – for final estimation."))
           ),
           hr(),
           h4('Maximum Allowed Outbreak Size'),
           fluidRow(
             column(4, numericInput(inputId="maxOutbreakSize", label=NULL, min=100, step=10, value=5000)),
             column(8, p("The maximum outbreak size allowed per catchment area."))
           )  
  )
}

campaignParamTab <- function(){
  tabPanel("Vaccine Campaign *", value="vacTab",
           h4("Regimens by Contact:"),
           fluidRow(
             column(4, numericInput(inputId="numContact", label="Contacts per Case", min=1, step=1, value=90)),
             column(8, p("The 'Regimens by Contact' defines the average number of contacts of a case. It is used to estimate the number of regimens needed based on contacts under ring vaccination."))
           ),
           hr(),
           h4("Campaign Start Threshold:"),
           fluidRow(
             column(4, numericInput(inputId="tho_cases", label="Cases", min=1, step=1, value=NA),
                    numericInput(inputId="tho_days", label="Days", min=1, step=1, value=NA)),
             column(8, p("The campaign start threshold defines how many ‘Cases’ in a sliding window of lengths ‘Days’ are required to trigger an outbreak response to start a campaign."))
           ),
           hr(),
           h4("Vaccination:"),
           fluidRow(
             column(4, numericInput(inputId="vaccov_pop", label="Immunization coverage", min=0, max=1, step=0.05, value=NA)),
             column(8, p("The proportion of the population that is vaccinated."))
           ),
           #hr(),
           #h4("Vaccination:"),
           #tags$label("for"="vaceff_2dose", "Numbers of regimens required", class="input-label"), # label
           fluidRow(
             column(4, checkboxInput(inputId="vaceff_2dose", label="Two-dose vaccination?", value = FALSE)),
             column(8, p("The number of regimens required to achieve full protection."))
           ),
           fluidRow(
             column(4, 
                    numericInput(inputId="vaceff_dose1", label="Efficacy 1st dose", min=0, max=1, step=0.05, value=NA),
                    conditionalPanel(condition="input.vaceff_2dose",
                                     numericInput(inputId="vaceff_dose2", label="Efficacy 2nd dose", min=0, max=1, step=0.05, value=NA)
                    )
             ),
             column(8, p("The proportion of vaccine recipients achieving protection after the first and second dose, respectively."))
             ),
           hr(),
           h4("Administration Delay:"),
           fluidRow(
             column(4,            
                    numericInput(inputId="vacdel_dose1", label="1st dose", min=1, step=1, value=NA),
                    conditionalPanel(condition="input.vaceff_2dose",
                                     numericInput(inputId="vacdel_dose2", label="2nd dose", min=1, step=1, value=NA)
                    )
             ),
             column(8, p("Delay in days to administer the overall vaccination regimen. The delay for the first dose is the time from campaign start to administration of the first dose. For the second dose, it is the time between administration of the first and second dose."))
           ),
           hr(),
           h4("Protection delay:"),
           fluidRow(
             column(4, numericInput(inputId="vacprot", label="Days", min=1, step=1, value=7)),
             column(8, p("The protection delay is the time needed for the immune system to develop maximal protection."))
           )
  )
}

resultTab <- function(){
  tabPanel(title="Vaccine Estimates", value="resTab",
           h4("Campaign Simulation Estimates"),
           #downloadButton("downloadData", "Download Estimations"),
           span(strong(textOutput("warnIncomplete")), style="color:red"),
           p("Number of cases averted and regimens required to address outbreak response. Because each replicate simulation required a different number of regimens, the reported numbers represent moments of the distribution of regimens required across those replicates."),
           uiOutput("CampaignStats"),
           hr(),
           h4('Total number of spillover cases across catchment areas'),
           fluidRow(
             column(6, uiOutput("DistSpillover")),
             column(6, p("Simulated spillover size per year across replicates (gray bars) and yearly spillover size distribution (red curve)."))
           ),
           hr(),
           h4("Average campaigns per week"),
           fluidRow(
             column(6, uiOutput("DistCampaignStart")),
             column(6, p("Seasonal timing of spillovers (red curve) and outbreak response campaigns (gray bars). The height of the bars represents the average number of outbreak response campaigns triggered in a given week."))
           ),
           hr(),
           h4("Seasonal timing of cases"),
           fluidRow(
             column(6, uiOutput("TimmingCases")),
             column(6, p("Fitted seasonal timing of spillovers (red curve) and simulated incidence of cases deriving from spillover (dark gray) and human-to-human transmission (light gray). The height of the bars represents the average  number of cases in a given week."))
           )
  )
}

ui = fluidPage(
    tags$style(type='text/css', ".control-label { font-size: 11px; }"),
    img(src='cepi_logo.jpg', align="right"),
    h1("CEPI Vaccine Campaign Simulator"),
    p(strong("Instructions:"),"Choose disease, adjust default parameters in tabs marked with * and press the GO button to see the results in the 'Vaccine Estimates' tab."),
    fluidRow(
      column(4, 
             #tags$label("for"="disease", "ddisease", class="input-label"), # label
             selectInput(inputId="disease", label=NULL, c("CHOOSE DISEASE"="", 
                                                "Lassa"="lassa", 
                                                "MERS"="mers", 
                                                "Nipah"="niv", 
                                                "RVF" ="rvf"), width="100%", selectize=FALSE)
             ),
      column(8,
             conditionalPanel(condition="input.disease!=''",
                              actionButton(inputId="goButton", label="GO!", width="100%", icon=icon("refresh"),
                                           style="color: #fff; background-color: #007bff; border-color: #fff")) # #337ab7 2e6da4
             )),

    conditionalPanel(condition="input.disease==''", {
      tabsetPanel(id="introTab", type="tabs",
                  homeParamTab()
      )
    }),
    conditionalPanel(condition="input.disease!=''", {
           tabsetPanel(id="paramTab", type="tabs",
                homeParamTab(),
                epiParamTab(),
                demoParamTab(),
                campaignParamTab(),
                otherParamTab(),
                resultTab()
            )
      })
)

# ui <- secure_app(ui)
