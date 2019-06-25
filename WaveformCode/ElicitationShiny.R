
{
  {
    if(file.exists('CheckforDefaultsScript.R')){
      source('CheckforDefaultsScript.R')
    }else{
      pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
      source("LibrariesAndSettings.R" , print.eval  = TRUE )
      DP_LoadPatientIndex()
      DP_ChooseDataReps()
      FilestoProcess <- DP_ChooseECGstoProcess() 
      HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
    }
    listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
    set.seed( 1 )
  }
}

# Prelimiaries 
{
  t_observation = seq(0.25  , 10 , 0.005)
  
}

{ui <- navbarPage(title = 'Heart-Rhythm Models Elicitaion Tool',
                  tabsetPanel( navbarMenu(title = 'Rhythms',
                                          tabPanel('Home' , tags$h1('Model Description')),
                                          tabPanel(title = 'Regular', 
                                                   wellPanel(tags$h1('Regular')) ,
                                                   fluidRow(
                                                     column(4 , wellPanel(sliderInput(inputId = "Reg_pi1" , 
                                                                                      label = 'Choose a Number for \\pi_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 1 ,
                                                                                      value = 0) ,
                                                                          sliderInput(inputId = "Reg_pi2" , 
                                                                                      label = 'Choose a Number for \\pi_2' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 1 ,
                                                                                      value = 1) ,
                                                                          #sliderInput(inputId = "PAC_pi3" , 
                                                                          #            label = 'Choose a Number for \\pi_3' , 
                                                                          #           min = 1 ,
                                                                          #            step = 0.001,
                                                                          #            max = 2 ,
                                                                          #            value = 0.1) ,
                                                                          sliderInput(inputId = "Reg_mu1" , 
                                                                                      label = 'Choose a Number for \\mu_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 3 ,
                                                                                      value = 0.75) ,
                                                                          sliderInput(inputId = "Reg_mu2" , 
                                                                                      label = 'Choose a Number for \\mu_2 - \\mu_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.25 ,
                                                                                      value = 0.05),
                                                                          sliderInput(inputId = "Reg_mu3" , 
                                                                                      label = 'Choose a Number for \\mu_3 - \\mu_2' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.25 ,
                                                                                      value = 0.05),
                                                                          sliderInput(inputId = "Reg_sigma1" , 
                                                                                      label = 'Choose a Number for \\sigma_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.1 ,
                                                                                      value = 0.01),
                                                                          sliderInput(inputId = "Reg_sigma2" , 
                                                                                      label = 'Choose a Number for \\sigma_2' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.1 ,
                                                                                      value = 0.01 ),
                                                                          sliderInput(inputId = "Reg_sigma3" , 
                                                                                      label = 'Choose a Number for \\sigma_3' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.1 ,
                                                                                      value = 0.01 )
                                                     ) ),
                                                     column(4 , plotOutput( outputId = 'Reg_RRDenisty' , inline = T)),
                                                     column(4 , plotOutput( outputId = 'Reg_RRTimes', inline = T))
                                                   ),
                                                   tags$h2('Simulated ECG'),
                                                   tags$hr(),
                                                   fluidRow(actionButton(inputId = 'Reg_UpdateECG' , label = 'Update')),
                                                   fluidRow(plotOutput( outputId = 'Reg_ECG', inline = F))
                                          ),
                                          tabPanel(title = 'Bigeminy', 
                                                   wellPanel(tags$h1('Bigeminy')) ,
                                                   fluidRow(
                                                     column(4 , wellPanel(sliderInput(inputId = "mu1" , 
                                                                                      label = 'Choose a Number for \\mu_2' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 2 ,
                                                                                      value = 0.5) ,
                                                                          sliderInput(inputId = "mu2" , 
                                                                                      label = 'Choose a Number for \\mu_2 - \\mu_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 2 ,
                                                                                      value = 0.2),
                                                                          sliderInput(inputId = "sigma1" , 
                                                                                      label = 'Choose a Number for \\sigma_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.1 ,
                                                                                      value = 0.01),
                                                                          sliderInput(inputId = "sigma2" , 
                                                                                      label = 'Choose a Number for \\sigma_2 - \\sigma1' , 
                                                                                      min = -0.05 ,
                                                                                      step = 0.001,
                                                                                      max = 0.05 ,
                                                                                      value = 0 ))),
                                                     column(4 , plotOutput( outputId = 'RRDenisty' , inline = T)),
                                                     column(4 , plotOutput( outputId = 'RRTimes', inline = T)
                                                     )
                                                   ),
                                                   tags$h2('Simulated ECG'),
                                                   tags$hr(),
                                                   fluidRow(actionButton(inputId = 'UpdateECG' , label = 'Update')),
                                                   fluidRow(plotOutput( outputId = 'ECG', inline = F))
                                          ),
                                          tabPanel(title = 'PACs', 
                                                   wellPanel(tags$h1('PACs')) ,
                                                   fluidRow(
                                                     column(4 , wellPanel(sliderInput(inputId = "PAC_pi1" , 
                                                                                      label = 'Choose a Number for \\pi_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 1 ,
                                                                                      value = 0.8) ,
                                                                          #sliderInput(inputId = "PAC_pi2" , 
                                                                          #            label = 'Choose a Number for \\pi_2' , 
                                                                          #            min = 0 ,
                                                                          #            step = 0.001,
                                                                          #            max = 1 ,
                                                                          #            value = 0.8) ,
                                                                          #sliderInput(inputId = "PAC_pi3" , 
                                                                          #            label = 'Choose a Number for \\pi_3' , 
                                                                          #           min = 1 ,
                                                                          #            step = 0.001,
                                                                          #            max = 2 ,
                                                                          #            value = 0.1) ,
                                                                          sliderInput(inputId = "PAC_mu1" , 
                                                                                      label = 'Choose a Number for \\mu_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 2 ,
                                                                                      value = 0.6) ,
                                                                          sliderInput(inputId = "PAC_mu2" , 
                                                                                      label = 'Choose a Number for \\mu_2 - \\mu_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 2 ,
                                                                                      value = 0.1),
                                                                          sliderInput(inputId = "PAC_mu3" , 
                                                                                      label = 'Choose a Number for \\mu_3 - \\mu_2' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 2 ,
                                                                                      value = 0.1),
                                                                          sliderInput(inputId = "PAC_sigma1" , 
                                                                                      label = 'Choose a Number for \\sigma_1' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.1 ,
                                                                                      value = 0.01),
                                                                          sliderInput(inputId = "PAC_sigma2" , 
                                                                                      label = 'Choose a Number for \\sigma_2' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.1 ,
                                                                                      value = 0.01 ),
                                                                          sliderInput(inputId = "PAC_sigma3" , 
                                                                                      label = 'Choose a Number for \\sigma_3' , 
                                                                                      min = 0 ,
                                                                                      step = 0.001,
                                                                                      max = 0.1 ,
                                                                                      value = 0.01 )
                                                     ) ),
                                                     column(4 , plotOutput( outputId = 'PACs_RRDenisty' , inline = T)),
                                                     column(4 , plotOutput( outputId = 'PACs_RRTimes', inline = T))
                                                   ),
                                                   tags$h2('Simulated ECG'),
                                                   tags$hr(),
                                                   fluidRow(actionButton(inputId = 'PACs_UpdateECG' , label = 'Update')),
                                                   fluidRow(plotOutput( outputId = 'PACs_ECG', inline = F))
                                          ),tabPanel(title = 'Atrial Fibrillation', 
                                                     wellPanel(tags$h1('Atrial Fibrillation')) ,
                                                     fluidRow(
                                                       column(4 , wellPanel(sliderInput(inputId = "AFib_pi1" , 
                                                                                        label = 'Choose a Number for \\pi_1' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 1 ,
                                                                                        value = 0.33333) ,
                                                                            sliderInput(inputId = "AFib_pi2" , 
                                                                                        label = 'Choose a Number for \\pi_2' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 1 ,
                                                                                        value = 0.33333) ,
                                                                            #sliderInput(inputId = "PAC_pi3" , 
                                                                            #            label = 'Choose a Number for \\pi_3' , 
                                                                            #           min = 1 ,
                                                                            #            step = 0.001,
                                                                            #            max = 2 ,
                                                                            #            value = 0.1) ,
                                                                            sliderInput(inputId = "AFib_mu1" , 
                                                                                        label = 'Choose a Number for \\mu_1' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 2 ,
                                                                                        value = 0.6) ,
                                                                            sliderInput(inputId = "AFib_mu2" , 
                                                                                        label = 'Choose a Number for \\mu_2 - \\mu_1' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 3 ,
                                                                                        value = 0.5),
                                                                            sliderInput(inputId = "AFib_mu3" , 
                                                                                        label = 'Choose a Number for \\mu_3 - \\mu_2' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 1 ,
                                                                                        value = 0.5),
                                                                            sliderInput(inputId = "AFib_sigma1" , 
                                                                                        label = 'Choose a Number for \\sigma_1' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 0.3 ,
                                                                                        value = 0.2),
                                                                            sliderInput(inputId = "AFib_sigma2" , 
                                                                                        label = 'Choose a Number for \\sigma_2' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 0.3 ,
                                                                                        value = 0.2 ),
                                                                            sliderInput(inputId = "AFib_sigma3" , 
                                                                                        label = 'Choose a Number for \\sigma_3' , 
                                                                                        min = 0 ,
                                                                                        step = 0.001,
                                                                                        max = 0.3 ,
                                                                                        value = 0.2 )
                                                       ) ),
                                                       column(4 , plotOutput( outputId = 'AFib_RRDenisty' , inline = T)),
                                                       column(4 , plotOutput( outputId = 'AFib_RRTimes', inline = T))
                                                     ),
                                                     tags$h2('Simulated ECG'),
                                                     tags$hr(),
                                                     fluidRow(actionButton(inputId = 'AFib_UpdateECG' , label = 'Update')),
                                                     fluidRow(plotOutput( outputId = 'AFib_ECG', inline = F))
                                          )                                          )
                  )
)
  
  
  server <- function( input , output ){
    # Regular elements
    {
      X_Reg <- reactive({CreateDefaultX_PACS( max(0,input$Reg_pi1 ) ,
                                              input$Reg_pi2,
                                              max(0,(1-input$Reg_pi1 - input$Reg_pi2) ) ,
                                              input$Reg_mu1 ,
                                              input$Reg_mu1 + input$Reg_mu2 , 
                                              input$Reg_mu1 + input$Reg_mu2 + input$Reg_mu3 , 
                                              input$Reg_sigma1,
                                              input$Reg_sigma2,
                                              input$Reg_sigma3) })
      # PACS RRDensity plot
      output$Reg_RRDenisty <- renderPlot({
        x <- seq(X_Reg()[4] - (3*X_Reg()[7]),X_Reg()[6] + (3*X_Reg()[9]),0.001)  
        f_x <- FM_EvaluateDenistyEstimate(x , X_Reg())
        
        p1_Reg <- ggplot(data.frame(RR = x , f_x = f_x) , aes(RR , f_x)) + geom_line(col = 'blue') +ggtitle('Distribution of RRTimes')
        p1_Reg
      },width = 450,height = 450)
      output$Reg_RRTimes <- renderPlot({
        
        RRTimes <- FM_SampleGMM(X_Reg() , 1000)
        
        p2_Reg <- ggplot(data.frame(t = cumsum(RRTimes) , RRTimes = RRTimes) , aes(t , RRTimes)) + geom_point(col =rgb(0,0,1,0.1)) + ylim(c(0,2)) +ggtitle('RRTimes')
        p2_Reg 
        
      },width = 450,height = 450)
      observeEvent( c(input$Reg_UpdateECG , X_Reg()) , {
        output$Reg_ECG <- renderPlot({
          
          
          RRTimes <- FM_SampleGMM(X_Reg() , 20)
          t = cumsum(RRTimes[,1])
          ECG <- PER_CreateECGReg( t , t_observation , RRTimes )
          p3_Reg <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
          p3_Reg <- p3_Reg + theme(
            panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                            size = 2, linetype = "solid"),
            panel.grid.major = element_line(size = 1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25)), 
            panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25))
          ) 
          
          p3_Reg <-p3_Reg + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
          p3_Reg
        }
        ,width = 1800,height = 200)
        
      }) 
    }
    # Bigeminy elements
    {X <- reactive({CreateDefaultX(input$mu1 , input$mu1 + input$mu2 , input$sigma1 , max(0.0001,input$sigma1 + input$sigma2 )) })
      
      output$RRDenisty <- renderPlot({
        x <- seq(X()[4] - (3*X()[7]),X()[6] + (3*X()[9]),0.001)  
        f_x <- FM_EvaluateDenistyEstimate(x , X())
        
        p1 <- ggplot(data.frame(RR = x , f_x = f_x) , aes(RR , f_x)) + geom_line(col = 'blue') +ggtitle('Distribution of RRTimes')
        p1
      },width = 450,height = 450)
      output$RRTimes <- renderPlot({
        
        RRTimes <- FM_SampleGMMBigeminy(X() , 1000)
        
        p2 <- ggplot(data.frame(t = cumsum(RRTimes) , RRTimes = RRTimes) , aes(t , RRTimes)) + geom_point(col =rgb(0,0,1,0.1)) + ylim(c(0,2)) +ggtitle('RRTimes')
        p2 
        
      }
      ,width = 450,height = 450)
      
      observeEvent( c(input$UpdateECG , X()) , {
        output$ECG <- renderPlot({
          
          
          RRTimes <- FM_SampleGMMBigeminy(X() , 20)
          t = cumsum(RRTimes)
          ECG <- PER_CreateECG( t , t_observation , RRTimes )
          p3 <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
          p3 <- p3 + theme(
            panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                            size = 2, linetype = "solid"),
            panel.grid.major = element_line(size = 1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25)), 
            panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25))
          ) 
          
          p3 <-p3 + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
          p3
        }
        ,width = 1800,height = 200)
        
      })
      
    }
    # PACs elements
    {
      X_Pacs <- reactive({CreateDefaultX_PACS( max(0,(1-input$PAC_pi1)/2 ) ,
                                               input$PAC_pi1,
                                               max(0,(1-input$PAC_pi1)/2 ) ,
                                               input$PAC_mu1 ,
                                               input$PAC_mu1 + input$PAC_mu2 , 
                                               input$PAC_mu1 + input$PAC_mu2 + input$PAC_mu3 , 
                                               input$PAC_sigma1,
                                               input$PAC_sigma2,
                                               input$PAC_sigma3) })
      # PACS RRDensity plot
      output$PACs_RRDenisty <- renderPlot({
        x <- seq(X_Pacs()[4] - (3*X_Pacs()[7]),X_Pacs()[6] + (3*X_Pacs()[9]),0.001)  
        f_x <- FM_EvaluateDenistyEstimate(x , X_Pacs())
        
        p1_Pacs <- ggplot(data.frame(RR = x , f_x = f_x) , aes(RR , f_x)) + geom_line(col = 'blue') +ggtitle('Distribution of RRTimes')
        p1_Pacs
      },width = 450,height = 450)
      output$PACs_RRTimes <- renderPlot({
        
        RRTimes <- PER_SampleGMMPACs(X_Pacs() , 1000)[1:1000,1]
        
        p2_PACs <- ggplot(data.frame(t = cumsum(RRTimes) , RRTimes = RRTimes) , aes(t , RRTimes)) + geom_point(col =rgb(0,0,1,0.1)) + ylim(c(0,2)) +ggtitle('RRTimes')
        p2_PACs 
        
      },width = 450,height = 450)
      observeEvent( c(input$PACs_UpdateECG , X_Pacs()) , {
        output$PACs_ECG <- renderPlot({
          
          
          RRTimes <- PER_SampleGMMPACs(X_Pacs() , 20)[1:20,]
          t = cumsum(RRTimes[,1])
          ECG <- PER_CreateECGPACS( t , t_observation , RRTimes )
          p3_Pacs <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
          p3_Pacs <- p3_Pacs + theme(
            panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                            size = 2, linetype = "solid"),
            panel.grid.major = element_line(size = 1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25)), 
            panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25))
          ) 
          
          p3_Pacs <-p3_Pacs + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
          p3_Pacs
        }
        ,width = 1800,height = 200)
        
      }) 
    }
    #Atrial Frillation Elements
    {
      X_AFib <- reactive({CreateDefaultX_PACS( max(0,input$AFib_pi1 ) ,
                                               input$AFib_pi2,
                                               max(0,(1-input$AFib_pi1 - input$AFib_pi2) ) ,
                                               input$AFib_mu1 ,
                                               input$AFib_mu1 + input$AFib_mu2 , 
                                               input$AFib_mu1 + input$AFib_mu2 + input$AFib_mu3 , 
                                               input$AFib_sigma1,
                                               input$AFib_sigma2,
                                               input$AFib_sigma3) })
      # PACS RRDensity plot
      output$AFib_RRDenisty <- renderPlot({
        x <- seq(X_AFib()[4] - (3*X_AFib()[7]),X_AFib()[6] + (3*X_AFib()[9]),0.001)  
        f_x <- FM_EvaluateDenistyEstimate(x , X_AFib())
        
        p1_AFib <- ggplot(data.frame(RR = x , f_x = f_x) , aes(RR , f_x)) + geom_line(col = 'blue') +ggtitle('Distribution of RRTimes')
        p1_AFib
      },width = 450,height = 450)
      output$AFib_RRTimes <- renderPlot({
        
        RRTimes <- FM_SampleGMM(X_AFib() , 1000)
        RRTimes[RRTimes < 0.3] = 0.3
        
        p2_PACs <- ggplot(data.frame(t = cumsum(RRTimes) , RRTimes = RRTimes) , aes(t , RRTimes)) + geom_point(col =rgb(0,0,1,0.1)) + ylim(c(0,2)) +ggtitle('RRTimes')
        p2_PACs 
        
      },width = 450,height = 450)
      observeEvent( c(input$AFib_UpdateECG , X_AFib()) , {
        output$AFib_ECG <- renderPlot({
          
          
          RRTimes <- FM_SampleGMM(X_AFib() , 20)
          RRTimes[RRTimes < 0.3] = 0.3
          t = cumsum(RRTimes[,1])
          ECG <- PER_CreateECGAFib( t , t_observation , RRTimes )
          p3_AFib <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
          p3_AFib <- p3_AFib + theme(
            panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                            size = 2, linetype = "solid"),
            panel.grid.major = element_line(size = 1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25)), 
            panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25))
          ) 
          
          p3_AFib <-p3_AFib + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
          p3_AFib
        }
        ,width = 1800,height = 200)
        
      }) 
    }
  }
}
shinyApp(ui = ui , server = server)
