
# Prelimiaries 
{
  source("LibrariesAndSettings.R" , print.eval  = TRUE )
  t_observation = seq(0.25  , 10 , 0.005)
  p_t <- seq(0.5 , 0.95,0.0005)
}

{
  ui <-fluidPage(title = 'Heart-Rhythm Models Elicitation Tool',
                 withMathJax(),
                 tabsetPanel( navbarMenu(title = 'Rhythms',
                                         {tabPanel(title = 'Home' ,
                                                   withMathJax(), 
                                                   wellPanel(tags$h1('Model Description')),
                                                   wellPanel(tags$h2('Bigeminy')),
                                                   tags$p('Bigeminy is a heart rhythm in which there are repeated beats, one long and one short. The model for bigeminy is a mixtures model of two Gaussian distribtutions. The probability density function (pdf) is given by'),
                                                   withMathJax(tags$p('$$f(RR , x) = 0.5\\Phi(RR , \\mu_{(1)} , \\sigma_{(1)}) + 0.5\\Phi(RR , \\mu_{(2)} , \\sigma_{(2)}).$$')) ,
                                                   withMathJax(tags$p('Here $$\\Phi(RR , \\mu , \\sigma)$$ is a Gaussian probability density function, evaluated at RR, with mean and standard deviation given by $$(\\mu , \\sigma).$$  ')) ,
                                                   withMathJax(tags$p('The four parameters  $$(\\mu_{(1)} , \\sigma_{(1)} , \\mu_{(2)} , \\sigma_{(2)})$$ control the average and the range of RR-times for the short and long beats. The first two control the average and range for the short beats and the second two the average and range of the second beat.')),
                                                   imageOutput(outputId = 'BigeminyDensity'),
                                                   wellPanel(tags$h2('PACs')),
                                                   tags$p('Premature-atrial-complexes (PACs), are heartbeats that arise within the atria of the heart which are outside the NSR. PACs are early (that is, premature) electrical impulses that are generated within the cardiac atria, but not from the sinus node. PACs momentarily interrupt the normal sinus rhythm with a fast beat followed by a compensatory pause. The model for PACs is a mixture of three Gaussian distributions. One compnent is the fast beats, one the normal and one the slow.'),
                                                   withMathJax(tags$p('$$f(RR , x) = \\frac{(1-\\pi_{(1)})}{2}\\Phi(RR , \\mu_{(1)} , \\sigma_{(1)}) + \\pi_{(1)}\\Phi(RR , \\mu_{(2)} , \\sigma_{(2)})+ \\frac{(1-\\pi_{(1)})}{2}\\Phi(RR , \\mu_{(3)} , \\sigma_{(3)}).$$')),
                                                   tags$p('The first component defines the distribution of the fast beats. The second component the distrbution of the normal beats and the third the distribution of the slow beats. The seven parameters to control the RR-times dstribution for PAcs are then'),
                                                   withMathJax(tags$p('$$(\\pi_{(1)} ,  \\mu_{(1)} ,  \\mu_{(2)} ,  \\mu_{(3)} , \\sigma_{(1)} , \\sigma_{(2)} , \\sigma_{(3)}).$$')),
                                                   tags$p('The first parameter controls the proportion of PACs to normal beats. The second to fourth control the average speed of the fast, normal and slow beats. The fith to seventh control the width of the fast, normal and slow beats.'),
                                                   imageOutput(outputId = 'PACsDensity'),
                                                   wellPanel(tags$h2('Regular')),
                                                   tags$p('Sinus rhythm is any rhythm where the depolarisation begins in the sinus node. Normal sinus rhythm is where all other components of the ECG are also not unusual. The density of the RR-times distribution is given by'),
                                                   withMathJax(tags$p('$$f(RR , x) = P(RR, \\mu ,\\sigma^2, s , k ).$$')),
                                                   tags$p('Here P is a pearson distribution with expectation, variance, skewness and kurtosis given by'),
                                                   withMathJax(tags$p('$$(\\mu ,\\sigma^2, s , k).$$')),
                                                   wellPanel(tags$h2('Atrial Fibrillation (AFib)')),
                                                   tags$p('Atrial Fibrillation (AFib) is an abnormal heart rhythm characterized by irregular beating (fluttering or fibrillation) of the atria and irregularly-irregular ventricular contractions.  The density of the RR-times distribtion is given by'),
                                                   withMathJax(tags$p('$$f(RR , x) = (\\pi_{(1)})\\Phi(RR , \\mu_{(1)} , \\sigma_{(1)}) + \\pi_{(2)}\\Phi(RR , \\mu_{(2)} , \\sigma_{(2)})+ (1 - \\pi_{(1)} - \\pi_{(2)} )\\Phi(RR , \\mu_{(3)} , \\sigma_{(3)}).$$')),
                                                   withMathJax(tags$p('$$(\\pi_{(1)} ,  \\pi_{(2)}, \\mu_{(1)} ,  \\mu_{(2)} ,  \\mu_{(3)} , \\sigma_{(1)} , \\sigma_{(2)} , \\sigma_{(3)}).$$')),
                                                   tags$p('The model is restricted so that $$\\mu_{(1)} <  \\mu_{(2)} <  \\mu_{(3)}$$ to make the specification easier. It is the same form as the Regular model to give flexibility in the shape of the RR-times distribution. The first two parameters controls the proportion of beats from component one and two. The second to fourth control the average speed of the fast, normal and slow beats. The fith to seventh control the width of the fast, normal and slow beats.'),
                                                   imageOutput(outputId = 'AFibDensity')
                                         )},# Front page
                                         # {tabPanel(title = 'Regular',
                                         #           wellPanel(tags$h1('Regular')) ,withMathJax(),
                                         #           fluidRow(
                                         #             column(4 , wellPanel(sliderInput(inputId = "Reg_pi1" , 
                                         #                                              label = '$$\\pi_{(1)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 1 ,
                                         #                                              value = 0 ) ,
                                         #                                  sliderInput(inputId = "Reg_pi2" , 
                                         #                                              label = '$$\\pi_{(2)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 1 ,
                                         #                                              value = 1) ,
                                         #                                  #sliderInput(inputId = "PAC_pi3" , 
                                         #                                  #            label = '$$ \\pi_{(3)} $$' , 
                                         #                                  #           min = 1 ,
                                         #                                  #            step = 0.001,
                                         #                                  #            max = 2 ,
                                         #                                  #            value = 0.1) ,
                                         #                                  sliderInput(inputId = "Reg_mu1" , 
                                         #                                              label = '$$\\mu_{(1)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 3 ,
                                         #                                              value = 0.75) ,
                                         #                                  sliderInput(inputId = "Reg_mu2" , 
                                         #                                              label = '$$\\mu_{(2)} - \\mu_{(1)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 0.25 ,
                                         #                                              value = 0.05),
                                         #                                  sliderInput(inputId = "Reg_mu3" , 
                                         #                                              label = '$$\\mu_{(3)} - \\mu_{(2)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 0.25 ,
                                         #                                              value = 0.05),
                                         #                                  sliderInput(inputId = "Reg_sigma1" , 
                                         #                                              label = '$$\\sigma_{(1)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 0.1 ,
                                         #                                              value = 0.01),
                                         #                                  sliderInput(inputId = "Reg_sigma2" , 
                                         #                                              label = '$$\\sigma_{(2)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 0.1 ,
                                         #                                              value = 0.01 ),
                                         #                                  sliderInput(inputId = "Reg_sigma3" , 
                                         #                                              label = '$$\\sigma_{(3)}$$' , 
                                         #                                              min = 0 ,
                                         #                                              step = 0.001,
                                         #                                              max = 0.1 ,
                                         #                                              value = 0.01 )
                                         #             ) ),
                                         #             column(4 , plotOutput( outputId = 'Reg_RRDenisty' , inline = T)),
                                         #             column(4 , plotOutput( outputId = 'Reg_RRTimes', inline = T))
                                         #           ),
                                         #           tags$h2('Simulated ECG'),
                                         #           tags$hr(),
                                         #           fluidRow(actionButton(inputId = 'Reg_UpdateECG' , label = 'Update')),
                                         #           fluidRow(plotOutput( outputId = 'Reg_ECG', inline = F))
                                         # )},# Regular
                                         {tabPanel(title = 'Regular',
                                                   wellPanel(tags$h1('Regular')) ,withMathJax(),
                                                   fluidRow(
                                                     column(4 , wellPanel(
                                                                          #sliderInput(inputId = "Reg_pi1" , 
                                                                          #            label = '$$\\pi_{(1)}$$' , 
                                                                          #            min = 0 ,
                                                                          #            step = 0.001,
                                                                          #            max = 1 ,
                                                                          #            value = 0 ) ,
                                                                          #sliderInput(inputId = "Reg_pi2" , 
                                                                          #            label = '$$\\pi_{(2)}$$' , 
                                                                          #            min = 0 ,
                                                                          #            step = 0.001,
                                                                          #            max = 1 ,
                                                                          #            value = 1) ,
                                                                          #sliderInput(inputId = "PAC_pi3" , 
                                                                          #            label = '$$ \\pi_{(3)} $$' , 
                                                                          #           min = 1 ,
                                                                          #            step = 0.001,
                                                                          #            max = 2 ,
                                                                          #            value = 0.1) ,
                                                                          sliderInput(inputId = "Reg_Mean" , 
                                                                                      label = '$$\\mu$$' , 
                                                                                      min = 0.4 ,
                                                                                      step = 0.001,
                                                                                      max = 2 ,
                                                                                      value = 1) ,
                                                                          sliderInput(inputId = "Reg_StandardDeviation" , 
                                                                                      label = '$$\\sigma$$' , 
                                                                                      min = 0.00000001 ,
                                                                                      step = 0.001,
                                                                                      max = 0.5 ,
                                                                                      value = 0.05),
                                                                          uiOutput(outputId = 'Skewness_Slider'),
                                                                          #sliderInput(inputId = "Reg_Skewness" , 
                                                                          #            label = '$$Skewness$$' , 
                                                                          #            min = -3 ,
                                                                          #            step = 0.01,
                                                                          #            max = 3 ,
                                                                          #            value = 0),
                                                                          sliderInput(inputId = "Reg_Kurtosis" , 
                                                                                      label = '$$Kurtosis$$' , 
                                                                                      min = 1 ,
                                                                                      step = 0.05,
                                                                                      max = 100 ,
                                                                                      value = 1.8)
                                                     ) ),
                                                   column(4 , plotOutput( outputId = 'Reg_RRDenisty' , inline = T)),
                                                   column(4 , plotOutput( outputId = 'Reg_RRTimes', inline = T))
                                         ),
                                           tags$h2('Simulated ECG'),
                                           tags$hr(),
                                           fluidRow(actionButton(inputId = 'Reg_UpdateECG' , label = 'Update')),
                                           fluidRow(plotOutput( outputId = 'Reg_ECG', inline = F))
                 )},# Regular
                 {tabPanel(title = 'Bigeminy', 
                           wellPanel(tags$h1('Bigeminy')) ,
                           fluidRow(
                             column(4 , wellPanel(sliderInput(inputId = "mu1" , 
                                                              label = '$$ \\mu_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 2 ,
                                                              value = 0.5) ,
                                                  sliderInput(inputId = "mu2" , 
                                                              label = '$$ \\mu_{(2)} - \\mu_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 2 ,
                                                              value = 0.2),
                                                  sliderInput(inputId = "sigma1" , 
                                                              label = '$$ \\sigma_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 0.1 ,
                                                              value = 0.01),
                                                  sliderInput(inputId = "sigma2" , 
                                                              label = '$$ \\sigma_{(2)} - \\sigma_{(1)}$$' , 
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
                 )},# Bigeminy
                 {tabPanel(title = 'PACs', 
                           wellPanel(tags$h1('PACs')) ,
                           fluidRow(
                             column(4 , wellPanel(sliderInput(inputId = "PAC_pi1" , 
                                                              label = '$$ \\pi_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 1 ,
                                                              value = 0.8) ,
                                                  #sliderInput(inputId = "PAC_pi2" , 
                                                  #            label = '$$ \\pi_{(2)} $$' , 
                                                  #            min = 0 ,
                                                  #            step = 0.001,
                                                  #            max = 1 ,
                                                  #            value = 0.8) ,
                                                  #sliderInput(inputId = "PAC_pi3" , 
                                                  #            label = '$$ \\pi_{(3)} $$' , 
                                                  #           min = 1 ,
                                                  #            step = 0.001,
                                                  #            max = 2 ,
                                                  #            value = 0.1) ,
                                                  sliderInput(inputId = "PAC_mu1" , 
                                                              label = '$$ \\mu_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 2 ,
                                                              value = 0.6) ,
                                                  sliderInput(inputId = "PAC_mu2" , 
                                                              label = '$$ \\mu_{(2)}  - \\mu_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 2 ,
                                                              value = 0.1),
                                                  sliderInput(inputId = "PAC_mu3" , 
                                                              label = '$$ \\mu_{(3)}  - \\mu_{(2)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 2 ,
                                                              value = 0.1),
                                                  sliderInput(inputId = "PAC_sigma1" , 
                                                              label = '$$ \\sigma_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 0.1 ,
                                                              value = 0.01),
                                                  sliderInput(inputId = "PAC_sigma2" , 
                                                              label = '$$ \\sigma_{(2)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 0.1 ,
                                                              value = 0.01 ),
                                                  sliderInput(inputId = "PAC_sigma3" , 
                                                              label = '$$ \\sigma_{(3)} $$' , 
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
                 )},# PACs
                 {tabPanel(title = 'Atrial Fibrillation', 
                           wellPanel(tags$h1('Atrial Fibrillation')) ,
                           fluidRow(
                             column(4 , wellPanel(sliderInput(inputId = "AFib_pi1" , 
                                                              label = '$$ \\pi_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 1 ,
                                                              value = 0.33333) ,
                                                  sliderInput(inputId = "AFib_pi2" , 
                                                              label = '$$ \\pi_{(2)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 1 ,
                                                              value = 0.33333) ,
                                                  #sliderInput(inputId = "PAC_pi3" , 
                                                  #            label = '$$ \\pi_{(3)} $$' , 
                                                  #           min = 1 ,
                                                  #            step = 0.001,
                                                  #            max = 2 ,
                                                  #            value = 0.1) ,
                                                  sliderInput(inputId = "AFib_mu1" , 
                                                              label = '$$ \\mu_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 2 ,
                                                              value = 0.6) ,
                                                  sliderInput(inputId = "AFib_mu2" , 
                                                              label = '$$ \\mu_{(2)}  - \\mu_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 3 ,
                                                              value = 0.5),
                                                  sliderInput(inputId = "AFib_mu3" , 
                                                              label = '$$ \\mu_{(3)}  - \\mu_{(2)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 1 ,
                                                              value = 0.5),
                                                  sliderInput(inputId = "AFib_sigma1" , 
                                                              label = '$$ \\sigma_{(1)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 0.3 ,
                                                              value = 0.2),
                                                  sliderInput(inputId = "AFib_sigma2" , 
                                                              label = '$$ \\sigma_{(2)} $$' , 
                                                              min = 0 ,
                                                              step = 0.001,
                                                              max = 0.3 ,
                                                              value = 0.2 ),
                                                  sliderInput(inputId = "AFib_sigma3" , 
                                                              label = '$$ \\sigma_{(3)} $$' , 
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
                 )},# AFib
                 {tabPanel(
                   title = 'P-waves' , 
                   wellPanel(tags$h1('P-wave Model')),
                   fluidRow(
                     column(4 , wellPanel(sliderInput(inputId = "Pwave_AL" , 
                                                      label = 'Choose left P-amplitude' , 
                                                      min = -40 ,
                                                      step = 0.001,
                                                      max = 40 ,
                                                      value = 6) ,
                                          sliderInput(inputId = "Pwave_AR" , 
                                                      label = 'Choose right P-amplitude' , 
                                                      min = -40 ,
                                                      step = 0.001,
                                                      max = 40 ,
                                                      value = 10.21) ,
                                          sliderInput(inputId = "Pwave_LL" , 
                                                      label = 'Choose left P-time before R' , 
                                                      min =  0.00001 ,
                                                      step = 0.001,
                                                      max = 0.5 ,
                                                      value = 0.225) ,
                                          sliderInput(inputId = "Pwave_LR" , 
                                                      label = 'Choose right P-time before R' , 
                                                      min = 0.00001 ,
                                                      step = 0.001,
                                                      max = 0.5 ,
                                                      value = 0.175),
                                          sliderInput(inputId = "Pwave_WL" , 
                                                      label = 'Choose left P-width' , 
                                                      min = 0 ,
                                                      step = 0.00001,
                                                      max = 0.1 ,
                                                      value = 0.03),
                                          sliderInput(inputId = "Pwave_WR" , 
                                                      label = 'Choose right P-width' , 
                                                      min = 0 ,
                                                      step = 0.00001,
                                                      max = 0.1 ,
                                                      value = 0.03)
                     ) ),
                     column(4 , plotOutput( outputId = 'Pwave_RRDenisty' ))
                   ) 
                 )  }, # P-waves
                 {tabPanel(
                   title = 'Heart-rhythm Discrepancy' , 
                   wellPanel(tags$h1('Heart-rhythm Discrepancy')),
                   fluidRow(
                     column(4 , wellPanel(sliderInput( inputId = "HeartrhythmDisc" , 
                                                       label = 'Choose dispersion parameter.' , 
                                                       min = 0 ,
                                                       step = 1,
                                                       max = 10000 ,
                                                       value = 5000),
                                          actionButton(inputId = 'HRDis_UpdateRRTimes' , label = 'Update RRtimes')
                     ) ),
                     
                     column(4 , plotOutput( outputId = 'HRDis_RRDenisty' , inline = T)),
                     column(4 , plotOutput( outputId = 'HRDis_RRTimes', inline = T))
                   ),
                   tags$h2('Simulated ECG'),
                   tags$hr(),
                   fluidRow(actionButton(inputId = 'HRDis_UpdateECG' , label = 'Update')),
                   fluidRow(plotOutput( outputId = 'HRDis_ECG', inline = F)))
                 })
  )
  )
  
  
  server <- function( input , output ){
    # Description page
    {{
      output$BigeminyDensity <- renderImage({
        
        list(src = 'www/BigeminyRRTimesDistribution.png',
             contentType = 'image/png',
             width = 400,
             height = 400,
             alt = "This is alternate text")
      }, deleteFile = FALSE)
    }
      {
        output$PACsDensity <- renderImage({
          
          list(src = 'www/PACsRRTimesDistribution.png',
               contentType = 'image/png',
               width = 400,
               height = 400,
               alt = "This is alternate text")
        }, deleteFile = FALSE)
      }
      {
        output$RegularDensity <- renderImage({
          
          list(src = 'www/RegularRRTimesDistribution.png',
               contentType = 'image/png',
               width = 400,
               height = 400,
               alt = "This is alternate text")
        }, deleteFile = FALSE)
      }
      {
        output$AFibDensity <- renderImage({
          
          list(src = 'www/AFibRRTimesDistribution.png',
               contentType = 'image/png',
               width = 400,
               height = 400,
               alt = "This is alternate text")
        }, deleteFile = FALSE)
      }}
    # Regular elements
    {
      output$Skewness_Slider <- renderUI({
        sliderInput(inputId = "Reg_Skewness" , 
                    label = withMathJax('$$Skewness$$') , 
                    min = -round(sqrt((input$Reg_Kurtosis -1) - 0.01),2) ,
                    step = 0.01,
                    max = +round(sqrt((input$Reg_Kurtosis -1) + 0.01),2) ,
                    value = 0)
      })
      X_Reg <- reactive({c(input$Reg_Mean,
                           input$Reg_StandardDeviation,
                           input$Reg_Skewness,
                           input$Reg_Kurtosis)
      })
      # PACS RRDensity plot
      output$Reg_RRDenisty <- renderPlot({
        x <- seq(X_Reg()[1] - (3*X_Reg()[2]),X_Reg()[1] + (3*X_Reg()[2]),0.0001)  
        f_x <- FM_PearonsRegularDenisty( X_Reg(),x)
        
        p1_Reg <- ggplot(data.frame(RR = x , f_x= log(f_x)) , aes(RR , f_x)) +
          geom_line(col = 'blue') +ggtitle('Distribution of RRTimes') +
          xlab('rr')+
          ylab('log(f(rr))')
        p1_Reg
      },width = 450,height = 450)
       output$Reg_RRTimes <- renderPlot({
         
         RRTimes <- FM_CleanRRTimes(FM_SamplePearonsRegular(X_Reg() , 1000))
         
         p2_Reg <- ggplot(data.frame(t = cumsum(RRTimes) , RRTimes = RRTimes) , aes(t , RRTimes)) + geom_point(col =rgb(0,0,1,0.1)) + ylim(c(0,2)) +ggtitle('RRTimes')
         p2_Reg 
         
       },width = 450,height = 450)
       observeEvent( c(input$Reg_UpdateECG , X_Reg()) , {
         output$Reg_ECG <- renderPlot({
           
           
           RRTimes <- FM_CleanRRTimes( FM_SamplePearonsRegular(X_Reg() , 100) )
           t = cumsum(RRTimes)
           ECG <- PER_CreateECGReg( t , t_observation , RRTimes )
           p3_Reg <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) +
             geom_line(col =rgb(0,0,0,0.9) , size = 0.7) +
             theme(
             panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                             size = 2, linetype = "solid"),
             panel.grid.major = element_line(size = 1, linetype = 'solid',
                                             colour = rgb(1,0,0,0.25)), 
             panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                             colour = rgb(1,0,0,0.25))
           ) 
           
           p3_Reg <-p3_Reg + 
             scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + 
             scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
           p3_Reg
         }
         ,width = 1500,height = 178)
         
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
        ,width = 1500,height = 178)
        
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
        ,width = 1500,height = 178)
        
      }) 
    }
    #Atrial Fibrillation Elements
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
        ,width = 1500,height = 178)
        
      }) 
    }
    # Heart rhythm discreapncy elements
    {
      { 
        #X_PwaveDiscrepancy <- reactive({CreateDefaultX_PACS( max(0,input$AFib_pi1 ) ,
        #                                                     input$AFib_pi2,
        #                                                     max(0,(1-input$AFib_pi1 - input$AFib_pi2) ) ,
        #                                                     input$AFib_mu1 ,
        #                                                     input$AFib_mu1 + input$AFib_mu2 , 
        #                                                     input$AFib_mu1 + input$AFib_mu2 + input$AFib_mu3 , 
        #                                                     input$AFib_sigma1,
        #                                                     input$AFib_sigma2,
        #                                                     input$AFib_sigma3) })
        # PACS RRDensity plot
        X_PwaveDiscrepancy <-  reactive({CreateDefaultX_PACS( 0 ,
                                                              1,
                                                              0 ,
                                                              0.1 ,
                                                              0.8 , 
                                                              0.1 , 
                                                              0.1,
                                                              0.035,
                                                              0.1) })
        xdis <- reactive({seq(X_PwaveDiscrepancy()[4] - (3*X_PwaveDiscrepancy()[7]),X_PwaveDiscrepancy()[6] + (3*X_PwaveDiscrepancy()[9]),0.001)})  
        G_0 <- reactive({ function(N){ FM_SampleGMM( X = X_PwaveDiscrepancy() , N ) } })
        RRTimesDis <- eventReactive(c(input$HRDis_UpdateRRTimes , input$HeartrhythmDisc , X_PwaveDiscrepancy()) , {
          FM_SampleDP(input$HeartrhythmDisc , l = 0.005 , n = 10000 , function(N){ FM_SampleGMM( X = X_PwaveDiscrepancy() , N )})
        })
      } # reactive variables
      output$HRDis_RRDenisty <- renderPlot({
        
        kdeestimate = kde(RRTimesDis())
        DaseDesnity <- FM_EvaluateDenistyEstimate(kdeestimate$eval.points , X_PwaveDiscrepancy())
        p1_HRDis <- ggplot(data.frame(RR = kdeestimate$eval.points , f_x = kdeestimate$estimate) , aes(RR , f_x)) + geom_line(col = 'blue') +ggtitle('Distribution of RRTimes')
        p1_HRDis + geom_line(data = data.frame(RR = kdeestimate$eval.points , f_x = DaseDesnity) , aes(RR , f_x)  , col = 'black')
      },width = 450,height = 450)
      output$HRDis_RRTimes <- renderPlot({
        RRTimes2 <- sample(RRTimesDis() , 1000)
        RRTimes2[RRTimes2 < 0.3] = 0.3
        
        p2_HRDis <- ggplot(data.frame(t = cumsum(RRTimes2) , RRTimes = RRTimes2) , aes(t , RRTimes)) + geom_point(col =rgb(0,0,1,0.1)) + ylim(c(0,2)) +ggtitle('RRTimes')
        p2_HRDis 
        
      },width = 450,height = 450)
      observeEvent( c(input$HRDis_UpdateECG , input$HeartrhythmDisc) , {
        output$HRDis_ECG <- renderPlot({
          RRTimes <- sample(RRTimesDis() , 20)
          RRTimes[RRTimes < 0.3] = 0.3
          t = cumsum(RRTimes)
          ECG <- PER_CreateECGReg( t , t_observation , RRTimes )
          p3_HRDis <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
          p3_HRDis <- p3_HRDis + theme(
            panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                            size = 2, linetype = "solid"),
            panel.grid.major = element_line(size = 1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25)), 
            panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                            colour = rgb(1,0,0,0.25))
          ) 
          
          p3_HRDis <-p3_HRDis + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
          p3_HRDis
        }
        ,width = 1500,height = 178)
        
      })
    }
    # Pwave Elements
    {X_Pwave <- reactive({
      #X_Pwave <- PER_CreatePwaveECGXReg()
      PER_CreatePwaveECGXReg(t = 1, X1 = input$Pwave_AL , X2 = input$Pwave_AR , X5 = input$Pwave_WL ,  X6 = input$Pwave_WR , X3 = input$Pwave_LL , X4 = input$Pwave_LR)
    })
      
      p_t <- p_t
      output$Pwave_RRDenisty <- renderPlot({
        p1_Pwave <- ggplot(
          data.frame(t = p_t, V = ECGSim_WrapperSingleBeat_m2( X_Pwave() , p_t) )
          , aes(t,V) ) + 
          geom_line(col = 'blue') +
          ggtitle('P-wave Morphology')+
          ylim(c(-20,20))
        p1_Pwave <- p1_Pwave + geom_line(data = data.frame(t = p_t , V = ECGSim_RightPwave( X_Pwave() , p_t)) , aes(t , V) ) +
          geom_line(data = data.frame(t = p_t , V = ECGSim_LeftPwave( X_Pwave() , p_t)) , aes(t , V) )
        print(p1_Pwave)
      },width = 450,height = 450)
      
    }
  }
}
shinyApp(ui = ui , server = server)
