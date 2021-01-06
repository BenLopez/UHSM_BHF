
# Prelimiaries 
{
  pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
  source("LibrariesAndSettings.R" , print.eval  = TRUE )
  DP_LoadPatientIndex()
  DP_ChooseDataReps()
  FilestoProcess <- c( "ECGI" ,  "ECGII" , "ECGIII")#DP_ChooseECGstoProcess() 
  HoursBeforeandAfter <- setNames(list(5 , 1) , c('numberhoursbefore' , 'numberhoursafter')) #DP_SelectHoursBeforeandAfter() 
  if(sum(listAllPatients == 'z1007') > 0){
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)}
  set.seed( 1 )
}

{
  # User Interface
  ui <-fluidPage(title = 'ECG Examination',
                 tags$h1('UHSM ECG Examination Application'),
                 selectInput(inputId = 'PatinetId' ,
                             label = 'Select patient to view' , 
                             choices = listAllPatients , selected = listAllPatients[[1]]  ),
                 uiOutput(outputId = 'Slider'),
                 #sliderInput(inputId = "SelectTime" , 
                 #           label = 'Select a Time' , 
                 #           min = 0 ,
                 #           step = 0.03,
                 #           max = 100 ,
                 #            value = 50 ,
                 #            width = 1500),
                 fluidRow(title = 'Title' , textOutput(outputId = 'title' )),
                 fluidRow(title = 'ECGs' , plotOutput(outputId = 'ECGI'  )) ,
                 #fluidRow(title = 'break' , tags$h2(' ')),
                 #fluidRow(title = 'break' , tags$h2(' ')),
                 #plotOutput(outputId = 'ECGII'  ) ,
                 #plotOutput(outputId = 'ECGIII'  ) ,
                 fluidRow(title = 'RRTimes',plotOutput(outputId = 'RRTimesPlot'  ))
  )
  
  # Actions
  server <- function( input , output ){
    
    # Load ECG and RPeak data 
    {
      ECGs <- reactive(DP_LoadReducedECGs( path , input$PatinetId , numberrep = numberrep , FilestoProcess = FilestoProcess))
      RPeaksStruct <- reactive(DP_LoadRpeaksfile(path , input$PatinetId ))
      MetaData <- reactive(DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode =  input$PatinetId))
      rangeoftimes   <- reactive(max(RPeaksStruct()$RRCombined$t) - min(RPeaksStruct()$RRCombined$t))
      #timetoview <- reactive( min(RPeaksStruct()$RRCombined$t) + (input$SelectTime/100)*rangeoftimes() )
    }
    
    output$Slider <- renderUI({
      sliderInput(inputId = "SelectTime" , 
                  label = 'Select a Time' , 
                  min = min(RPeaksStruct()$RRCombined$t) ,
                  step = 10,
                  max = max(RPeaksStruct()$RRCombined$t) ,
                  value = RPeaksStruct()$RRCombined$t[1000] ,
                  width = 1440)
    })
    # Time to view
    {
      timetoview <- reactive( input$SelectTime )
    }
    {
      output$title <- renderText(input$PatinetId) 
    }
    # RR Times plot
    output$RRTimesPlot <- renderPlot({
      
      if(nrow(MetaData()) ==0){
        RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeaksStruct() ,
                                            MetaData = DP_CreateDummyMetaData(PatIndex2017 , input$PatinetId))
        RRPlot <- BC_PlotAddViewingRegionLines(RRPlot , c(timetoview() , timetoview() + 10))  
      }else{
        #RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeaksStruct() , MetaData = DP_CreateDummyMetaData(PatIndex2017 , input$PatinetId))
        RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeaksStruct() , MetaData = MetaData() )
        RRPlot <- BC_PlotAddViewingRegionLines(RRPlot , c(timetoview() , timetoview() + 10))
      }
      print(RRPlot)
    },width = 1500,height = 178 )
    # ECGI Times plot
    output$ECGI <- renderPlot({ 
      ECGIPlot <- BC_PlotCreateECGPlots(RPeaksStruct() ,  ECGs()$ECGI  , timestart = timetoview() , ECGindex = 1, timeindex = 10)
      ECGIIPlot <- BC_PlotCreateECGPlots(RPeaksStruct() ,  ECGs()$ECGII  , timestart = timetoview() , ECGindex = 2, timeindex = 10)
      ECGIIIPlot <- BC_PlotCreateECGPlots(RPeaksStruct() ,  ECGs()$ECGIII  , timestart = timetoview() , ECGindex = 3, timeindex = 10)
      print(grid.arrange(ECGIPlot , ECGIIPlot , ECGIIIPlot , nrow = 3 ,ncol = 1  ))
    },width = 1500,height = 300)
    output$ECGII <- renderPlot({ 
      ECGIIPlot <- BC_PlotCreateECGPlots(RPeaksStruct() ,  ECGs()$ECGII  , timestart = timetoview() , ECGindex = 2, timeindex = 10)
      #ECGIIPlot <- BC_PlotECGAddTitle( ECGIIPlot, timetoview() )
      print(ECGIIPlot)
    },width = 1500,height = 178)
    output$ECGIII <- renderPlot({ 
      ECGIIIPlot <- BC_PlotCreateECGPlots(RPeaksStruct() ,  ECGs()$ECGIII  , timestart = timetoview() , ECGindex = 3, timeindex = 10)
      #ECGIIIPlot <- BC_PlotECGAddTitle( ECGIIIPlot  , timetoview() )
      print(ECGIIIPlot)
    },width = 1500,height = 178)
  }
  
  # ECGII Times plot
  
  # ECGIII Times plot
}


shinyApp(ui = ui , server = server)
