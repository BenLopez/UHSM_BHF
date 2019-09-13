check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Usage example
packages<-c('shinycssloaders',
            'shiny',
            'shinythemes',
            'ggplot2',
            'plotly',
            'mvtnorm',
            'expm',
            'metRology')
check.packages(packages)

library(shinycssloaders)
library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(mvtnorm)
library(expm)
library("metRology")

RenalFailure = read.csv("RenalFailure.csv")
PatientIndex = read.csv("PatientIndex.csv")
FlowSheet = read.csv("FlowSheet.csv")
Biochemistry = read.csv("Biochemistry.csv")

inputm = "0.55,-0.2"
inputC = "0.01,0.001"
inputG = "1,1,0,1"
inputFF = "1,0"
inputdeltamu = 0.80
inputdeltamunew = 0.12
inputdeltabeta = 0.90
inputdeltabetanew = 0.12
inputn = 20
inputS = 0.1
inputd = 2
inputdel = 0.95
inputdelnew = 0.80

ui <- navbarPage(
  theme = shinytheme("cerulean"),
  "DLM",
  navbarMenu("DLM",
             tabPanel("DLM",
                      sidebarLayout(
                        sidebarPanel(
                          fluidRow(
                            
                            column(5, 
                                   numericInput(inputId = "k", label = "Enter k", value = 6, min = 1)),
                            
                            column(7, 
                                   numericInput(inputId = "Observation", label = "Enter Observation", value = 0, min = 0))
                            
                          ),
                          
                          fluidRow(
                            
                            column(7, 
                                   selectInput(inputId = 'Id', label = 'Select Patient', unique(RenalFailure$NewPseudoId)))
                          ),
                          fluidRow(
                            textOutput("FirstAKI1UO"),
                            textOutput("FirstFilter"),
                            textOutput("Gender"),
                            textOutput("Age"),
                            textOutput("Weight"),
                            textOutput("Height"),
                            textOutput("EuroSCORE"),
                            textOutput("ProcDetails"),
                            textOutput("Urgency")
                          )
                          
                        ),
                        mainPanel(
                          verbatimTextOutput("JointProbability"),
                          plotOutput("scat"), 
                          tableOutput("df"),
                          plotlyOutput("scat3"),
                          plotlyOutput("scat4"),
                          plotlyOutput("scat6")
                        )
                      )
             )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  output$df <-renderTable({
    
    SubsetRenalFailureId = subset(RenalFailure, NewPseudoId == input$Id) 
    SubsetRenalFailureId2 = SubsetRenalFailureId[-1, ]
    PatientIndexId = subset(PatientIndex, NewPseudoId == input$Id)
    SubsetRenalFailureId2$Urine60 = SubsetRenalFailureId2$Urine60/PatientIndexId$Weight
    Y = SubsetRenalFailureId2$Urine60
    Yc = Y + 0.1
    y = log(Yc)
    
    Cdiag <- as.numeric(unlist(strsplit(inputC,",")))
    C <- array(diag(Cdiag, nrow = 2), dim = c(2,2,length(y) + 1))
    mvec <- as.numeric(unlist(strsplit(inputm,",")))
    m <- array(mvec, dim = c(2,length(y) + 1))
    d <- array(inputd, dim = length(y) + 1)
    n <- array(inputn, dim = length(y) + 1)
    Wmu <- array(0, dim = length(y))
    Wbeta <- array(0, dim = length(y))
    S <- array(inputS, dim = length(y) + 1)
    PSI <- array(0, dim = length(y) + 1)
    a <- array(mvec, dim = c(2,input$k +1 ,length(y) + 1))
    R <- array(diag(Cdiag, nrow = 2), dim = c(2,2,input$k +1 ,length(y) + 1))
    f <- array(rep(NA,1), dim = c(input$k,length(y)))
    Q <- array(rep(NA,1), dim = c(input$k,length(y)))
    LowerInterval <- array(rep(NA,1), dim = c(input$k,length(y)))
    UpperInterval <- array(rep(NA,1), dim = c(input$k,length(y)))
    A <- array(rep(NA,2), dim = c(2,length(y)))
    e <- array(rep(NA,1), dim = c(input$k,length(y)))
    P <- array(rep(NA,1), dim = c(input$k,length(y)))
    Qij <- array(NA, dim=c(input$k, input$k, length(y)))
    ActualError <- array(NA, dim = c(input$k,length(y)))
    BayesFactor <- array(NA, dim = length(y))
    LocalBayesFactor <- array(1, dim = length(y) + 1)
    l <- array(0, dim = length(y) + 1)
    Hmu <- array(0, dim = length(y) + 1)
    Hbeta <- array(0, dim = length(y) + 1)
    
    Gvec <- as.numeric(unlist(strsplit(inputG,",")))
    G = matrix(Gvec, nrow = 2, byrow =TRUE)
    FF = as.numeric(unlist(strsplit(inputFF,",")))
    
    W <- array(diag(0, nrow = 2), dim = c(2,2,length(y)))
    H <- array(diag(0, nrow = 2), dim = c(2,2,length(y) + 1))
      
    for(i in 1:length(y)){
      Wbeta[i] <- (1/inputdeltabeta - 1)*C[2,2,i]
      Wmu[i] <- (1/inputdeltamu - 1)*C[1,1,i] 
      W[1,1,i] <- Wmu[i] + Wbeta[i]
      W[1,2,i] <- Wbeta[i]
      W[2,1,i] <- Wbeta[i]
      W[2,2,i] <- Wbeta[i]
      for(j in 1:input$k){
        a[,j+1,i] <- G%*%a[,j,i]
        R[,,j+1,i] <- G%*%R[,,j,i]%*%t(G) + W[,,i] + H[,,i]
        for(t in 1:input$k){
          Qij[j,t,i] <- R[1,1,j+1,i] + (t-j)*R[1,2,j+1,i]
        }
        f[j,i] <- FF%*%a[,j+1,i]
        Q[j,i] <- FF%*%R[,,j+1,i]%*%FF + S[i] + PSI[i]
        P[j,i] <- pt.scaled(log(0.3), df = n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        LowerInterval[j,i] = qt.scaled(0.025, df=n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        UpperInterval[j,i] = qt.scaled(0.975, df=n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        if(UpperInterval[j,i] > log(10)){
          UpperInterval[j,i] = log(10)
        }
      }
      A[,i] <- R[,,2,i]%*%FF/Q[1,i]
      e[,i] <- y[i:(i+input$k -1)] - f[,i]
      BayesFactor[i] <- dt.scaled(e[1,i]/sqrt(Q[1,i]), df = n[i])/dt.scaled(e[1,i]/sqrt(Q[1,i]), df =n[i], sd = 3)
      if(is.na(BayesFactor[i])){
        BayesFactor[i] <- 1 #If y[i] is missing then we do not know anything about the Bayes' factor. But if the previous two Bayes' factors are less than one and then there is a missing observation we do not want the cumulative Bayes' factor to start over incase we miss a change in parameter values
      }
      LocalBayesFactor[i+1] <- BayesFactor[i]*min(1, LocalBayesFactor[i])
      if(LocalBayesFactor[i+1] < 1){
        l[i+1] = l[i] + 1
      }
      if(is.na(e[1,i]) || BayesFactor[i] < exp(-2)){
        n[i+1] <- n[i]
        d[i+1] <- d[i]
        m[,i+1] <- m[,i]
      } else {
        n[i+1] <- inputdel*n[i] + 1
        d[i+1] <- inputdel*d[i] + S[i]*e[1,i]^2/Q[1,i]
        m[,i+1] <- a[,2,i] + A[,i]*e[1,i]
      }
      S[i+1] = d[i+1]/n[i+1]
      C[,,i+1] <- (S[i+1]/S[i])*(R[,,2,i] - A[,i]%*%t(A[,i])*Q[1,i])
      a[,1,i+1] <- m[,i+1]
      R[,,1,i+1] <- C[,,i+1]
      if(LocalBayesFactor[i+1] < exp(-2) || l[i+1] > 2){
        Hbeta[i+1] <- (1/inputdeltabetanew - 1/inputdeltabeta)*C[2,2,i+1]
        Hmu[i+1] <- (1/inputdeltamunew - 1/inputdeltamu)*C[1,1,i+1] 
        H[1,1,i+1] <- Hmu[i+1] + Hbeta[i+1]
        H[1,2,i+1] <- Hbeta[i+1]
        H[2,1,i+1] <- Hbeta[i+1]
        H[2,2,i+1] <- Hbeta[i+1]
        LocalBayesFactor[i+1] = 1
        l[i+1] = 0
      }
      if(BayesFactor[i] < exp(-2)){
        PSI[i+1] <- ((inputdelnew - inputdel)*(d[i] - n[i]*S[i]*e[1,i]^2/Q[1,i]))/((inputdelnew*n[i] + 1)*(inputdel*n[i] + 1))
      }
      ActualError[,i] <- Y[i:(i+input$k -1)] - exp(f[,i]) + 0.1
    }
    
    df = data.frame((input$Observation+1):(input$Observation+input$k),exp(f[,input$Observation+1]) -0.1, exp(LowerInterval[,input$Observation+1]) - 0.1, exp(UpperInterval[,input$Observation+1]) - 0.1, Y[(input$Observation+1):(input$Observation+input$k)], ActualError[(input$Observation+1):(input$Observation+input$k)], format(P[,input$Observation+1], digits = 5))
    colnames(df) = c("Time (Hours)","Forecasts", "Lower Interval", "Upper Interval", "y", "e", "P")
    return(df)
    
  })
  
  output$scat <- renderPlot({
    
    SubsetRenalFailureId = subset(RenalFailure, NewPseudoId == input$Id) 
    SubsetRenalFailureId2 = SubsetRenalFailureId[-1, ]
    PatientIndexId = subset(PatientIndex, NewPseudoId == input$Id)
    SubsetRenalFailureId2$Urine60 = SubsetRenalFailureId2$Urine60/PatientIndexId$Weight
    Y = SubsetRenalFailureId2$Urine60
    Yc = Y + 0.1
    y = log(Yc)
    
    Cdiag <- as.numeric(unlist(strsplit(inputC,",")))
    C <- array(diag(Cdiag, nrow = 2), dim = c(2,2,length(y) + 1))
    mvec <- as.numeric(unlist(strsplit(inputm,",")))
    m <- array(mvec, dim = c(2,length(y) + 1))
    d <- array(inputd, dim = length(y) + 1)
    n <- array(inputn, dim = length(y) + 1)
    Wmu <- array(0, dim = length(y))
    Wbeta <- array(0, dim = length(y))
    S <- array(inputS, dim = length(y) + 1)
    PSI <- array(0, dim = length(y) + 1)
    a <- array(mvec, dim = c(2,input$k +1 ,length(y) + 1))
    R <- array(diag(Cdiag, nrow = 2), dim = c(2,2,input$k +1 ,length(y) + 1))
    f <- array(rep(NA,1), dim = c(input$k,length(y)))
    Q <- array(rep(NA,1), dim = c(input$k,length(y)))
    LowerInterval <- array(rep(NA,1), dim = c(input$k,length(y)))
    UpperInterval <- array(rep(NA,1), dim = c(input$k,length(y)))
    A <- array(rep(NA,2), dim = c(2,length(y)))
    e <- array(rep(NA,1), dim = c(input$k,length(y)))
    P <- array(rep(NA,1), dim = c(input$k,length(y)))
    Qij <- array(NA, dim=c(input$k, input$k, length(y)))
    ActualError <- array(NA, dim = c(input$k,length(y)))
    BayesFactor <- array(NA, dim = length(y))
    LocalBayesFactor <- array(1, dim = length(y) + 1)
    l <- array(0, dim = length(y) + 1)
    Hmu <- array(0, dim = length(y) + 1)
    Hbeta <- array(0, dim = length(y) + 1)
    
    Gvec <- as.numeric(unlist(strsplit(inputG,",")))
    G = matrix(Gvec, nrow = 2, byrow =TRUE)
    FF = as.numeric(unlist(strsplit(inputFF,",")))
    
    W <- array(diag(0, nrow = 2), dim = c(2,2,length(y)))
    H <- array(diag(0, nrow = 2), dim = c(2,2,length(y) + 1))
    
    for(i in 1:length(y)){
      Wbeta[i] <- (1/inputdeltabeta - 1)*C[2,2,i]
      Wmu[i] <- (1/inputdeltamu - 1)*C[1,1,i] 
      W[1,1,i] <- Wmu[i] + Wbeta[i]
      W[1,2,i] <- Wbeta[i]
      W[2,1,i] <- Wbeta[i]
      W[2,2,i] <- Wbeta[i]
      for(j in 1:input$k){
        a[,j+1,i] <- G%*%a[,j,i]
        R[,,j+1,i] <- G%*%R[,,j,i]%*%t(G) + W[,,i] + H[,,i]
        for(t in 1:input$k){
          Qij[j,t,i] <- R[1,1,j+1,i] + (t-j)*R[1,2,j+1,i]
        }
        f[j,i] <- FF%*%a[,j+1,i]
        Q[j,i] <- FF%*%R[,,j+1,i]%*%FF + S[i] + PSI[i]
        P[j,i] <- pt.scaled(log(0.3), df = n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        LowerInterval[j,i] = qt.scaled(0.025, df=n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        UpperInterval[j,i] = qt.scaled(0.975, df=n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        if(UpperInterval[j,i] > log(10)){
          UpperInterval[j,i] = log(10)
        }
      }
      A[,i] <- R[,,2,i]%*%FF/Q[1,i]
      e[,i] <- y[i:(i+input$k -1)] - f[,i]
      BayesFactor[i] <- dt.scaled(e[1,i]/sqrt(Q[1,i]), df = n[i])/dt.scaled(e[1,i]/sqrt(Q[1,i]), df =n[i], sd = 3)
      if(is.na(BayesFactor[i])){
        BayesFactor[i] <- 1 #If y[i] is missing then we do not know anything about the Bayes' factor. But if the previous two Bayes' factors are less than one and then there is a missing observation we do not want the cumulative Bayes' factor to start over incase we miss a change in parameter values
      }
      LocalBayesFactor[i+1] <- BayesFactor[i]*min(1, LocalBayesFactor[i])
      if(LocalBayesFactor[i+1] < 1){
        l[i+1] = l[i] + 1
      }
      if(is.na(e[1,i]) || BayesFactor[i] < exp(-2)){
        n[i+1] <- n[i]
        d[i+1] <- d[i]
        m[,i+1] <- m[,i]
      } else {
        n[i+1] <- inputdel*n[i] + 1
        d[i+1] <- inputdel*d[i] + S[i]*e[1,i]^2/Q[1,i]
        m[,i+1] <- a[,2,i] + A[,i]*e[1,i]
      }
      S[i+1] = d[i+1]/n[i+1]
      C[,,i+1] <- (S[i+1]/S[i])*(R[,,2,i] - A[,i]%*%t(A[,i])*Q[1,i])
      a[,1,i+1] <- m[,i+1]
      R[,,1,i+1] <- C[,,i+1]
      if(LocalBayesFactor[i+1] < exp(-2) || l[i+1] > 2){
        Hbeta[i+1] <- (1/inputdeltabetanew - 1/inputdeltabeta)*C[2,2,i+1]
        Hmu[i+1] <- (1/inputdeltamunew - 1/inputdeltamu)*C[1,1,i+1] 
        H[1,1,i+1] <- Hmu[i+1] + Hbeta[i+1]
        H[1,2,i+1] <- Hbeta[i+1]
        H[2,1,i+1] <- Hbeta[i+1]
        H[2,2,i+1] <- Hbeta[i+1]
        LocalBayesFactor[i+1] = 1
        l[i+1] = 0
      }
      if(BayesFactor[i] < exp(-2)){
        PSI[i+1] <- ((inputdelnew - inputdel)*(d[i] - n[i]*S[i]*e[1,i]^2/Q[1,i]))/((inputdelnew*n[i] + 1)*(inputdel*n[i] + 1))
      }
      ActualError[,i] <- Y[i:(i+input$k -1)] - exp(f[,i]) + 0.1
    }
    
    df = data.frame((input$Observation+1):(input$Observation+input$k),exp(f[,input$Observation+1]) -0.1, exp(LowerInterval[,input$Observation+1]) - 0.1, exp(UpperInterval[,input$Observation+1]) - 0.1, Y[(input$Observation+1):(input$Observation+input$k)], ActualError[(input$Observation+1):(input$Observation+input$k)])
    colnames(df) = c("X","f", "LowerInterval", "UpperInterval", "y", "e")

    dfActualUO = data.frame(1:(input$k+input$Observation), Y[1:(input$k+input$Observation)])
    colnames(dfActualUO) = c("ActualUO", "y")
    
    ggplot()+
      geom_point(data=df, aes(x=X, y =f, colour = "Forecasts"))+
      geom_point(data=dfActualUO, aes(x=ActualUO, y = y, colour = "UO"))+
      geom_point(data=df, aes(x=X, y =LowerInterval), shape=95)+
      geom_errorbar(data = df, aes(x = X, ymin = LowerInterval, ymax = UpperInterval, color = "95% PI"), width = 0)+
      geom_point(data=df, aes(x=X, y =UpperInterval), shape=95)+
      geom_hline(aes(yintercept = 0.5, colour = "AKI1UO"))+
      geom_hline(aes(yintercept = 0.3, colour = "Ralib"))+
      scale_colour_manual(values = c("UO" = "black", "AKI1UO" = "blue", "Ralib" = "red", "95% PI" = "black", "Forecasts" = "red"),
                          guide = guide_legend(override.aes = list(
                            linetype = c('blank', 'solid','blank', 'solid', 'blank'),
                            shape = c(124,NA, 16,NA, 16),
                            size     = c(5, 0.5, 1.5, 0.5, 1.5)),
                            title = ""))+
      labs(title="Urine Output")+
      labs(x = "Time(Hours)")+
      labs(y = "Urine Output (ml/kg)")+
      theme_bw()
    
  })
  
  output$JointProbability <-renderPrint({
    
    SubsetRenalFailureId = subset(RenalFailure, NewPseudoId == input$Id) 
    SubsetRenalFailureId2 = SubsetRenalFailureId[-1, ]
    PatientIndexId = subset(PatientIndex, NewPseudoId == input$Id)
    SubsetRenalFailureId2$Urine60 = SubsetRenalFailureId2$Urine60/PatientIndexId$Weight
    Y = SubsetRenalFailureId2$Urine60
    Yc = Y + 0.1
    y = log(Yc)
    
    Cdiag <- as.numeric(unlist(strsplit(inputC,",")))
    C <- array(diag(Cdiag, nrow = 2), dim = c(2,2,length(y) + 1))
    mvec <- as.numeric(unlist(strsplit(inputm,",")))
    m <- array(mvec, dim = c(2,length(y) + 1))
    d <- array(inputd, dim = length(y) + 1)
    n <- array(inputn, dim = length(y) + 1)
    Wmu <- array(0, dim = length(y))
    Wbeta <- array(0, dim = length(y))
    S <- array(inputS, dim = length(y) + 1)
    PSI <- array(0, dim = length(y) + 1)
    a <- array(mvec, dim = c(2,input$k +1 ,length(y) + 1))
    R <- array(diag(Cdiag, nrow = 2), dim = c(2,2,input$k +1 ,length(y) + 1))
    f <- array(rep(NA,1), dim = c(input$k,length(y)))
    Q <- array(rep(NA,1), dim = c(input$k,length(y)))
    LowerInterval <- array(rep(NA,1), dim = c(input$k,length(y)))
    UpperInterval <- array(rep(NA,1), dim = c(input$k,length(y)))
    A <- array(rep(NA,2), dim = c(2,length(y)))
    e <- array(rep(NA,1), dim = c(input$k,length(y)))
    P <- array(rep(NA,1), dim = c(input$k,length(y)))
    Qij <- array(NA, dim=c(input$k, input$k, length(y)))
    ActualError <- array(NA, dim = c(input$k,length(y)))
    BayesFactor <- array(NA, dim = length(y))
    LocalBayesFactor <- array(1, dim = length(y) + 1)
    l <- array(0, dim = length(y) + 1)
    Hmu <- array(0, dim = length(y) + 1)
    Hbeta <- array(0, dim = length(y) + 1)
    
    Gvec <- as.numeric(unlist(strsplit(inputG,",")))
    G = matrix(Gvec, nrow = 2, byrow =TRUE)
    FF = as.numeric(unlist(strsplit(inputFF,",")))
    
    W <- array(diag(0, nrow = 2), dim = c(2,2,length(y)))
    H <- array(diag(0, nrow = 2), dim = c(2,2,length(y) + 1))
    
    for(i in 1:length(y)){
      Wbeta[i] <- (1/inputdeltabeta - 1)*C[2,2,i]
      Wmu[i] <- (1/inputdeltamu - 1)*C[1,1,i] 
      W[1,1,i] <- Wmu[i] + Wbeta[i]
      W[1,2,i] <- Wbeta[i]
      W[2,1,i] <- Wbeta[i]
      W[2,2,i] <- Wbeta[i]
      for(j in 1:input$k){
        a[,j+1,i] <- G%*%a[,j,i]
        R[,,j+1,i] <- G%*%R[,,j,i]%*%t(G) + W[,,i] + H[,,i]
        for(t in 1:input$k){
          Qij[j,t,i] <- R[1,1,j+1,i] + (t-j)*R[1,2,j+1,i]
        }
        f[j,i] <- FF%*%a[,j+1,i]
        Q[j,i] <- FF%*%R[,,j+1,i]%*%FF + S[i] + PSI[i]
        P[j,i] <- pt.scaled(log(0.3), df = n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        LowerInterval[j,i] = qt.scaled(0.025, df=n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        UpperInterval[j,i] = qt.scaled(0.975, df=n[i], mean = f[j,i], sd = sqrt(Q[j,i]))
        if(UpperInterval[j,i] > log(10)){
          UpperInterval[j,i] = log(10)
        }
      }
      A[,i] <- R[,,2,i]%*%FF/Q[1,i]
      e[,i] <- y[i:(i+input$k -1)] - f[,i]
      BayesFactor[i] <- dt.scaled(e[1,i]/sqrt(Q[1,i]), df = n[i])/dt.scaled(e[1,i]/sqrt(Q[1,i]), df =n[i], sd = 3)
      if(is.na(BayesFactor[i])){
        BayesFactor[i] <- 1 #If y[i] is missing then we do not know anything about the Bayes' factor. But if the previous two Bayes' factors are less than one and then there is a missing observation we do not want the cumulative Bayes' factor to start over incase we miss a change in parameter values
      }
      LocalBayesFactor[i+1] <- BayesFactor[i]*min(1, LocalBayesFactor[i])
      if(LocalBayesFactor[i+1] < 1){
        l[i+1] = l[i] + 1
      }
      if(is.na(e[1,i]) || BayesFactor[i] < exp(-2)){
        n[i+1] <- n[i]
        d[i+1] <- d[i]
        m[,i+1] <- m[,i]
      } else {
        n[i+1] <- inputdel*n[i] + 1
        d[i+1] <- inputdel*d[i] + S[i]*e[1,i]^2/Q[1,i]
        m[,i+1] <- a[,2,i] + A[,i]*e[1,i]
      }
      S[i+1] = d[i+1]/n[i+1]
      C[,,i+1] <- (S[i+1]/S[i])*(R[,,2,i] - A[,i]%*%t(A[,i])*Q[1,i])
      a[,1,i+1] <- m[,i+1]
      R[,,1,i+1] <- C[,,i+1]
      if(LocalBayesFactor[i+1] < exp(-2) || l[i+1] > 2){
        Hbeta[i+1] <- (1/inputdeltabetanew - 1/inputdeltabeta)*C[2,2,i+1]
        Hmu[i+1] <- (1/inputdeltamunew - 1/inputdeltamu)*C[1,1,i+1] 
        H[1,1,i+1] <- Hmu[i+1] + Hbeta[i+1]
        H[1,2,i+1] <- Hbeta[i+1]
        H[2,1,i+1] <- Hbeta[i+1]
        H[2,2,i+1] <- Hbeta[i+1]
        LocalBayesFactor[i+1] = 1
        l[i+1] = 0
      }
      if(BayesFactor[i] < exp(-2)){
        PSI[i+1] <- ((inputdelnew - inputdel)*(d[i] - n[i]*S[i]*e[1,i]^2/Q[1,i]))/((inputdelnew*n[i] + 1)*(inputdel*n[i] + 1))
      }
      ActualError[,i] <- Y[i:(i+input$k -1)] - exp(f[,i]) + 0.1
    }
    
    JointProbability = rep(0,length(y))
    
    for(i in 1:length(y)){
      Joint = Qij[,,i]
      diag(Joint) = Q[,i]
      Joint[lower.tri(Joint)] <- t(Joint)[lower.tri(Joint)]
      set.seed(2131)
      JointProbability[i] = pmvt(lower = -Inf*rep(1,input$k), upper = rep(log(0.3),input$k), sigma = Joint, df = n[i], delta = f[,i], type = "shifted")[1]
    }
    
    JointProbability = round(JointProbability, 6)
    ProlongedHighRisk = rep(0,(length(y) + 1))

    for(i in 2:length(y)){
      if(JointProbability[i] > 0.8 & ProlongedHighRisk[i-1] == 0){
        ProlongedHighRisk[i] = 1
      } else if(JointProbability[i] > 0.8 & ProlongedHighRisk[i-1] != 0)
        ProlongedHighRisk[i] = ProlongedHighRisk[i-1] + 1
    }
    
    lst=list()
    if(JointProbability[input$Observation+1] >0.8){
      lst[[1]] = paste("The probability that the next", input$k, "urine outputs are all below 0.3 is", JointProbability[input$Observation + 1])
      lst[[2]] = paste("This patient has been at high risk for", ProlongedHighRisk[input$Observation + 1] - 1, "hour")
      if(!is.na(PatientIndexId$FirstFilter)){
        lst[[3]] = paste("Decision for RRT made at hour",which(SubsetRenalFailureId2$Time == as.character(PatientIndexId$FirstFilter)))
      }
      return(lst)
    } else {
      lst[[1]] = paste("The probability that the next", input$k, "urine outputs are all below 0.3 is", JointProbability[input$Observation + 1])
      if(!is.na(PatientIndexId$FirstFilter)){
        lst[[2]] = paste("Decision for RRT made at hour",which(SubsetRenalFailureId2$Time == as.character(PatientIndexId$FirstFilter)))
      }
      return(lst)
    }
    
  })
  
  output$Gender = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Gender: ", PatientIndex$Gender[i])
  })
  
  output$Age = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Age: ", PatientIndex$Age[i])
  })
  
  output$Weight = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Weight(kg): ", PatientIndex$Weight[i])
  })
  
  output$Height = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Height(cm): ", PatientIndex$Height[i])
  })
  
  output$ProcDetails = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("ProcDetails: ", PatientIndex$ProcDetails[i])
  })
  
  output$Urgency = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Urgency: ", PatientIndex$Urgency[i])
  })
  
  output$scat3 <- renderPlotly({
    
    SubsetFlowSheeti = subset(FlowSheet, NewPseudoId == input$Id)
    Difftime = c(1, as.numeric(difftime(SubsetFlowSheeti$Result.DT[2:length(SubsetFlowSheeti$Result.DT)], SubsetFlowSheeti$Result.DT[1:(length(SubsetFlowSheeti$Result.DT)-1)], units="hours")))
    Times = cumsum(Difftime)
    RealTimes = Times[which(Times<=(input$Observation+1))]
    RealTimes = RealTimes - 1
    SubsetFlowSheetiRealTime = SubsetFlowSheeti[1:length(which(Times<=(input$Observation+1))),]
    
    ggplot(data=SubsetFlowSheetiRealTime, aes(x = RealTimes))+
      geom_point(aes(y = as.numeric(as.character(CVP)), colour = "CVP"))+
      scale_colour_manual("", 
                          breaks = c("CVP"),
                          values = c("black"))+
      labs(title="CVP")+
      labs(x = "Time(hours)")+
      labs(y = "CVP")+
      theme_bw()
    
  })
  
  output$scat4 <- renderPlotly({
    
    Biochemistryi = subset(Biochemistry, NewPseudoId == input$Id)
    Biochemistryi = Biochemistryi[order(Biochemistryi$PostOpUsandEsTime),]
    
    FirstUreaCreatinineTime = as.character(Biochemistryi$PostOpUsandEsTime[1]) #This will need to be used to find the corresponding time and hence position in the RenalFailure table. So that the first Urea and Creatine are plotted at time, say 6, if the first measurements are made at the sixth hour of UOs for that person
    SubsetRenalFailurei = subset(RenalFailure, NewPseudoId == input$Id)
    FirstUreaCreatininePosition = which(SubsetRenalFailurei$Time == FirstUreaCreatinineTime)
    
    if(length(FirstUreaCreatininePosition)>0){
      FirstUreaCreatininePosition = FirstUreaCreatininePosition
    } else {
      FirstUreaCreatininePosition = as.numeric(difftime(FirstUreaCreatinineTime, SubsetRenalFailurei$Time[1], units="hours")) + 1
    }
    
    if(length(Biochemistryi$PostOpUsandEsTime)>1){
      Difftime = c(FirstUreaCreatininePosition, as.numeric(difftime(Biochemistryi$PostOpUsandEsTime[2:length(Biochemistryi$PostOpUsandEsTime)], Biochemistryi$PostOpUsandEsTime[1:(length(Biochemistryi$PostOpUsandEsTime)-1)], units="hours")))
    } else {
      Difftime = FirstUreaCreatininePosition
    }
    
    Times = cumsum(Difftime)
    RealTimes = Times[which(Times<=(input$Observation+1))]
    RealTimes = RealTimes - 1
    BiochemistryiRealTime = Biochemistryi[1:length(which(Times<=(input$Observation+1))),]
    
    if(length(RealTimes)>0){
      ggplot(data=BiochemistryiRealTime, aes(x = RealTimes))+
        geom_point(aes(y=Urea, colour = "Urea"))+
        geom_point(aes(y=Creatinine, colour = "Creatinine"))+
        scale_colour_manual("", 
                            breaks = c("Urea", "Creatinine"),
                            values = c("red", "black"))+
        labs(title="Urea and Creatinine")+
        labs(x = "Time(hours)")+
        labs(y = "Urea and Creatinine")+
        theme_bw()
    } else {
      ggplot()+theme_bw()
    }
    
  })
  
  output$scat6 <- renderPlotly({
    
    SubsetFlowSheeti = subset(FlowSheet, NewPseudoId == input$Id)
    Difftime = c(1, as.numeric(difftime(SubsetFlowSheeti$Result.DT[2:length(SubsetFlowSheeti$Result.DT)], SubsetFlowSheeti$Result.DT[1:(length(SubsetFlowSheeti$Result.DT)-1)], units="hours")))
    Times = cumsum(Difftime)
    Replace = which(is.na(SubsetFlowSheeti$ReliableART.M))
    SubsetFlowSheeti$ReliableART.M[Replace] = as.numeric(as.character(SubsetFlowSheeti$NBP.M[Replace]))
    RealTimes = Times[which(Times<=(input$Observation+1))]
    RealTimes = RealTimes - 1
    SubsetFlowSheetiRealTime = SubsetFlowSheeti[1:length(which(Times<=(input$Observation+1))),]
    
    ggplot(data=SubsetFlowSheetiRealTime, aes(x = RealTimes))+
      geom_point(aes(y=ReliableART.M, colour = "ReliableART.M"))+
      scale_colour_manual("", 
                          breaks = c("ReliableART.M"),
                          values = c("black"))+
      labs(title="ReliableART.M")+
      labs(x = "Time(hours)")+
      labs(y = "ReliableART.M")+
      theme_bw()
    
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)