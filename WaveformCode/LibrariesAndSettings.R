# A script to load appropriate libraries and settings.
# Load in libraries

library(ggplot2)
library(plotly)
library(peakPick)
library(signal)
library(smoother)
library(pdist)
library(abind)
library(wavelets)
library(FNN)
library(zoo)
library(deSolve)
library(fields)
library(grid)
library(gridExtra)

# Load Data
options(digits = 15)
options(digits.secs=3)

source("FourierSourceFunctions.R")
source("SourceFunctions.R")
source("PeakExtractionSourceCode.R")
source('EmulationSourceCode.R')
source('BayesLinearDynamicUpdateSourceCode.R')
source('DiscreteAnalysisSourceFunctions.R')
print('Libraries, data options and source functions loaded.')

