# A script to load appropriate libraries and settings.
# Load in libraries

library(ggplot2)
library(plotly)
library(peakPick)
library(signal)
library(smoother)
library(pdist)
library(abind)

# Load Data
options(digits = 15)
options(digits.secs=3)

source("FourierSourceFunctions.R")
source("SourceFunctions.R")
source("PeakExtractionSourceCode.R")
source('EmulationSourceCode.R')

print('Libraries, data options and source functions loaded.')

