# A script to load appropriate libraries and settings.
# Load in libraries

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Usage example
packages<-c("ggplot2",
            "plotly", 
            "peakPick", 
            "signal", 
            "smoother",
            "pracma",
            "pdist",
            "abind",
            "wavelets",
            "FNN",
            "zoo",
            "deSolve",
            "fields",
            "grid",
            "gridExtra",
            "R.matlab",
            "lhs",
            "lattice",
            "moments",
            "mvtnorm",
            "mclust",
            "MASS",
            "latex2exp",
            "ggpubr",
            "ks",
            "plot3D",
            "entropy",
            "dplyr",
            'miscTools',
            'rmutil',
            'matrixStats',
            'gridExtra',
            'GGally',
            'lattice',
            'pROC',
            'mice',
            'shiny',
            'xtable')
check.packages(packages)

library( ggplot2 )#
library( plotly )#
library( peakPick )#
library( signal )#
library( smoother )#
library( pdist )#
library( abind )#
library( wavelets )#
library( FNN )#
library( zoo )#
library( deSolve )#
library( fields )#
library( grid )#
library( gridExtra )#
library( R.matlab )#
#library( sn )
library( lhs )#
library( lattice )#
library( moments )#
library( mvtnorm )#
library( mclust )#
library( MASS )#
library( latex2exp )#
library( ggpubr )#
library( ks )#
library( plot3D )#
library( miscTools )
library( matrixStats )
library(gridExtra)
library(rmutil)
library(GGally)
library(lattice)
library(pROC)
library(mice)
library(shiny)
library(xtable)
# Load Data
options(digits = 15)
options(digits.secs=3)

source( "FourierSourceFunctions.R" )
source( "WaveSplitSourceFunctions.R" )
source( "PeakExtractionSourceFunctions.R" )
source( 'EmulationSourceCode.R' )
source( 'BayesLinearDynamicUpdateSourceCode.R' )
source( 'DiscreteAnalysisSourceFunctions.R' )
source( 'ASWFSourceFunctions.R' )
source( 'AccessandLoadingDataSourceFunctions.R' )
source( 'mitdbSourceFunctions.R' )
source( 'AFDetectionSourceFunctions.R')
source( 'ECGSimulatorSourceFunctions.R')
source( 'HistoryMatchingSourceFunctions.R')
source( 'PlottingSourceFunctions.R' )
source( 'VectorisedOperationsSourceFunctions.R' )
source( 'DenistyFunctionsSourceFunctions.R' )
source( 'DenistyFunctionsSourceFunctions.R' )
source( 'BayesClassifierSourceFunctions.R' )
source( 'KernelDensityEstimationSourceFunctions.R' )
source( 'CorrelationDiscrepancySourceFunctions.R' )
source('BayesLinearUpdateSourceFunctions.R' )
source('PWaveHistoryMatchingSourceFunctions.R' )
source('BayesLinearBayesForecastingSourceFunctions.R')
source('PreOpModelSourceFunctions.R')
source('ForwardModellingSourceFunctions.R')
source('CriticalValueEmulationsSourceFunctions.R')
source('StepWiseModelSelectionSourceFunctions.R')
source('PriorElicitationRhythmSourceFunctions.R')
print('Libraries, data options and source functions loaded.')

