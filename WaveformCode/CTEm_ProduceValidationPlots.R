{Xstar <- Validation_X[ , -3]
EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = EmulatorParametersMean  )

V_std <- var((Validation_ImCrit[,1] - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) ))
R_sq <- (1- var(Validation_ImCrit[,1] - EmulatorOutput$E_D_fX)/var(Validation_ImCrit[,1]))


BE_PlotStdResiduals(Validation_ImCrit[,1] , EmulatorOutput)
}
{
Xstar <- Validation_X[ , -3]
EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = EmulatorParametersMax  )

V_std <- var((Validation_ImCrit[,2] - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) ))
R_sq <- (1- var(Validation_ImCrit[,2] - EmulatorOutput$E_D_fX)/var(Validation_ImCrit[,2]))


BE_PlotStdResiduals(Validation_ImCrit[,2] , EmulatorOutput) + ggtitle(paste0('R Squared = ' , round(R_sq , 2) , ' Variance Std Resid = ' , round(V_std , 2) ) )}
