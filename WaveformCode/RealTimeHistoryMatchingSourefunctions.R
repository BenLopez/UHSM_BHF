{Sigma = matrix( c(1 , 0 , 0 , 1) , 2 , 2 )
q <- 0.95

CTEm_VolumeHyperSphere <- function(n,r){
  return((pi^(n/2))/(gamma((n/2) + 1))*r^n)
}

CTEm_VolumeofElipse <- function(n , r , a  ){
  # Volume on ball radius c where a = c vec(b)
  V = CTEm_VolumeHyperSphere(n,r)*cumprod(a)[n]
  return(V)
}
CTEm_CredibleVolume <- function(Sigma , q){
  
  n <- dim(Sigma)[1]
  r <- sqrt(qchisq(q , n))
  a <- eigen(Sigma)$Values
  V = CTEm_VolumeofElipse(n , r , a )
  return(V)
}
CTEm_CredibleVolume(Sigma , 0.95)
RTHM_EstimateNumberofPoints <- function(x , ztest , fx ,V_e , PwaveEmulatorParameters = BE_CreateDefaultEmulationClass(), doplots = 0,  qq = 0.05){
  
  UBztest <- ztest + 3*sqrt(V_e)
  LBztest <- ztest - 3*sqrt(V_e)
  
  if(sum(abs(fx - LBztest)<0.0075)==0 & sum(abs(fx - UBztest)<0.0075) ==0){
    indexs <- sort(c((which.min(abs(fx - LBztest))),(which.min(abs(fx - UBztest)) )))
  }
  if(sum(abs(fx - LBztest)<0.0075)==0 & sum(abs(fx - UBztest)<0.0075) != 0){
    indexs <- sort(c(min(which(abs(fx - UBztest)<0.0075)),max(which(abs(fx - UBztest)<0.0075) )))
  }
  if(sum(abs(fx - LBztest)<0.0075)!=0 & sum(abs(fx - UBztest)<0.0075) == 0){
    indexs <- sort(c(min(which(abs(fx - LBztest)<0.0075)),max(which(abs(fx - LBztest)<0.0075) )))
  }
  if(sum(abs(fx - LBztest)<0.0075)!=0 & sum(abs(fx - UBztest)<0.0075) != 0){
    indexs <- sort(c(min(which(abs(fx - LBztest)<0.0075)),max(which(abs(fx - UBztest)<0.0075) )))
  }
  
  index1 <- indexs[1]
  index2 <- indexs[2]   
  
  rangex <- sort(c(x[index1] , x[index2] ))
  
  muxstar <- mean(rangex)
  sigmaxstar <- (abs(diff(rangex))/6)^2
  if(doplots == 1){
    pa <- ggplot(data.frame(x , fx=rollmean(fx,40,na.pad = T)) ) +
      geom_line( aes(x , fx, color="F(x)") ) +
      xlab(TeX('x')) +
      ylab(TeX('z')) + 
      ggtitle(TeX('Finding C[x]'))  +
      geom_line(data = data.frame(x, z = (ztest+0*fx)) , aes(x,z,color = 'z') ) +
      geom_line(data = data.frame(x, z= ztest + 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) +
      geom_line(data = data.frame(x, z= ztest - 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) +
      geom_line(data = data.frame(x= 0*x + x[index1]  , z = fx),aes(x,z,color = 'C[x]'))+
      geom_line(data = data.frame(x= 0*x + x[index2]  , z = fx),aes(x,z,color = 'C[x]'))+
      scale_color_manual(values = c(
        "F(x)" = 'blue',
        'z' = 'black',
        'Credible Interval z' = 'red',
        'C[x]' = 'lightblue'  )) +
      labs(color = 'Legend')
  }
  
  ftest <- rollmean(fx,40,na.pad = T)[(x>rangex[1])&(x<rangex[2])] 
  tmp <- fx[(x>rangex[1])&(x<rangex[2])] 
  ftest[is.na(  ftest)] <- tmp[is.na(  ftest)]
  xtest <- x[(x>rangex[1])&(x<rangex[2])] 
  
  if(doplots == 1){
    pb <- ggplot(data.frame(x = xtest , fx = ftest) ) +
      geom_line( aes(x , fx, color="F(x)") ) +
      xlab(TeX('x')) +
      ylab(TeX('z')) + 
      ggtitle(TeX('Final Wave'))  +
      geom_line(data = data.frame(x = xtest, z = (ztest+0*ftest)) , aes(x,z,color = 'z') ) +
      geom_line(data = data.frame(x = xtest, z= ztest + 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) +
      geom_line(data = data.frame(x = xtest, z= ztest - 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) +
      geom_line(data = data.frame(x = 0*xtest + x[index1]  , z = ftest),aes(x,z,color = 'C[x]'))+
      geom_line(data = data.frame(x = 0*xtest + x[index2]  , z = ftest),aes(x,z,color = 'C[x]'))+
      scale_color_manual(values = c(
        "F(x)" = 'blue',
        'z' = 'black',
        'Credible Interval z' = 'red',
        'C[x]' = 'lightblue'  )) +
      labs(color = 'Legend')
  }
  
  H <- PwaveEmulatorParameters$MeanFunction(xtest)
  if(dim(H)[1]<dim(H)[2]){return(list(nproposed=dim(H)[2] + 2,muxstar=muxstar,sigmaxstar=sigmaxstar ))}
  Beta = solve(DP_AddNugget(t(H)%*%H , 0.0000000001*diag(dim(H)[2]) ))%*%t(H)%*%ftest
  sigma = as.numeric(var(ftest - H%*%Beta))
  nproposed <- dim(H)[2] + 2
  
  if(length(xtest)<nproposed){
    return(list(nproposed=nproposed,muxstar=muxstar,sigmaxstar=sigmaxstar ))
  }
  xtmp <- x[round(as.matrix(seq(index1 , index2 , abs(index2 - index1)/nproposed)))]
  M_VW <- as.numeric(sigma)
  
  
  while( (M_VW/V_e) > qq){
    
    if(length(xtest) < nproposed){
      disp('Not enough points')
      break}
      
    xtmp <- as.matrix(x[round(as.matrix(seq(index1 , index2 , abs(index2 - index1)/nproposed)))][1:nproposed])
    KXX <-  sigma*CF_ExponentialFamily( xtmp , xtmp , PwaveEmulatorParameters$CorrelationLength(xtmp, nproposed) , p=2) 
    KXstarX <- sigma*CF_ExponentialFamily( xtest , xtmp , PwaveEmulatorParameters$CorrelationLength(xtmp, nproposed) , p=2) 
    KXstarXstar <- sigma*CF_ExponentialFamily( xtest , xtest , PwaveEmulatorParameters$CorrelationLength(xtmp, nproposed) , p=2)
    
#   if(kappa(DP_AddNugget(KXX , 1e-5*sigma*diag(dim(KXX)[1]))) > 1e16){break}
    V_D <- KXstarXstar - KXstarX%*%solve(DP_AddNugget(KXX , 1e-5*sigma*diag(dim(KXX)[1])) )%*%t(KXstarX)
#   if(kappa(V_D) > 1e16){break}
    M_VW <- quantile(diag(V_D),0.8)
    nproposed <- nproposed + 1
  }
  
  if(doplots == 1){
    if(nproposed > (dim(H)[2] + 2) ){
      nproposed <- nproposed -1
    }
    pb <- pb + 
      geom_point(data = data.frame(x= xtmp[1:nproposed] , z = rollmean(fx,40,na.pad = T)[round(as.matrix(seq(index1 , index2 , abs(index2 - index1)/(nproposed)) )[1:(nproposed)])]),
                 aes(x , z , color = 'Training Points in C[x]'))+
      scale_color_manual(values = c(
        "F(x)" = 'blue',
        'z' = 'black',
        'Credible Interval z' = 'red',
        'C[x]' = 'lightblue',
        'Training Points in C[x]' = 'purple')) +
      labs(color = 'Legend')
      }
  if(doplots ==1){
    x11(20,14)
  grid.arrange(pa,pb,nrow = 2 , ncol = 1)
  return( list(nproposed=(nproposed),muxstar=muxstar,sigmaxstar=sigmaxstar,pa=pa,pb=pb ) )}else{
    return( list(nproposed=(nproposed),muxstar=muxstar,sigmaxstar=sigmaxstar) )  }
}
}
