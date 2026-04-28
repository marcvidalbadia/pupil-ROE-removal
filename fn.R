pup.med <- function(y, ant=0.1, post=0.2, sp=30, method=c("t-Student","Gaussian","Kalman")) {
  ##' pup.med - R function to remove blink artifacts and impute missing data
  ##' y: a pupillary time series containing blink artifacts
  ##' ant: time before artifact onset (sec.)
  ##' post: time after artifact onset (sec.)
  ##' sp: sampling rate of x
  
  if (!(inherits(y, "numeric")))
    stop("Argument y not a numeric object")
  l <- length(y)
  original <- y
  if(y[1]<2) { y[1] <- median(y[1:sp*2], na.rm=T)
  cat("Missing first observation - median imputation") }
  
  get.intv <- function(outlier, l.=l, ant.=ant, post.=post, sp.=sp) {
    if (outlier<(l.-sp.*(ant.+post.)*2) & outlier>(sp.*(ant.+post.)*2)) {
      intv <- (outlier-sp.*ant.):(outlier+sp.*post.)
    } else { intv <- outlier }
    return(intv)
  }
  
  if (length(which(y<3))!=0) {
    outliers <- attributes(
      imputeFin::impute_AR1_Gaussian(y,remove_outliers = T, verbose=F, 
                                     outlier_prob_th = 0.05))$index_outliers
    for (indx in 1:length(outliers))  y[get.intv(outliers[indx])] <- NA
    missing <- which(is.na(y))
  } else { 
    outliers <- 0
    missing <- which(is.na(y))
  }
  if (method=="t-Student") {
    y <- imputeFin::impute_AR1_t(y,remove_outliers = F, verbose=F)
  } else if (method=="Gaussian") {
    y <- imputeFin::impute_AR1_Gaussian(y,remove_outliers = F, verbose=F)
  } else if (method=="Kalman") {
    y <-  imputeTS::na_kalman(y)
  }
  if(is.na(y[l])) {
    cat("Missing last observation - median imputation")
    y[l] <- median(y[(l-sp*2):l], na.rm=T);
    if (method=="t-Student") {
      y <- imputeFin::impute_AR1_t(y,remove_outliers = F, verbose=F)
    } else if (method=="Gaussian") {
      y <- imputeFin::impute_AR1_Gaussian(y,remove_outliers = F, verbose=F)
    } else if (method=="Kalman") {
      y <-  imputeTS::na_kalman(y)
    }
  }
  
  # --- old
  gap_sec <- 0.2   # 200 ms separation between blinks
  nblinks <- outliers[c(1,which(diff(outliers)>sp * gap_sec)+1)]
  ipblinks <- missing[c(1,which(diff(missing)>sp * gap_sec)+1)]
  # --- old
  # --- Blink detection (separate) ---     depends on signal only
  thr <- quantile(original, 0.05, na.rm = TRUE)
  dy  <- c(0, diff(original))
  blink_idx <- which(original < thr | dy < -3 * sd(dy, na.rm = TRUE))
  ipblinks <- blink_idx[c(1, which(diff(blink_idx) > sp * 0.2) + 1)]
  # --- Blink detection (separate) ---     depends on signal only
  b <- rep(0,l)
  b[ipblinks] <- 1
  ratey <- c()
  #sy <- seq(0,l,sp)
  sy <- seq(1, l, sp)
  for (k in 1:(length(sy)-1)) ratey <- append(ratey,length(which(b[sy[k]:sy[k+1]]==1)))
  blink.rate <- mean(ratey)
  
  attributes(y) <- NULL
  pupmed <- list(original, y,outliers,missing,nblinks,blink.rate)
  names(pupmed) <- c("Originaldata","Pupildata","Outliers","Missing","Estimated_blinks","Blink_rate")
  return(pupmed)
}

which.peaks = function(x,partial=FALSE,decreasing=FALSE){
  if (decreasing){
    if (partial){
      which(diff(c(FALSE,diff(x)>0,TRUE))>0)
    }else {
      which(diff(diff(x)>0)>0)+1
    }
  } else {
    if (partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    } else {
      which(diff(diff(x)>=0)<0)+1
    }
  }
}



artifact.corrector <- function(x, W, sd.factor=3) {
  ##' x: reconstructed pupil time series with pup.med function
  ##' W: critical frequencies of the filter. See signal::butter function
  ##' sd.factor: hyperparameter to detect artifact onset
  
  bf <- signal::butter(3, W=W, type="pass")
  y1 <- signal::filtfilt(bf,x)
  dy1 <- diff(y1)
  
  minimas <- which.peaks(dy1, partial = FALSE, decreasing = TRUE)
  maximas <- which.peaks(dy1, partial = FALSE, decreasing = FALSE)
  
  vmin <- rep(NA,length(dy1))
  vmin[minimas] <- dy1[minimas] 
  
  minart <- which(vmin < median(dy1) - sd.factor * sd(dy1)) # detected artifacts
  
  for (sm in 1:length(minart)) {
    maximas <- which.peaks(y1, partial = FALSE, decreasing = FALSE)
    minimas <- which.peaks(y1, partial = FALSE, decreasing = TRUE)
    
    sminart <- minimas[which(minimas > minart[sm])][1] 
    
    if(length(sminart)!=0 & !is.na(sminart)){
      grid <- sort(c(maximas,minimas))
      
      pre  <- grid[which(grid==sminart)-1]
      post <- grid[which(grid==sminart)+1]
      
      if(length(pre)==0 || is.na(pre)) pre <- 1
      if(length(post)==0 || is.na(post)) post <- length(x)  # <-- FIX
      
      molne1 <- which(y1[sminart:post] > y1[pre])[1]
      if (!is.na(molne1)) post <- sminart + molne1
      
      molne2 <- which(y1[pre:sminart] < y1[post])[1]
      if (!is.na(molne2)) pre <- pre + molne2
      
      int <- pre:post
      n <- length(int)
      
      if(!is.na(y1[post]) & !is.na(y1[pre])) {
        if(y1[post] > y1[pre]) {
          x[int] <- (x[int]-x[int[n]]) -(y1[int]-y1[int[n]]) + x[int[n]]
        } else {
          x[int] <- (x[int]-x[int[1]]) -(y1[int]-y1[int[1]]) + x[int[1]]
        }
      }
    }
  }
  
  attr(x, 'Number of corrected artifacts') <- length(minart)
  attr(x, 'Artifact onsets') <- minart
  
  return(x)
}


pup.artifact <- function(y,
                         sd.factor.low=3,
                         sd.factor.high=3,
                         Nf=15,
                         LPF=NA) {
  ##' y: reconstructed pupil time series with pup.med function
  ##' sd.factor.low: hyperparameter to detect low frequency artifacts
  ##' sd.factor.high: hyperparameter to detect high frequency artifacts
  ##' Nf: Nyquist frequency 
  ##' LPF: final smoothing
  
  if (!(inherits(y, "numeric")))
    stop("Argument y not a numeric object")
  
  N <- length(y)
  y <- y - (my <- mean(y))
  
  # Low frequency
  art.onset.low <- NA
  if (!is.na(sd.factor.low)) {
    art.onset.low <- c()
    for (turbi in seq(0.03/Nf,4/Nf,0.015/Nf)) {
      y <- artifact.corrector(y[1:N], W=c(0,turbi), sd.factor=sd.factor.low)
      if (length(attributes(y)$`Artifact onsets`)!=0) 
        art.onset.low <- append(attributes(y)$`Artifact onsets`, art.onset.low)
    }
  }
  
  # High frequency
  art.onset.high <- NA
  if (!is.na(sd.factor.high)) {
    y <- y - (my2 <- mean(y))
    art.onset.high <- c(0)
    
    for (turbi in seq(0.5/Nf,4/Nf,0.015/Nf)) {
      y <- artifact.corrector(y[1:N], W=c(0.25/Nf,turbi), sd.factor=sd.factor.high)
      if (length(attributes(y)$`Artifact onsets`)!=0) 
        art.onset.high <- append(attributes(y)$`Artifact onsets`, art.onset.high)
    }
    
    y <- y + my + my2
  } else { 
    y <- y + my 
  }
  
  # Smoothing
  if (!is.na(LPF)) {
    bf <- signal::butter(3, c(0,LPF/Nf), type="pass")
    y <- y - (my <- mean(y))
    y <- signal::filtfilt(bf,y) + my
  } 
  
  attr(y, "Artifact onsets low freq.") <- art.onset.low 
  attr(y, "Artifact onsets high freq.") <- art.onset.high
  
  return(y)
}

