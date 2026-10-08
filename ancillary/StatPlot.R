StatPlot <- function(dataframe,mainheading,ylabel,responseno,ylow,yhigh,xlow,xhigh,lwd_scale=1) {
  
  ########################## StatPlot ###############################
  #
  # Script to calculate summary statistics, standardise a response variable from
  # a single brain region and plot it against preference score with adjustable
  # x and y axes.
  #  
  # The variables in the input dataframe "dataframe" must have columns 1 through 16
  # labelled as follows, as in accompanying dataframe  "LeftIMMTest":
  #  1          TrCond
  #  2          TrUntr
  #  3          Region
  #  4            Side
  #  5           Batch
  #  6           Chick
  #  7          Sample
  #  8           TrApp
  #  9           FApp1
  #  10          NApp1
  #  11          NApp2
  #  12          FApp2
  #  13           FApp
  #  14           NApp
  #  15         TotApp
  #  16           Pref
  #  
  # Arguments:
  # dataframe: dataframe from single brain region (left IMM, right IMM, left PPN, right PPN)
  # mainheading: heading for graph in double quotes ("")
  # ylabel: label for y axis in double quotes ("")
  # responseno: index number of column containing the variable to be analysed.
  #   (17 in the accompanying example)
  # ylow: lower bound of y-axis (this may have to be adjusted after a pilot plot)
  # yhigh: upper bound of y-axis (ditto)
  # xlow: lowest x coordinate to be plotted
  # xhigh: highest x coordinate to be plotted
  # lwd_scale: factor applied to all line widths (added for the figures of
  #   CHARM-Vis; default 1)
  #
  # Example console command to analyse the accompanying dataframe LeftIMMTest:
  # StatPlot1(LeftIMMTest,"Left IMM","M-CPEB-3",17,0.6,1.6,50,100)
  #
  # The position of the "Untrained" label in the plot may need to be adjusted (with variable "RelPosn"
  # below) depending on the position of the mean value for untrained chicks.
  # 
  # Output is sent to the console and default plotting device.
  #  
  #  
  # Brian McCabe
  # 27 June 2026
  
  # Load library for lme function
  require(nlme)
  
  # Copy input dataframe to working dataframe and append response variable
  All <- dataframe
  All$UnCorrected <- All[[responseno]]
  
  # Find maximum and minimum preference scores
  MaxPref <- max(All$Pref, na.rm=TRUE)
  MinPref <- min(All$Pref, na.rm=TRUE)
  
  lowxforint <- ifelse(MinPref > 50, 49.8, MinPref)    # Lower bound for regression plot
  highxforint <- MaxPref + 0.2                         # Upper bound for regression plot
  
  # Fit Preference score, extract slope and calculate y-intercepts
  linfit <- lme(UnCorrected ~ Pref, random=~1|Batch, na.action=na.exclude, data=All)
  slope <- summary(linfit)$tTable[2,1]
  Int <- summary(linfit)$tTable[1,1]
  Int50 <- summary(linfit)$tTable[1,1] + slope*50
  IntMaxPref <- summary(linfit)$tTable[1,1] + slope*MaxPref
  
  # Calculate correlation coefficient and variance of residuals
  t <- summary(linfit)$tTable[2,4]               # t statistic
  resdf <- summary(linfit)$tTable[2,3]           # Residual degrees of freedom
  r <- t/sqrt(resdf+t^2)                         # Correlation coefficient
  tProb <- summary(linfit)$tTable[2,5]           # Probability of t
  varFitres <- as.numeric(VarCorr(linfit)[2,1])  # Variance of residuals from regression
  
  # Calculate standard errors of intercepts at 50 and maximum preference score  
  All$PrefLess50 <- All$Pref - 50
  linfitLess50 <- lme(UnCorrected ~ PrefLess50, random=~1|Batch, na.action=na.exclude, data=All)
  seInt50 <- summary(linfitLess50)$tTable[1,2]
  
  All$PrefLessMaxPref <- All$Pref - MaxPref
  linfitLessMaxPref <- lme(UnCorrected ~ PrefLessMaxPref, random=~1|Batch, na.action=na.exclude, data=All)
  seMaxInt <- summary(linfitLessMaxPref)$tTable[1,2]
  
  # Correct for variation between batches by fitting Batch to response and adding back
  # overall mean to residuals
  BatchFitted <- lme(UnCorrected ~ 1, random = ~1|Batch, na.action=na.exclude, data = All)
  Corrected <- residuals(BatchFitted) + mean(All$UnCorrected, na.rm=TRUE)
  
  # Append corrected data to working dataframe
  All <- cbind(All,Corrected)
  
  # Create separate dataframes for trained and untrained chicks
  Tr <- subset(All, TrUntr == "Trained")
  Untr <- subset(All, TrUntr == "Untrained")
  
  # Calculate total variance of trained chicks after removing effect of Batch
  full <- lme(UnCorrected ~ 1, random=~1|Batch, na.action=na.exclude, data=Tr)
  varTr <- as.numeric(VarCorr(full)[2,1])       # Variance
  varTrdf <- summary(full)$tTable[1,3]          # Degrees of freedomn for variance
  
  # Items of text for display of correlation coefficient
  s1 <- "r = "
  s2 <- as.character(round(r, digits = 3))
  s3 <- ', p = '
  s4 <- ifelse(tProb < 0.001, 
               format(tProb, scientific = TRUE, digits = 3),
               as.character(round(tProb, digits = 3)))
  rtext <- paste0(s1,s2,s3,s4)
  
  #Plot regression with preference score after removing effect of Batch
  par(mfrow=c(2,2))                           # Two rows, two plots    
  par(fig=c(0.2,1,0,1))                       # Specify frame for regression plot
  par(cex.lab=1.5,cex.axis=1.5,cex.main=1.5)  # Magnifications
  
  plot(Tr$Pref, Tr$Corrected, main=mainheading, ylim=c(ylow,yhigh), xlim=c(xlow,xhigh), frame.plot=FALSE, # Plot response against preference score
       xlab="Preference score", ylab="", xaxt = "n", yaxt="n", pch=16, cex=1.8,
       mgp=c(2,1,-0.71), cex.lab=1.8, cex.main=2)
  axis(side=1,lwd=2*lwd_scale, pos=ylow)
  clip(lowxforint, highxforint, ylow, yhigh)
  abline(a=Int, b=slope, xpd=F, lwd=2*lwd_scale)        # Fit regression line 
  arrows(50,ylow,50,Int50,length=0,angle=0,code=2,lty=2,lwd=2*lwd_scale)                 # Dashed line to y = In50
  arrows(MaxPref,ylow,MaxPref,IntMaxPref,length=0,angle=0,code=2,lty=2,lwd=2*lwd_scale)  # Dashed line to y = IntMaxPref
  Position <- ylow + (Int50-ylow)/2
  text(70, Position, rtext, cex=2)
    
  # Find residuals of regression plot and run correlations between preference score and
  # (i) absolute residuals (ii) squared residuals
  model1 <- lm(Corrected ~ Pref, data=Tr, na.action=na.exclude)
  resabs <- abs(model1$residuals)   # Absolute values of residuals
  ressq <- (model1$residuals)^2     # Squared values of residuals
  prefscores <- model1$model[,2]    # Preference scores
  rresabs <- cor(resabs,prefscores) # Correlation coefficient, absolute values
  rressq <- cor(ressq,prefscores)   # Correlation coefficient, squared values
  
  # Find mean, SEM and degrees of freedom for untrained chicks after removing effect
  # of Batch
  nUntr <- sum(!is.na(Untr$Corrected))              # Number of observations
  DFavg <- nUntr-1                                  # Degrees of freedom
  avg <- mean(Untr$Corrected, na.rm=TRUE)           # Mean
  sem <- sd(Untr$Corrected, na.rm=TRUE)/sqrt(nUntr) # SEM
  varUntr <- var(Untr$Corrected, na.rm=TRUE)        # Variance
  xtemp <- rep(10,nUntr)
   
  # Plot mean and SEM for untrained chicks
  par(fig=c(0,1,0,1),new=T)                         # Frame for plot
  
  plot(13, avg, frame.plot=FALSE, ylim=c(ylow,yhigh), xlim=c(0,100), xaxt="n",
       xlab="", ylab="", cex=2, lwd=2*lwd_scale)
  points(x=xtemp, y=Untr$Corrected, cex=1.5, lwd=2*lwd_scale)
  arrows(13, avg-sem, 13, avg+sem, length=0.075, angle=90, code=3, lwd=2*lwd_scale) # Plot +/- SEM
  axis(side=2,lwd=2*lwd_scale)
  title(ylab=ylabel, line=2.7, cex.lab=1.8)
  RelPosn = avg + sem + 0.05
  text(10,yhigh,"Untrained",cex=2)
  
  # Adjust length of horixontal dashed lines for y intercepts at MaxPref and Pref = 50
  # (depends on MinPref) and draw y intercepts
  if(MinPref <= 40) {xint <- 43}
  if(MinPref >= 50) {xint <- 24}
  if(MinPref > 40 & MinPref < 50) {xint <- 35}
  arrows(ylow,Int50,xint,Int50,length=0,angle=0,code=2,lty=2,lwd=2*lwd_scale)
  arrows(ylow,IntMaxPref,MaxPref,IntMaxPref,length=0,angle=0,code=2,lty=2,lwd=2*lwd_scale)
  
  # Plot SEs of intercepts
  col1 <- rgb(0.5,0.5,0.5,0.5)                # 50% transparency
  Batch <- Tr$Batch
  rect(-5,Int50-seInt50,0,Int50+seInt50,col=col1,lwd=lwd_scale)
  rect(-5,IntMaxPref-seMaxInt,0,IntMaxPref+seMaxInt,col=col1,lwd=lwd_scale)
  
  # t-tests comparing intercepts with untrained mean, degrees of freedom adjusted
  # for different variances
  DiffIMP_Mean <- IntMaxPref - avg        # Intercept at max. pref. vs untrained mean
  SEDIMP_Mean <- sqrt(sem^2 + seMaxInt^2)
  tIMP_Mean <- DiffIMP_Mean/SEDIMP_Mean
  DFIMP_Mean <- (sem^2 + seMaxInt^2)^2/((sem^4/DFavg) + (seMaxInt^4/resdf))
  ptIMP_Mean <- 2*(1 - pt(abs(tIMP_Mean),DFIMP_Mean))
  
  DiffInt50_Mean <- Int50 - avg           # Intercept at pref. 50 vs untrained mean
  SEInt50_Mean <- sqrt(sem^2 + seInt50^2)
  tInt50_Mean <- DiffInt50_Mean/SEInt50_Mean
  DFInt50_Mean <- (sem^2 + SEInt50_Mean^2)^2/((sem^4/DFavg) + (SEInt50_Mean^4/resdf))
  ptInt50_Mean <- 2*(1 - pt(abs(tInt50_Mean),DFInt50_Mean))
  
  # Calculate (residual regression variance)/(variance for untrained chicks)
  # and test whether residual regression variance is significantly lower.
  FRegr_Untr <- varFitres/varUntr         # (residual regression variance)/(variance for untrained chicks)
  FProb1 <- pf(FRegr_Untr,resdf,DFavg)    # P, 1-tailed, prediction F < 1
  FTr_Untr <- varTr/varUntr               # (total variance trained)/(variance untrained)
  FProb2 <- 1-pf(FTr_Untr,varTrdf,DFavg)  # P, 1-tailed, prediction F > 1
  
  # Output statistics
  Value <- c(avg,sem,DFavg,Int50,seInt50,IntMaxPref,seMaxInt,resdf,r,t,tProb,
             DiffIMP_Mean,SEDIMP_Mean,tIMP_Mean,DFIMP_Mean,ptIMP_Mean,
             DiffInt50_Mean,SEInt50_Mean,tInt50_Mean,DFInt50_Mean,ptInt50_Mean,
             varFitres,varUntr,FRegr_Untr,resdf,DFavg,FProb1,varTr,varTrdf,FTr_Untr,
             FProb2,rresabs,rressq)
  
  Variable <- c("avg","sem","DFavg","Int50","seInt50","IntMaxPref","seMaxInt",
                "resdf","r","t","tProb","DiffIMP_Mean","SEDIMP_Mean","tIMP_Mean",
                "DFIMP_Mean","ptIMP_Mean","DiffInt50_Mean","SEInt50_Mean","tInt50_Mean",
                "DFInt50_Mean","ptInt50_Mean","varFitres","varUntr","FRegr_Untr",
                "resdf","DFavg","FProb1","varTr","varTrdf","FTr_Untr","FProb2",
                "rresabs","rressq")
  
  Notes <- c("Mean untrained chicks",
             "SEM untrained chicks",
             "DF untrained chicks",
             "Intercept at preference score 50",
             "SE intercept at preference score 50",
             "Intercept maximum preference",
             "SE intercept maximum preference",
             "Residual DF regression",
             "Correlation coefficient",
             "t correlation",
             "P correlation",
             "Difference between untrained mean and intercept at maximum preference",
             "SE difference between untrained mean and intercept at maximum preference",
             "t difference between untrained mean and intercept at maximum preference",
             "DF difference between untrained mean and intercept at maximum preference",
             "P difference between untrained mean and intercept at maximum preference",
             "Difference between mean and intercept at preference score 50",
             "SE difference between mean and intercept at preference score 50",
             "t diff between untrained mean and intercept at preference score 50",
             "DF diff between untrained mean and intercept at preference score 50",
             "P diffbetween untrained mean and intercept at preference score 50",
             "Residual variance regression",
             "Variance untrained chicks",
             "(Residual variance regression)/(Variance untrained chicks)",
             "DF residual variance regression",
             "DF variance untrained chicks",
             "P (Residual variance regression)/(Variance untrained chicks)",
             "Total variance trained chicks",
             "DF total variance trained chicks",
             "(Total variance trained chicks)/(Variance untrained chicks)",
             "P (Total variance trained chicks)/(Variance untrained chicks)",
             "Correlation absolute residuals regression vs preference score",
             "Correlation squared residuals regression vs preference score")
  
  dfout <- data.frame(Variable,Value,Notes)
  
  return(dfout)
  
}
