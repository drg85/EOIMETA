# The following code is taken verbatim from iatal et al method,
# Source: https://github.com/iatal/fragility_ma/tree/master/app/functions_fragility
# License: MIT



#Revman Meta-analysis
rev_ma <- function(data,method,random,measure){

  #data should have as columns: EVENTS_1, EVENTS_2, TOTAL_1, TOTAL_2
  #method can be: "Inverse", "MH" or "PETO"
  #random can be: "YES" or "NO"
  #measure can be: "RR","RD" or "OR"

  dbin <- data
  if(nrow(dbin)==0) return(NULL)
  else{
    meta1 <- metabin(event.e=EVENTS_1, n.e=TOTAL_1, event.c=EVENTS_2, n.c=TOTAL_2,
                     data = dbin,
                     method = method,
                     sm = measure,
                     RR.Cochrane = TRUE,control = list(maxiter = 1000, tol = 1e-10, stepadj = 0.5))

    #For fragility, we calculate the exponnential of the CI, even if measure = RD
    #In order to have allways a comparison of CI to 1
    res <- ifelse(random=="NO",list(exp(c(meta1$lower.fixed,meta1$upper.fixed))),
                  list(exp(c(meta1$lower.random,meta1$upper.random))))

    return(unlist(res))
  }

}


#Fragility Function for statistically significant meta-analyses
############################################################################
frag_ma <- function(data,method,random,measure,with_progress = FALSE, progress = NULL){

  #data should have as columns: EVENTS_1, EVENTS_2, TOTAL_1, TOTAL_2
  #method can be: "Inverse", "MH" or "PETO"
  #random can be: "YES" or "NO"
  #measure can be: "RR","RD" or "OR"

  frag <- 0
  ddoi <- data
  nb_test <- 0

  init <- rev_ma(ddoi,method,random,measure)

  #Significant result < 1
  #We need to either add events to first arm, or supress events from second arm

  if(init[2] < 1){

    ci_up <- init[2]

    while(ci_up<1){
      #We test modification of upper bound of confidence interval when
      #adding events to first arm of each study

      Ladd <- sapply(1:nrow(ddoi),function(i){
        dch <- ddoi
        dch$EVENTS_1[i] <- dch$EVENTS_1[i] + 1
        if(dch$EVENTS_1[i]>dch$TOTAL_1[i]) return(0)
        else return(rev_ma(dch,method,random,measure)[2])

      })

      #We test modification of upper bound of confidence interval when
      #suppressing events to second arm of each study
      Lsupp <- sapply(1:nrow(ddoi),function(i){
        dch <- ddoi
        dch$EVENTS_2[i] <- dch$EVENTS_2[i] - 1
        if(dch$EVENTS_2[i]<0) return(0)
        else return(rev_ma(dch,method,random,measure)[2])
      })

      #We keep only the addition or suppression of event modifying the most
      #the upper bound of the confidence interval
      madd <- which.max(Ladd)
      msupp <- which.max(Lsupp)
      if(length(msupp)==0){
        ddoi$EVENTS_1[madd] <- ddoi$EVENTS_1[madd] + 1
        ci_up <- Ladd[madd]
      }
      else{
        if(Ladd[madd] > Lsupp[msupp]) {
          ddoi$EVENTS_1[madd] <- ddoi$EVENTS_1[madd] + 1
          ci_up <- Ladd[madd]
        }
        else {
          ddoi$EVENTS_2[msupp] <- ddoi$EVENTS_2[msupp] - 1
          ci_up <- Lsupp[msupp]
        }
      }
      #increment of fragility
      frag <- frag + 1

      nb_test <- nb_test + 2*nrow(ddoi)

      #Progress function for app
      if(with_progress) {
        progress$inc(2*nrow(ddoi),
                     detail = ifelse(frag < 50,
                                     paste(c("Calculating ",nb_test," new pooled treatment effects"),
                                           collapse = ""),
                                     paste(c("Wow! Fragility Index superior to 50! Still calculating ",
                                             nb_test," new pooled treatment effects"),
                                           collapse = ""))
        )
      }

    }
  }

  #Significant result > 1
  #We need to either add events to second arm, or supress events from first arm

  if(init[1] > 1){

    ci_low <- init[1]

    while(ci_low>1){
      #We test modification of upper bound of confidence interval when
      #adding events to first arm of each study
      Ladd <- sapply(1:nrow(ddoi),function(i){
        dch <- ddoi
        dch$EVENTS_2[i] <- dch$EVENTS_2[i] + 1
        if(dch$EVENTS_2[i]>dch$TOTAL_2[i]) return(Inf)
        else return(rev_ma(dch,method,random,measure)[1])
      })
      #We test modification of upper bound of confidence interval when
      #suppressing events to second arm of each study
      Lsupp <- sapply(1:nrow(ddoi),function(i){
        dch <- ddoi
        dch$EVENTS_1[i] <- dch$EVENTS_1[i] - 1
        if(dch$EVENTS_1[i]<0) return(Inf)
        else return(rev_ma(dch,method,random,measure)[1])
      })

      #We keep only the addition or suppression of event modifying the most
      #the upper bound of the confidence interval
      madd <- which.min(Ladd)
      msupp <- which.min(Lsupp)
      if(length(msupp)==0){
        ddoi$EVENTS_2[madd] <- ddoi$EVENTS_2[madd] + 1
        ci_low <- Ladd[madd]
      }
      else{
        if(Ladd[madd] < Lsupp[msupp]) {
          ddoi$EVENTS_2[madd] <- ddoi$EVENTS_2[madd] + 1
          ci_low <- Ladd[madd]
        }
        else {
          ddoi$EVENTS_1[msupp] <- ddoi$EVENTS_1[msupp] - 1
          ci_low <- Lsupp[msupp]
        }
      }

      #increment of fragility
      frag <- frag + 1
      nb_test <- nb_test + 2*nrow(ddoi)

      #Progress function for app
      if(with_progress) {
        progress$inc(2*nrow(ddoi),
                     detail = ifelse(frag < 50,
                                     paste(c("Calculating ",nb_test," new pooled treatment effects"),
                                           collapse = ""),
                                     paste(c("Wow! Fragility Index superior to 50! Still calculating ",
                                             nb_test," new pooled treatment effects"),
                                           collapse = ""))
        )
      }

    }
  }

  return(list(frag,ddoi))

}


#Fragility function for statistically non-significant results
frag_ma_ns <- function(data,method,random,measure,with_progress = FALSE, progress = NULL){

  #data should have as columns: EVENTS_1, EVENTS_2, TOTAL_1, TOTAL_2
  #method can be: "Inverse", "MH" or "PETO"
  #random can be: "YES" or "NO"
  #measure can be: "RR","RD" or "OR"

  ddoi <- data
  frag <- 0
  ddoiA <- ddoi
  ddoiB <- ddoi
  nb_test <- 0

  init <- rev_ma(ddoi,method,random,measure)

  #We test if it is possible to have a CI greater than 1 or a CI lower than 1 in the worst
  #case scenario
  dtestA0 <- ddoi
  dtestA0$EVENTS_1 <- dtestA0$TOTAL_1
  dtestA0$EVENTS_2 <- 0

  dtestB0 <- ddoi
  dtestB0$EVENTS_1 <- 0
  dtestB0$EVENTS_2 <- dtestB0$TOTAL_2

  #If in the worst case scenarios we still overlap 1, fragility is Inf
  worst_case <- c(rev_ma(dtestA0,method,random,measure)[1]<1,
                  rev_ma(dtestB0,method,random,measure)[2]>1)

  if(worst_case[1] & worst_case[2]) return(list(Inf,NULL))

  #If not, we count the minimal number of modifs to stop overlapping 1
  else{
    ci_low <- init[1]
    ci_up  <- init[2]

    while(ci_low<1 & ci_up>1){

      #We test modification of lower bound of confidence interval when
      #adding events to first arm of each study
      #Only if it would be possible to overlap 1 in WCS
      if(!worst_case[1]) {

        Ladd <- sapply(1:nrow(ddoiA),function(i){
          dch <- ddoiA
          dch$EVENTS_1[i] <- dch$EVENTS_1[i] + 1
          if(dch$EVENTS_1[i]>dch$TOTAL_1[i]) return(0)
          else return(rev_ma(dch,method,random,measure)[1])
        })
        #We test modification of lower bound of confidence interval when
        #suppressing events to second arm of each study
        Lsupp <- sapply(1:nrow(ddoiA),function(i){
          dch <- ddoiA
          dch$EVENTS_2[i] <- dch$EVENTS_2[i] - 1
          if(dch$EVENTS_2[i]<0) return(0)
          else return(rev_ma(dch,method,random,measure)[1])
        })
        #We keep only the addition or suppression of event modifying the most
        #the lower bound of the confidence interval to become the closest to 1
        madd <- which.max(Ladd)
        msupp <- which.max(Lsupp)
        if(length(msupp)==0){
          ddoiA$EVENTS_1[madd] <- ddoiA$EVENTS_1[madd] + 1
          ci_low <- Ladd[madd]
        }
        else{
          if(Ladd[madd] > Lsupp[msupp]) {
            ddoiA$EVENTS_1[madd] <- ddoiA$EVENTS_1[madd] + 1
            ci_low <- Ladd[madd]
          }
          else {
            ddoiA$EVENTS_2[msupp] <- ddoiA$EVENTS_2[msupp] - 1
            ci_low <- Lsupp[msupp]
          }
        }
      }

      #We test modification of lower bound of confidence interval when
      #adding events to first arm of each study
      #Only if it would be possible to overlap 1 in WCS
      if(!worst_case[2]) {

        #We test modification of upper bound of confidence interval when
        #adding events to second arm of each study
        Ladd <- sapply(1:nrow(ddoiB),function(i){
          dch <- ddoiB
          dch$EVENTS_2[i] <- dch$EVENTS_2[i] + 1
          if(dch$EVENTS_2[i]>dch$TOTAL_2[i]) return(Inf)
          else return(rev_ma(dch,method,random,measure)[2])

        })
        #We test modification of upper bound of confidence interval when
        #suppressing events to first arm of each study
        Lsupp <- sapply(1:nrow(ddoiB),function(i){
          dch <- ddoiB
          dch$EVENTS_1[i] <- dch$EVENTS_1[i] - 1
          if(dch$EVENTS_1[i]<0) return(Inf)
          else return(rev_ma(dch,method,random,measure)[2])
        })
        #We keep only the addition or suppression of event modifying the most
        #the upper bound of the confidence interval to become the closest to 1
        madd <- which.min(Ladd)
        msupp <- which.min(Lsupp)
        if(length(msupp)==0){
          ddoiB$EVENTS_2[madd] <- ddoiB$EVENTS_2[madd] + 1
          ci_up <- Ladd[madd]
        }
        else{
          if(Ladd[madd] < Lsupp[msupp]) {
            ddoiB$EVENTS_2[madd] <- ddoiB$EVENTS_2[madd] + 1
            ci_up <- Ladd[madd]
          }
          else {
            ddoiB$EVENTS_1[msupp] <- ddoiB$EVENTS_1[msupp] - 1
            ci_up <- Lsupp[msupp]
          }
        }
      }

      #increment of fragility
      frag <- frag + 1
      nb_test <- nb_test + 4*nrow(ddoi)

      #Progress function for app
      if(with_progress) {
        progress$inc(4*nrow(ddoi),
                     detail = ifelse(frag < 50,
                                     paste(c("Calculating ", nb_test," new pooled treatment effects"),
                                           collapse = ""),
                                     paste(c("Wow! Fragility Index superior to 50! Still calculating ",
                                             nb_test," new pooled treatment effects"),
                                           collapse = ""))
        )
      }


    }

  }

  if(ci_low>=1) DR <- ddoiA
  else DR <- ddoiB

  return(list(frag,DR))

}
