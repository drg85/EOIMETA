
#' @name eoimeta
#' @title Ellipse of Insignificance calculator for Meta Analyses of 2x2 trials
#' @description Extension of Ellipse of insignificance (EOI) and Region of Attainable Redaction (ROAR) analysis meta for 2x2 trials to meta-analysis.

#' @param av List of positive cases in the experimental arm of individual trials
#' @param bv List of non-events in the experimental arm of individual trials
#' @param cv List of positive cases in the control arm of individual trials
#' @param dv List of non-events in the control arm of individual trials
#' @param v Optional alpha level to be tested (default p < 0.05)
#' @param atal Run Atal et al study specific fragility algorithm (default = 0, 1 for on)
#' @param wiv Used generic weighed inverse variance correction for included studies (default = 1)

#' @return Ellipse of insignificance statistics and refined fragility measures for any 2x2 binary outcome trial




#' @importFrom rootSolve multiroot

#' @importFrom meta metabin
#'
#' @importFrom R6 R6Class

#' @importFrom Carlson elliptic_E
#'
t
#' @importFrom CVXR Maximize Minimize Problem Variable log_det norm2 psolve
#' @importFrom uniformly runif_in_ellipsoid runif_in_sphere runif_on_ellipsoid runif_on_sphere

#' @importFrom fitConic fitConic
#' @importFrom rcdd makeV scdd
#' @importFrom graphics abline curve lines par polypath
#'

#' @export


eoimetafix <- function(av,bv,cv,dv, v = 0.05,atal = 0, wiv = 1){
  atalv <- atal
  alpha <- v  # Use v as the alpha level
  v <- qchisq(alpha, df = 1, lower.tail = FALSE)



  #EOI META ANALYSIS CMH
 # eoimeta <- metabin(av,av + bv,cv,cv + dv,
  #                   common = T,random = T,
   #                  method = 'MH',sm = "RR") # Mantel Haenszel weighting

 if(wiv ==1){
   #Continuity correction
  av[av == 0] <- 0.5
  bv[bv == 0] <- 0.5
  cv[cv == 0] <- 0.5
  dv[dv == 0] <- 0.5
 }

  eoimeta <- metabin(av, av + bv, cv, cv + dv,
                     common = TRUE, random = TRUE,
                     method = "MH", sm = "RR",
                     incr = 0.5, allstudies = TRUE,
                     control = list(maxiter = 1000, tol = 1e-10, stepadj = 0.5))


  # Extract Fixed-Effect Relative Risk
  RR_fixed <- exp(eoimeta$TE.fixed)
  #RR_fixed_CI <- exp(c(eoimeta$lower.fixed, eoimeta$upper.fixed))

  # Extract Random-Effects Relative Risk
  RR_random <- exp(eoimeta$TE.random)
  #RR_random_CI <- exp(c(eoimeta$lower.random, eoimeta$upper.random))

  pvalcmh <- max(c(eoimeta$pval.common, eoimeta$pval.random))

  rgroup <- c(RR_fixed,RR_random)


  ap <- sum(av)
  bp <- sum(bv)
  cp <- sum(cv)
  dp <- sum(dv)


  if(wiv==1){

    logrr <- log((av*(cv +dv)/(cv*(av + bv))))
    selogrr <- sqrt(1/av + 1/(av+bv) + 1/cv + 1/(cv+dv))
    wi <- 1/(selogrr^2)
    top <- sum(wi*logrr)
    bot <- sum(wi)
    rrpg <- exp(top/bot)
    #ap <- round((bp*cp*rrpg)/(dp + cp*(1-rrpg)))
    xrr <- (cp*(ap+bp)*rrpg)/(cp+dp) - ap
    ap <- round(ap + xrr)
    bp <- round(bp - xrr)

  }


  RR_crude <- ap/(ap + bp) / (cp / (cp + dp))


  nps <- ap + bp + cp + dp

  chstat <- (nps * ((ap * dp - bp * cp)^2)) / ((ap + bp) * (cp + dp) * (ap + cp) * (bp + dp))

  pvalcrude <- pchisq(chstat, df=1, lower.tail=FALSE)

  #RRcheck
  dev <- max(abs(1 - rgroup/RR_crude))

  if (dev < 0.1) {
    het <- 0
  } else if ( dev >= 0.1 & dev < 0.2) {
    het <- 1

  } else {
    het <- 2
  }



  # creating a data frame
  eoiroar <- data.frame(
    study_number = c(1:length(av))
  )

  eoiroar$exp_pos_a <- NA  # Create an empty column to store results
  eoiroar$exp_neg_b <- NA  # Create an empty column to store results
  eoiroar$con_pos_c <- NA  # Create an empty column to store results
  eoiroar$con_neg_d <- NA  # Create an empty column to store results
  eoiroar$RRpoint <- NA  # Create an empty column to store results
  eoiroar$dmin <- NA
  eoiroar$rmin <- NA  # Create an empty column to store results
  eoiroar$xpv <- NA  # Create an empty column to store results
  eoiroar$ypv <- NA  # Create an empty column to store results
  eoiroar$pval <- NA  # Create an empty column to store results

  for (i in 1:length(av)) {
    eoiroar$exp_pos_a[i] <- av[i]
    eoiroar$exp_neg_b[i] <- bv[i]
    eoiroar$con_pos_c[i] <- cv[i]
    eoiroar$con_neg_d[i] <- dv[i]

    eoiroar$RRpoint[i] <- (av[i]/(bv[i] + av[i]))/(cv[i]/(cv[i] + dv[i]))

    h <- suppressMessages(eoifunc(av[i],bv[i],cv[i],dv[i], v = alpha, vis = 0))
    eoiroar$dmin[i] <- h$dmin
    eoiroar$xpv[i] <- h$min_ellipse_x
    eoiroar$ypv[i] <- h$min_ellipse_y
    eoiroar$pval[i] <- h$p_value

    h2 <- suppressMessages(roarfunc(av[i],bv[i],cv[i],dv[i], v = alpha, vis = 0))
    #eoiroar$rmin[i] <- h2$rmin


    if (!is.null(h2$rmin) && length(h2$rmin) > 0) {
      eoiroar$rmin[i] <- h2$rmin
    } else {
      eoiroar$rmin[i] <- NA  # or some other default value
    }




  }

  h <- suppressMessages(eoifunc(ap,bp,cp,dp, v = alpha, vis = 0))
  dminc <- h$dmin
  exprecode <- h$exp_frac
  conrecode <- h$con_frac
  totalrecode <- h$total_frac

  h2 <- suppressMessages(roarfunc(ap,bp,cp,dp, v = alpha, vis = 0))
  rminc <- h2$rmin
  expredact <- h2$exp_frac
  conredact <- h2$con_frac
  totalredact <- h2$total_frac


  cat("\014")





  # Compute the sum of each column
  final_row <- data.frame(
    study_number = NA,
    exp_pos_a = sum(eoiroar$exp_pos_a, na.rm = TRUE),
    exp_neg_b = sum(eoiroar$exp_neg_b, na.rm = TRUE),
    con_pos_c = sum(eoiroar$con_pos_c, na.rm = TRUE),
    con_neg_d = sum(eoiroar$con_neg_d, na.rm = TRUE),
    RRpoint = RR_crude ,  # RRpoint doesn't sum meaningfully, so set to NA
    dmin = dminc,  # Assuming dmin should also be NA
    rmin = rminc,   # Assuming rmin should also be NA
    xpv = NA,
    ypv = NA,
    pval = eoimeta$pval.random
  )

  # Append the final row to the dataframe
  eoiroar <- rbind(eoiroar, final_row)

  np <- ap + bp + cp + dp




  EVENTS_1 <- av
  TOTAL_1 <- av + bv
  EVENTS_2 <- cv
  TOTAL_2 <- cv + dv
  frameddata <- data.frame(EVENTS_1,TOTAL_1,EVENTS_2,TOTAL_2)
  #revman check
  sigcheck <- rev_ma(frameddata,"MH","YES","RR")
  atalvec <- NA
  atal_frag <- NA


  if(atalv == 1){
  #this part runs iatial's method for comparison!

 #  EVENTS_1 = av
 #  TOTAL_1 = av + bv
 #  EVENTS_2 = cv
 #  TOTAL_2 = cv + dv
 #  frameddata <- data.frame(EVENTS_1,TOTAL_1,EVENTS_2,TOTAL_2)
 #  #revman check
 # sigcheck <- rev_ma(frameddata,"MH","YES","RR")

 #check if sigcheck crosses unity

 if (min(sigcheck) < 1 && max(sigcheck) > 1) {
   gg <- suppressWarnings(frag_ma_ns(frameddata,"MH","YES","RR"))
 } else {
   gg <- suppressWarnings(frag_ma(frameddata,"MH","YES","RR"))
 }

atal_frag <- gg[[1]]
atalvec <- gg[[2]]
  }


  sigcheckloop <- vector("list", length(av))

#this part simulates individual flipping of studies
  for (i in 1:length(av)) {

    avt = av
    bvt = bv
    cvt = cv
    dvt = dv

    avt[i] <- round(av[i] - eoiroar$xpv[i])
    bvt[i] <- round(bv[i] + eoiroar$xpv[i])
    cvt[i] <- round(cv[i] + eoiroar$ypv[i])
    dvt[i] <- round(dv[i] - eoiroar$ypv[i])


    frameddata <- data.frame(
      EVENTS_1 = avt,
      TOTAL_1 = avt + bvt,
      EVENTS_2 = cvt,
      TOTAL_2 = cvt + dvt
    )
    sigcheckloop[[i]] <- rev_ma(frameddata,"MH","YES","RR")
    rm(avt)
    rm(bvt)
    rm(cvt)
    rm(dvt)

  }




  output_list <- list(
    totalsamp = np,
    alpha = alpha,
    rrcmh = rgroup,
    rcrude = RR_crude,
    RR_random = RR_random,
    dev = dev,
    het = het,
    pvalcmh = pvalcmh,
    pvalcrude = pvalcrude,
    ap = ap,
    bp = bp,
    cp = cp,
    dp = dp,
    np = ap + bp + cp + dp,
    ci95 = sigcheck,
    exprecode = exprecode,
    conrecode = conrecode,
    totalrecode = totalrecode,
    expredact = expredact,
    conredact = conredact,
    totalredact = totalredact,
    dmin = dminc,
    rmin = rminc,
    feedvec = c(ap,bp,cp,dp),
    atal_frag = atal_frag
  )

  cat("\n**Meta Summary:**", output_list$np, "\n")
  cat("**p-value (cmh):**", output_list$pvalcmh, "\n")
  cat("**p-value (crude):**", output_list$pvalcrude, "\n")
  cat("**RR random:**", output_list$RR_random, "\n")
  cat("**CMH Hazard Ratio (Exp/Con):**", output_list$rrcmh, "\n")
  cat("**Crude Hazard Ratio (Exp/Con):**", output_list$rcrude, "\n")
  cat("**95% Confidence bounds (Exp/Con):**", output_list$ci95, "\n")
  cat("**Confounding (%):**", 100*output_list$dev, "\n")

  if (het == 0){

    cat("\n**Low degree of confounding, grouping is possible:**\n")

  }

  if (het == 1){

    cat("\n**Moderate degree of confounding (10-20%), interpret grouping with caution :**\n")

  }

  if (het == 2){

    cat("\n**WARNING: High degree of confounding (>20%), grouping is unreliable :**\n")

  }


  cat("\n**Crude pooled EOI / ROAR figures:**\n")
  cat("  Total events (Experimental arm) :", output_list$ap, "\n")
  cat("  Total non-events (Experimental arm) :", output_list$bp, "\n")
  cat("  Total events (Control arm) :", output_list$cp, "\n")
  cat("  Total non-events (Control arm) :", output_list$dp, "\n")



  # exp_frac, con_frac, total_frac: Expected, control, and total fractions
  cat("\n**EOI Fractions:**\n")
  cat("  Experimental arm re-coding fraction (%):", 100*output_list$exprecode, "\n")
  cat("  Control arm re-coding fraction (%):", 100*output_list$conrecode, "\n")
  cat("  Total sample re-coding fraction (%):", 100*output_list$totalrecode, "\n")

  # exp_frac, con_frac, total_frac: Expected, control, and total fractions
  cat("\n**ROAR Fractions:**\n")
  cat("  Experimental arm re-coding fraction (%):", 100*output_list$expredact, "\n")
  cat("  Control arm re-coding fraction (%):", 100*output_list$conredact, "\n")
  cat("  Total sample re-coding fraction (%):", 100*output_list$totalredact, "\n")

  # exp_frac, con_frac, total_frac: Expected, control, and total fractions
  cat("\n**Summary stats:**\n")
  cat("  Resolved FECKUP vector (EOI, dmin):", output_list$dmin, "\n")
  cat("  Resolved FOCK vector (ROAR, rmin):", output_list$rmin, "\n")


  # Atal method fragility index
  cat("\n**Atal's method fragility index:**\n")
  cat("  Atal fragility:", output_list$atal_frag , "\n")




  result <- list(
    eoiroar = eoiroar,     # The data frame with the results
    output_list = output_list,  # The list with summary statistics
    sigcheckloop = sigcheckloop, #The list of flips statistics
    eoimeta = eoimeta,
    atalvec = atalvec #Atal's vector
  )

  return(invisible(result))

# return(invisible(output_list))
  #return(eoiroar)

}


