
#' @name eoifunc
#' @title Ellipse of Insignificance calculator for 2x2 trials
#' @description Ellipse of insignificance (EOI) analysis is a powerful refined fragility index to ascertaining robustness of results in dichotomous outcome trials, quantifying the degree of re-coding in both trial arms is required to lose (or gain) significance at a given level. It is powerful and flexible, and if the sensitivity and specificity data for a trial are also known, it can deduce whether reported results are robust. For full methods and citations, please see Grimes 2022 (Elife): "The ellipse of insignificance, a refined fragility index for ascertaining robustness of results in dichotomous outcome trials" - 10.7554/eLife.79573

#' @param a Number of positive cases in the experimental arm
#' @param b Number of negative cases in the experimental arm
#' @param c Number of positive cases in the control arm
#' @param d Number of negative cases in the control arm
#' @param v Optional alpha level to be tested (default p < 0.05)
#' @param vis Graphic summary on / off
#' @param ssvev An optional vector of experimental sensitivity, specificity and control sensitivity and specificity
#' @return Ellipse of insignificance statistics and refined fragility measures for any 2x2 binary outcome trial

#' @examples g <- eoifunc(700,300,500,500)
#' @examples g <- eoifunc(113,887,24,976, ssvec = c(0.75,0.9,0.68,0.99))

#' @importFrom rootSolve multiroot
#' @importFrom CVXR Maximize Minimize Problem Variable log_det norm2 psolve
#' @importFrom uniformly runif_in_ellipsoid runif_in_sphere runif_on_ellipsoid runif_on_sphere
#' @importFrom Carlson elliptic_E
#' @importFrom R6 R6Class
#' @importFrom fitConic fitConic
#' @importFrom rcdd makeV scdd
#' @importFrom graphics abline curve lines par polypath
#'
#'

# @section Attribution:
# Functions originally from the PlaneGeometry package (v1.6.0) by Michel van den Bergh
# have been adapted and included in this package. See plane_geometry_utils.R for details

#'
#' @export


eoifunc <- function(a,b,c,d, v = 0.05, vis = 1,ssvec = c(NA, NA, NA, NA)){
  alpha <- v  # Use v as the alpha level
  v <- qchisq(alpha, df = 1, lower.tail = FALSE)
  n <- a + b + c + d

    chstat <- (n * ((a * d - b * c)^2)) / ((a + b) * (c + d) * (a + c) * (b + d))
    pval <- pchisq(chstat, df=1, lower.tail=FALSE)


    tstat_test <- chstat - v
    vtest <- ifelse(tstat_test <= 0, 0, 1)

    # Calculations
      Ae <- (c + d) * ((c + d) * n + (a + b) * v)
      Be <- 2 * (a + b) * (c + d) * (n - v)
      Ce <- (a + b) * ((a + b) * n + (c + d) * v)
      De <- (c + d) * (2 * (b * c - a * d) * n + (a + b) * (b - a + d - c) * v)
      Ee <- (a + b) * (2 * (b * c - a * d) * n + (c + d) * (a - b + c - d) * v)
      Fe <- ((b * c - a * d)^2) * n - (a + b) * (a + c) * (b + d) * (c + d) * v

      # Find FECKUP vector points
      xp <- numeric()
      yp <- numeric()# Define the system of equations as functions
      # Define the equations

      # Define the equations
      # Define the equations
      equation1 <- function(xp, yp) {
        return ((2 * Ae * xp + Be * yp + De) * yp - xp * (Be * xp + 2 * Ce * yp + Ee))
      }

      equation2 <- function(xp, yp) {
        return (Ae * (xp^2) + Be * xp * yp + Ce * (yp^2) + De * xp + Ee * yp + Fe)
      }

      # Define a function that takes a vector of xp and yp and returns a vector of equations
      system_of_equations <- function(xy) {
        xp <- xy[1]
        yp <- xy[2]
        c(equation1(xp, yp), equation2(xp, yp))
      }

      # Number of solutions you want to find
      num_solutions <- 4

      # Initialize an empty list to store solutions
      solutions <- list()

      # Find all solutions
      attempts <- 0
      while (length(solutions) < num_solutions && attempts < 100) {
        # Generate random initial values
        initial_values <- c(runif(1, -n, n), runif(1, -n, n))

        # Use multiroot function to find roots
        result <- multiroot(system_of_equations, initial_values)

        # Calculate the distance between new solution and existing solutions
        distances <- sapply(solutions, function(sol) sum((sol - result$root)^2))

        # If the new solution is significantly different, store it
        if (all(distances > 1e-6)) {
          solutions <- c(solutions, list(result$root))
        }

        attempts <- attempts + 1
      }

      # Initialize a list to store the solution details
      solution_details <- list()

      # Display the solutions
      for (i in 1:length(solutions)) {
        #cat("Solution", i, ": xp =", solutions[[i]][1], ", yp =", solutions[[i]][2], "\n")
        sum_of_squares_i <- solutions[[i]][1]^2 + solutions[[i]][2]^2
        #cat("Sum of xp^2 + yp^2 for Solution", i, ":", sum_of_squares_i, "\n")

        # Store the solution details
        solution_details[[i]] <- list(xp = solutions[[i]][1], yp = solutions[[i]][2], sum_of_squares = sum_of_squares_i)

      }

      # Calculate the index with the smallest sum of squares
      dvmin <- which.min(sapply(solution_details, function(sol) sol$sum_of_squares))

      # Extract the values
      xpv <- as.numeric(solution_details[[dvmin]][1])
      ypv <- as.numeric(solution_details[[dvmin]][2])
      FKUP2 <- as.numeric(solution_details[[dvmin]][3])
      FKUP <- sqrt(FKUP2)


      # Find xi point
      A1 <- (c+d) * ((c+d) * n + (a+b) * v)
      B1 <- (c+d) * (2 * (b*c-a*d) * n + (a+b) * (b-a + d -c) * v)
      C1 <- ((b*c-a*d)^2) * n - (a+b) * (a+c) * (b+d) * (c+d) * v
      X1 <- (-B1 + sqrt(B1^2 - 4*A1*C1)) / (2*A1)
      X2 <- (-B1 - sqrt(B1^2 - 4*A1*C1)) / (2*A1)
      g <- c(X1, X2)
      gabs <- abs(g)
      xi <- (g[which.min(gabs)])


      # Find yi point
      A2 <- (a + b) * ((a + b) * n + (c + d) * v)
      B2 <- (a + b) * (2 * (b * c - a * d) * n + (a - b + c - d) * (c + d) * v)
      C2 <- ((b * c - a * d)^2) * n - (a + b) * (a + c) * (b + d) * (c + d) * v
      Y1 <- (-B2 + sqrt(B2^2 - 4*A2*C2)) / (2*A2)
      Y2 <- (-B2 - sqrt(B2^2 - 4*A2*C2)) / (2*A2)
      g2 <- c(Y1, Y2)
      gabs2 <- abs(g2)
      yi <- (g2[which.min(gabs2)])



      # Find all relevant terms
      dmin <- ceiling(abs(xpv) + abs(ypv))
      qvec <- c(1 - (a + b - abs(xi)) / (a + b),
                1 - (c + d - abs(yi)) / (c + d),
                1 - (n - abs(dmin)) / n)


      xm <- NA
      ym <- NA


      # Check if all values in ssvec are real numbers between 0 and 1
      if (!any(is.na(ssvec)) && all(is.numeric(ssvec)) && all(ssvec >= 0 & ssvec <= 1)) {


        sne <- ssvec[1]
        spe <- ssvec[2]
        snc <- ssvec[3]
        spc <- ssvec[4]

        xm <- abs((b*(1-spe)-a*(1-sne)) / (sne + spe - 1))
        ym <- abs((c*(1-snc)-d*(1-spc)) / (snc + spc - 1))
      }




      output_list <- list(
        totalsamp = n,
        hazardrat = a*(c+d)/(c*(a+b)),
        alpha = alpha,
        p_value = pval,
        signif = vtest,
        min_ellipse_x = xpv,
        min_ellipse_y = ypv,
        xi_min = xi,
        yi_min = yi,
        exp_frac = qvec[1],
        con_frac = qvec[2],
        total_frac = qvec[3],
        FECKUP = FKUP,
        dmin = dmin,
        ssvec = ssvec,
        xm_ss = xm,
        ym_ss = ym

      )

      if (vis==1) {

        #library(ggplot2)



        csize <- 1
        xtight <- sort(c(0-0.05*xi,1.05*xi))
        ytight <- sort(c(0-0.05*yi,1.05*yi))

        ell <- EllipseFromEquation(A = Ae, B = Be, C = Ce, D = De, E = Ee, F = Fe)
        origincirc <- Circle$new(c(0,0), csize)
        xicirc <- Circle$new(c(xi,0), csize)
        yicirc <- Circle$new(c(0,yi), csize)
        fkcirc <- Circle$new(c(xpv,ypv), csize)
        box <- ell$boundingbox()
        plot(NULL, xlim = xtight, ylim = ytight, xlab = "Experimental group re-coding (x)", ylab = "Control group re-coding (y)", main = "Ellipse of insignificance (region of interest)")
        draw(ell, col = "lightblue", border = "purple", lwd = 3)
        draw(origincirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
        draw(xicirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
        #draw(yicirc, col = "white", border = "black", lwd = 1)
        draw(yicirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
        draw(fkcirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
        segments(x0 = 0, y0 = 0, x1 = xi, y1 = 0, col = "darkgreen", lty = "dashed")
        segments(x0 = 0, y0 = 0, x1 = 0, y1 = yi, col = "darkgreen", lty = "longdash")
        segments(x0 = 0, y0 = 0, x1 = xpv, y1 = ypv, col = "red", lty = "dotted", lwd = 3)
        text(c(xpv/2),c(ypv/2),labels=c("FECKUP"))
        text(c(0),c(yi/2),labels=c("|yi|"))
        text(c(xi/2),c(0),labels=c("|xi|"))

      }






          cat("\n**Total Sample:**", output_list$totalsamp, "\n")
          cat("**p-value:**", output_list$p_value, "\n")
          cat("**Hazard Ratio (Exp/Con):**", output_list$hazardrat, "\n")


          cat("\n**FECK UP Vector co-ordinates:**\n")
          cat("  X-point:", output_list$min_ellipse_x, "\n")
          cat("  Y-point:", output_list$min_ellipse_y, "\n")
          cat("  FECKUP vector length", output_list$FECKUP, "\n")
          cat("  Resolved FECKUP length (dmin) ", output_list$dmin, "subjects" , "\n")


        cat("\n**Single arm re-coding summary:**\n")
        cat("  Minimal Exp Arm only recoding (xi) :", output_list$xi, "\n")
        cat("  Minimal Con Arm only recoding (yi) :", output_list$yi, "\n")



        # exp_frac, con_frac, total_frac: Expected, control, and total fractions
        cat("\n**Fractions:**\n")
        cat("  Experimental arm re-coding fraction (%):", 100*output_list$exp_frac, "\n")
        cat("  Control arm re-coding fraction (%):", 100*output_list$con_frac, "\n")
        cat("  Total sample re-coding fraction (%):", 100*output_list$total_frac, "\n")


        if (!any(is.na(ssvec)) && all(is.numeric(ssvec)) && all(ssvec >= 0 & ssvec <= 1)) {
          cat("\n**Sensitivity and Specificity calculations:**\n")
          cat(" Experimental Sensitiity and Specificty (%):", 100*c(sne,spe), "\n")
          cat(" Control Sensitiity and Specificty (%):", 100*c(snc,spc), "\n")
          cat(" Minimal number miscoded in Exp arm (xm) :", xm , "\n")
          cat(" Minimal number miscoded in Con arm (ym) :", ym , "\n")

          if (abs(xm) > abs(xi)) {
            cat("\n**WARNING: Minimal miscoded in Experimental arm exceeds xi - RESULT SUSPECT:**\n")

          }

          if (abs(ym) > abs(yi)) {
            cat("\n**WARNING: Minimal miscoded in Control arm exceeds yi - RESULT SUSPECT:**\n")

          }





          }




    return(invisible(output_list))





}


