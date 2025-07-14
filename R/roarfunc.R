
#' @name roarfunc
#' @title Region of Attainable redaction Analysis for 2x2 trials
#' @description Region of Attainable redaction (ROAR) Analysis is an extention of Ellipse of Insignificance analysis for binary 2x2 trials that calculates the degree of redaction required to change a result from insignificant to significant or vice versa. It allows estimation of the degree of 'hidden' results required in control and experimental arms to change from null to non-null. For full methods and citations, please see Grimes 2024 (Elife): "Region of Attainable Redaction, an extension of Ellipse of Insignificance analysis for gauging impacts of data redaction in dichotomous outcome trials " - 10.7554/eLife.93050

#' @param a Number of positive cases in the experimental arm
#' @param b Number of negative cases in the experimental arm
#' @param c Number of positive cases in the control arm
#' @param d Number of negative cases in the control arm
#' @param v Optional alpha level to be tested (default p < 0.05)
#' @param vis Graphic summary on / off
#' @return ROAR statistics and hypothetical redaction levels for any 2x2 binary outcome trial

#' @examples h <- roarfunc(700,300,500,500)
#' @examples h <- roarfunc(250,550,400,500, v = 0.01)

#' @importFrom rootSolve multiroot
#' @importFrom ggplot2 ggplot labs theme annotate aes geom_tile guides scale_fill_manual theme_classic element_text
#'
#' @export

roarfunc <- function(a,b,c,d, v = 0.05, vis=1){  #matched
  alpha <- v  # Use v as the alpha level
  v <- qchisq(alpha, df = 1, lower.tail = FALSE)
  n <- a + b + c + d

  chstat <- (n * ((a * d - b * c)^2)) / ((a + b) * (c + d) * (a + c) * (b + d))
  pval <- pchisq(chstat, df=1, lower.tail=FALSE)


  tstat_test <- chstat - v
  vtest <- ifelse(tstat_test <= 0, 0, 1)


  if ((a + b == 0) || (c + d == 0) || a == 0 || c == 0) {
    return(list(
      pval = NA,
      chstat = NA,
      ratiocheck = NA,
      RRbound_upp = NA,
      RRbound_low = NA,
      xe = NA,
      ye = NA,
      FOCK = NA,
      xc = NA,
      yc = NA
    ))
  }









suppressWarnings({

  ratiocheck <- (a / (a + b)) / (c / (c + d))
  RRbound_upp <- exp(log(ratiocheck) + 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)))
  RRbound_low <- exp(log(ratiocheck) - 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)))

  if (ratiocheck > 1) {
    coeffx3y2 <- 1
    coeffx3y <- 2 * c
    coeffx3 <- c^2
    coeffx2y3 <- 1
    coeffx2y2 <- n + 2 * (b + c) - v
    coeffx2y <- 4 * b * c + c^2 - 2 * a * d + 2 * c * n - (a + 2 * c + d) * v
    coeffx2 <- c * (2 * b * c - 2 * a * d + c * n) - (a + c) * (c + d) * v
    coeffxy3 <- 2 * b
    coeffxy2 <- b^2 + 4 * b * c - 2 * a * d + 2 * b * n - (a + 2 * b + d) * v
    coeffxy <- 2 * (b + c) * (b * c - a * d) + 4 * b * c * n - 2 * a * d * n - (a + 2 * b + d) * (a + 2 * c + d) * v
    coeffx <- (b * c - a * d) * (b * c - a * d + 2 * c * n) - (a + c) * (a + 2 * b + d) * (c + d) * v
    coeffy3 <- b^2
    coeffy2 <- b * (-2 * a * d + b * (2 * c + n)) - (a + b) * (b + d) * v
    coeffy <- (b * c - a * d) * (-a * d + b * (c + 2 * n)) - (a + b) * (b + d) * (a + 2 * c + d) * v
    coeffc <- (b * c - a * d)^2 * n - (a + b) * (a + c) * (b + d) * (c + d) * v

    gfull <- function(x, y) {
      return(coeffx3y2 * (x^3 * y^2) + coeffx3y * (x^3 * y) + coeffx3 * (x^3) +
               coeffx2y3 * (x^2 * y^3) + coeffx2y2 * (x^2 * y^2) + coeffx2y * (x^2 * y) +
               coeffx2 * (x^2) + coeffxy3 * (x * y^3) + coeffxy2 * (x * y^2) +
               coeffxy * (x * y) + coeffx * (x) + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }

    disfun <- function(x, y) {
      return(sqrt(x^2 + y^2))
    }


    dgdx <- function(x, y) {
      return(3*coeffx3y2 * (x^2 * y^2) + 3*coeffx3y * (x^2 * y) + 3*coeffx3 * (x^2) +
               2*coeffx2y3 * (x * y^3) + 2*coeffx2y2 * (x * y^2) + 2*coeffx2y * (x * y) +
               2*coeffx2 * (x) + coeffxy3 * (y^3) + coeffxy2 * ( y^2) +
               coeffxy * (y) + coeffx * (1))
    }

    dgdy <-  function(x, y) {
      return(2*coeffx3y2 * (x^3 * y) + coeffx3y * (x^3 * 1) + 0 +
               3*coeffx2y3 * (x^2 * y^2) + 2*coeffx2y2 * (x^2 * y) + coeffx2y * (x^2 * 1) +
               0 + 3*coeffxy3 * (x * y^2) + 2*coeffxy2 * (x * y) +
               coeffxy * (x * 1) + 0 + 3*coeffy3 * (y^2) + 2*coeffy2 * (y) +
               coeffy * (1) + 0)
    }





    dDdx <- function(x, y) {
      return(x/sqrt(x^2 + y^2))
    }


    dDdy <- function(x, y) {
      return(y/sqrt(x^2 + y^2))
    }

    eqsimp <- function(x, y) {
      return(y * dgdx(x, y) - x * dgdy(x, y))
    }


    # Define a function that takes a vector of xp and yp and returns a vector of equations
    system_of_equations <- function(xy) {
      x <- xy[1]
      y <- xy[2]
      c(gfull(x, y), eqsimp(x, y))
    }

    # Number of solutions you want to find
    num_solutions <- 25

    # Initialize an empty list to store solutions
    solutions <- list()

    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 500) {
      # Generate random initial values
      initial_values <- c(runif(1, -0.1*n, n), runif(1, -0.1*n, n))

      # Use multiroot function to find roots
      result <- multiroot(system_of_equations, initial_values)

      # Calculate if root is +
      rsign1 <- result$root[[1]]
      rsign2 <- result$root[[2]]

      # Calculate the distance between new solution and existing solutions
      distances <- sapply(solutions, function(sol) sum((sol - result$root)^2))

      # If the new solution is significantly different, store it
      if (all(distances > 1e-6) && rsign1 > 0 && rsign2 > 0) {
        solutions <- c(solutions, list(result$root))
      }

      attempts <- attempts + 1
    }

    # Initialize a list to store the solution details
    solution_details <- list()

    # Display the solutions
    for (i in 1:length(solutions)) {
      #cat("Solution", i, ": x =", solutions[[i]][1], ", y =", solutions[[i]][2], "\n")
      sum_of_squares_i <- solutions[[i]][1]^2 + solutions[[i]][2]^2
      #cat("Sum of x^2 + y^2 for Solution", i, ":", sum_of_squares_i, "\n")

      # Store the solution details
      solution_details[[i]] <- list(xp = solutions[[i]][1], yp = solutions[[i]][2], sum_of_squares = sum_of_squares_i)

    }

    # Calculate the index with the smallest sum of squares
    dvmin <- which.min(sapply(solution_details, function(sol) sol$sum_of_squares))

    # Extract the values
    xe <- as.numeric(solution_details[[dvmin]][1])
    ye <- as.numeric(solution_details[[dvmin]][2])
    FOCK <- sqrt(as.numeric(solution_details[[dvmin]][3]))

    #now we solve for xc and yc. first define...


    gyfull <- function(y) {
      return(0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }

    gxfull <- function(x) {
      return(0 + 0 + coeffx3 * (x^3) +
               0 + 0 + 0 +
               coeffx2 * (x^2) + 0 + 0 +
               0 + coeffx * (x) + 0 + 0 +
               0 + coeffc)
    }


    #search for 2 roots around xe,ye
    # Number of solutions you want to find
    num_solutions <- 3

    # Initialize an empty list to store solutions
    solutionsxc <- list()
    solutionsyc <- list()

    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 100) {

      # Use multiroot function to find roots
      resultxc <- multiroot(gxfull, xe)
      resultyc <- multiroot(gyfull, ye)

      # Calculate if root is +
      rsignx <- resultxc$root
      rsigny <- resultyc$root

      # Calculate the distance between new solution and existing solutions
      distancesxc <- sapply(solutionsxc, function(sol) sum((sol - resultxc$root)^2))
      distancesyc <- sapply(solutionsyc, function(sol) sum((sol - resultyc$root)^2))

      # If the new solution is significantly different, store it
      if (all(distancesxc > 1e-6) && rsignx > 0) {
        solutionsxc <- c(solutionsxc, list(resultxc$root))
      }

      # If the new solution is significantly different, store it
      if (all(distancesyc > 1e-6) && rsigny > 0) {
        solutionsyc <- c(solutionsyc, list(resultyc$root))
      }


      attempts <- attempts + 1
    }

    yc <- min(as.numeric(solutionsyc))
    xc <- min(as.numeric(solutionsxc))



  } else {
    coeffx3y2 <- 1
    coeffx3y <- 2*d
    coeffx3 <- d^2
    coeffx2y3 <- 1
    coeffx2y2 <- n + 2*(a + d) - v
    coeffx2y <- -2*b*c + d*(4*a + d + 2*n) - v*(b + c + 2*d)
    coeffx2 <- d*(-2 *b* c + d*(2 *a + n)) - (b + d)*(c + d)*v
    coeffxy3 <- 2*a
    coeffxy2 <- -2*b*c + a*(a + 4*d + 2*n) - (2*a + b + c)*v
    coeffxy <- 2*(a + d)*(-b*c + a*d) - 2*b*c*n + 4*a*d*n - (2*a + b + c)*(b + c + 2*d)*v
    coeffx <- (b*c - a*d)*(b*c - d*(a + 2*n)) - (2*a + b + c)*(b + d)*(c + d)*v
    coeffy3 <- a^2
    coeffy2 <- a*(-2*b*c + a*(2*d + n)) - (a + b)*(a + c)*v
    coeffy <- (b*c - a*d)*(b*c - a*(d + 2*n)) - (a + b)*(a + c)*(b + c + 2*d)*v
    coeffc <- ((b*c - a*d)^2)*n - (a + b)*(a + c)*(b + d)*(c + d)*v

    gfull <- function(x, y) {
      return(coeffx3y2 * (x^3 * y^2) + coeffx3y * (x^3 * y) + coeffx3 * (x^3) +
               coeffx2y3 * (x^2 * y^3) + coeffx2y2 * (x^2 * y^2) + coeffx2y * (x^2 * y) +
               coeffx2 * (x^2) + coeffxy3 * (x * y^3) + coeffxy2 * (x * y^2) +
               coeffxy * (x * y) + coeffx * (x) + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }

    disfun <- function(x, y) {
      return(sqrt(x^2 + y^2))
    }


    dgdx <- function(x, y) {
      return(3*coeffx3y2 * (x^2 * y^2) + 3*coeffx3y * (x^2 * y) + 3*coeffx3 * (x^2) +
               2*coeffx2y3 * (x * y^3) + 2*coeffx2y2 * (x * y^2) + 2*coeffx2y * (x * y) +
               2*coeffx2 * (x) + coeffxy3 * (y^3) + coeffxy2 * ( y^2) +
               coeffxy * (y) + coeffx * (1))
    }

    dgdy <-  function(x, y) {
      return(2*coeffx3y2 * (x^3 * y) + coeffx3y * (x^3 * 1) + 0 +
               3*coeffx2y3 * (x^2 * y^2) + 2*coeffx2y2 * (x^2 * y) + coeffx2y * (x^2 * 1) +
               0 + 3*coeffxy3 * (x * y^2) + 2*coeffxy2 * (x * y) +
               coeffxy * (x * 1) + 0 + 3*coeffy3 * (y^2) + 2*coeffy2 * (y) +
               coeffy * (1) + 0)
    }





    dDdx <- function(x, y) {
      return(x/sqrt(x^2 + y^2))
    }


    dDdy <- function(x, y) {
      return(y/sqrt(x^2 + y^2))
    }

    eqsimp <- function(x, y) {
      return(y * dgdx(x, y) - x * dgdy(x, y))
    }


    # Define a function that takes a vector of xp and yp and returns a vector of equations
    system_of_equations <- function(xy) {
      x <- xy[1]
      y <- xy[2]
      c(gfull(x, y), eqsimp(x, y))
    }

    # Number of solutions you want to find
    num_solutions <- 25

    # Initialize an empty list to store solutions
    solutions <- list()

    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 500) {
      # Generate random initial values
      initial_values <- c(runif(1, -0.1*n, n), runif(1, -0.1*n, n))

      # Use multiroot function to find roots
      result <- multiroot(system_of_equations, initial_values)

      # Calculate if root is +
      rsign1 <- result$root[[1]]
      rsign2 <- result$root[[2]]

      # Calculate the distance between new solution and existing solutions
      distances <- sapply(solutions, function(sol) sum((sol - result$root)^2))

      # If the new solution is significantly different, store it
      if (all(distances > 1e-6) && rsign1 > 0 && rsign2 > 0) {
        solutions <- c(solutions, list(result$root))
      }

      attempts <- attempts + 1
    }

    # Initialize a list to store the solution details
    solution_details <- list()

    # Display the solutions
    for (i in 1:length(solutions)) {
      #cat("Solution", i, ": x =", solutions[[i]][1], ", y =", solutions[[i]][2], "\n")
      sum_of_squares_i <- solutions[[i]][1]^2 + solutions[[i]][2]^2
      #cat("Sum of x^2 + y^2 for Solution", i, ":", sum_of_squares_i, "\n")

      # Store the solution details
      solution_details[[i]] <- list(xp = solutions[[i]][1], yp = solutions[[i]][2], sum_of_squares = sum_of_squares_i)

    }

    # Calculate the index with the smallest sum of squares
    dvmin <- which.min(sapply(solution_details, function(sol) sol$sum_of_squares))

    # Extract the values
    xe <- as.numeric(solution_details[[dvmin]][1])
    ye <- as.numeric(solution_details[[dvmin]][2])
    FOCK <- sqrt(as.numeric(solution_details[[dvmin]][3]))




    #now we solve for xc and yc. first define...


    gyfull <- function(y) {
      return(0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }

    gxfull <- function(x) {
      return(0 + 0 + coeffx3 * (x^3) +
               0 + 0 + 0 +
               coeffx2 * (x^2) + 0 + 0 +
               0 + coeffx * (x) + 0 + 0 +
               0 + coeffc)
    }


    #search for 2 roots around xe,ye
    # Number of solutions you want to find
    num_solutions <- 3

    # Initialize an empty list to store solutions
    solutionsxc <- list()
    solutionsyc <- list()

    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 100) {

      # Use multiroot function to find roots
      resultxc <- multiroot(gxfull, xe)
      resultyc <- multiroot(gyfull, ye)

      # Calculate if root is +
      rsignx <- resultxc$root
      rsigny <- resultyc$root

      # Calculate the distance between new solution and existing solutions
      distancesxc <- sapply(solutionsxc, function(sol) sum((sol - resultxc$root)^2))
      distancesyc <- sapply(solutionsyc, function(sol) sum((sol - resultyc$root)^2))

      # If the new solution is significantly different, store it
      if (all(distancesxc > 1e-6) && rsignx > 0) {
        solutionsxc <- c(solutionsxc, list(resultxc$root))
      }

      # If the new solution is significantly different, store it
      if (all(distancesyc > 1e-6) && rsigny > 0) {
        solutionsyc <- c(solutionsyc, list(resultyc$root))
      }


      attempts <- attempts + 1
    }

    yc <- min(as.numeric(solutionsyc))
    xc <- min(as.numeric(solutionsxc))

  }

    #stats on FOCK vector
    rmin <- floor(xe + ye)
    rhoe <- 1 - (a+b)/(a+b+xc)
    rhoc <- 1 - (c+d)/(c+d+yc)
    rhoall <- 1 - (n)/(n + rmin)

    })



     if (vis==1) {



       #tempplotfunc
       # Create a grid of x and y values
       xgrid <- seq(0, 1.05*xc, by = 1.05*xc/250)
       ygrid <- seq(0, 1.05*yc, by = 1.05*yc/250)
       fullgrid <- expand.grid(x = xgrid, y = ygrid)

       z <- matrix(0, nrow = length(xgrid), ncol = length(ygrid))
       for (i in 1:length(xgrid)) {
         for (j in 1:length(ygrid)) {
           if (gfull(xgrid[i], ygrid[j]) > 0) {
             z[i, j] <- 0
           } else {
             z[i, j] <- 1
           }
         }
       }


       nrows <- length(xgrid)
       ncols <- length(ygrid)

       angche <- 180*(atan(ye / xe))/pi


       # Create a data frame for ggplot
       df <- data.frame(
         x = rep(xgrid, each = ncols),
         y = rep(ygrid, times = nrows),
         z = as.vector(z)
       )








       csize <- 1

   plot_obj <- ggplot(df, aes(x = x, y = y)) +
       geom_tile(aes(fill = factor(z))) + guides(fill = "none") +
       scale_fill_manual(values = c("white", "lightblue")) +
       labs(x = "Experimental arm redaction", y = "Control arm redaction", fill = "Z Values", title="Region of retainable redaction (ROAR)") +
       theme_classic() + theme(plot.title = element_text(color = "black", size = 15, face = "bold")) +
       annotate("segment", x = 0, y = 0, xend = 0, yend = yc, colour = "darkgreen", linetype = "dashed", size = 2) +
       annotate("segment", x = 0, y = 0, xend = xc, yend = 0, colour = "blue", linetype = "dashed", size = 2) +
       annotate("segment", x = 0, y = 0, xend = xe, yend = ye, colour = "red", size = 2) +
       annotate("segment", x = 0, y = 0, xend = xe, yend = 0, colour = "orange", linetype = "dotted", size = 1) +
       annotate("segment", x =xe, y = 0, xend = xe, yend = ye, colour = "orange", linetype = "dotted", size = 1) +
       annotate("text",x = xe / 2, y = 1, label = "|xe|", size = 3, color = "black") +
       annotate("text",x = xe + 1, y = ye/2, label = "|ye|", size = 3, color = "black") +
       annotate("text",x = xc/2, y = 1, label = "|xc|", size = 3, color = "black") +
       annotate("text",x = 1, y = yc/2, label = "|yc|", size = 3, color = "black") +
       annotate("text",x = xe / 2, y = ye / 2, label = "FOCK", size = 4, color = "black", angle = angche)

     print(plot_obj)  # Explicitly print the plot

   }



output_list <- list(
  p_value = pval,
  signif = vtest,
  xe = xe,
  ye = ye,
  xe_only = xc,
  yc_only = yc,
  exp_frac = rhoe,
  con_frac = rhoc,
  total_frac = rhoall,
  FOCK = FOCK,
  min_full = rmin + n,
  rmin = rmin)







cat("\n**Total Sample:**", output_list$totalsamp, "\n")
cat("**p-value:**", output_list$p_value, "\n")



cat("\n**FOCK Vector co-ordinates:**\n")
cat("  X-point:", output_list$xe, "\n")
cat("  Y-point:", output_list$ye, "\n")
cat("  Resolved FOCK length: ", output_list$rmin, "subjects" , "\n")


cat("\n**Single arm redacted summary:**\n")
cat("  Minimal Exp Arm only redacted (xi) :", output_list$xe_only, "\n")
cat("  Minimal Con Arm only redacted (yi) :", output_list$yc_only, "\n")



# exp_frac, con_frac, total_frac: Expected, control, and total fractions
cat("\n**Fractions:**\n")
cat("  Experimental arm redacted fraction (%):", 100*output_list$exp_frac, "\n")
cat("  Control arm redacted fraction (%):", 100*output_list$con_frac, "\n")
cat("  Total sample redacted fraction (%):", 100*output_list$total_frac, "\n")




























return(invisible(output_list))


  }



