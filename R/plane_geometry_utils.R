# The following code is adapted from the PlaneGeometry R package,
# originally authored by Michel van den Bergh.
# Source: https://cran.r-project.org/src/contrib/Archive/PlaneGeometry/
# License: GPL-2

#draw tools


draw <- function(x, ...) {
  UseMethod("draw")
}


#' @method draw Triangle
draw.Triangle <- function(x, ...) {
  A <- x$A; B <- x$B; C <- x$C
  lines(rbind(A, B, C, A), ...)
  invisible()
}


#' @method draw Circle
draw.Circle <- function(x, npoints = 100L, ...) {
  path <- .circlePoints(
    seq(0, 2 * pi, length.out = npoints + 1L)[-1L],
    x$center, x$radius
  )
  polypath(path, ...)
}


#' @method draw Arc
draw.Arc <- function(x, npoints = 100L, ...) {
  alpha1 <- x$alpha1
  alpha2 <- x$alpha2
  if (x$degrees) {
    alpha1 <- alpha1 * pi / 180
    alpha2 <- alpha2 * pi / 180
  }
  path <- .circlePoints(
    seq(alpha1, alpha2, length.out = npoints),
    x$center, x$radius
  )
  lines(path, ...)
}


#' @method draw Ellipse
draw.Ellipse <- function(x, npoints = 100L, ...) {
  alpha <- x$alpha
  if (x$degrees) alpha <- alpha * pi / 180
  path <- .ellipsePoints(
    seq(0, 2 * pi, length.out = npoints + 1L)[-1L],
    x$center, x$rmajor, x$rminor, alpha
  )
  polypath(path, ...)
}


#' @method draw EllipticalArc
draw.EllipticalArc <- function(x, npoints = 100L, ...) {
  lines(x$path(npoints), ...)
}


#' @method draw Line
draw.Line <- function(x, ...) {
  extendA <- x$extendA
  extendB <- x$extendB
  if (extendA && extendB) {
    do <- x$directionAndOffset()
    theta <- do$direction
    offset <- do$offset
    if (sin(theta) != 0) {
      abline(a = offset / sin(theta), b = -1 / tan(theta), ...)
    } else {
      abline(v = offset / cos(theta), ...)
    }
  } else if (extendA) {
    do <- x$directionAndOffset()
    theta <- do$direction
    offset <- do$offset
    A <- x$A; B <- x$B
    bounds <- par("usr")
    if (sin(theta) != 0) {
      curve(offset / sin(theta) - x / tan(theta), add = TRUE, n = 2,
            from = B[1], to = ifelse(A[1] < B[1], bounds[1L], bounds[2L]), ...)
    } else {
      M <- if (A[2] < B[2]) c(A[1], bounds[3L]) else c(A[1], bounds[4L])
      lines(rbind(B, M), ...)
    }
  } else if (extendB) {
    do <- x$directionAndOffset()
    theta <- do$direction
    offset <- do$offset
    A <- x$A; B <- x$B
    bounds <- par("usr")
    if (sin(theta) != 0) {
      curve(offset / sin(theta) - x / tan(theta), add = TRUE, n = 2,
            from = A[1], to = ifelse(A[1] < B[1], bounds[2L], bounds[1L]), ...)
    } else {
      M <- if (A[2] < B[2]) c(A[1], bounds[4L]) else c(A[1], bounds[3L])
      lines(rbind(A, M), ...)
    }
  } else {
    lines(rbind(x$A, x$B), ...)
  }
  invisible()
}


#ELLIPSE

Ellipse <- R6Class(

  "Ellipse",

  private = list(
    .center = c(NA_real_, NA_real_),
    .rmajor = NA_real_,
    .rminor = NA_real_,
    .alpha = NA_real_,
    .degrees = NA
  ),

  active = list(
    #' @field center get or set the center
    center = function(value) {
      if (missing(value)) {
        private[[".center"]]
      } else {
        center <- as.vector(value)
        stopifnot(
          is.numeric(center),
          length(center) == 2L,
          !any(is.na(center)),
          all(is.finite(center))
        )
        private[[".center"]] <- center
      }
    },

    #' @field rmajor get or set the major radius of the ellipse
    rmajor = function(value) {
      if (missing(value)) {
        private[[".rmajor"]]
      } else {
        rmajor <- as.vector(value)
        rminor <- private[[".rminor"]]
        stopifnot(
          is.numeric(rmajor),
          length(rmajor) == 1L,
          !is.na(rmajor),
          is.finite(rmajor),
          rmajor > 0,
          rmajor >= rminor
        )
        private[[".rmajor"]] <- rmajor
      }
    },

    #' @field rminor get or set the minor radius of the ellipse
    rminor = function(value) {
      if (missing(value)) {
        private[[".rminor"]]
      } else {
        rminor <- as.vector(value)
        rmajor <- private[[".rmajor"]]
        stopifnot(
          is.numeric(rminor),
          length(rminor) == 1L,
          !is.na(rminor),
          is.finite(rminor),
          rminor > 0,
          rminor <= rmajor
        )
        private[[".rminor"]] <- rminor
      }
    },

    #' @field alpha get or set the angle of the ellipse
    alpha = function(value) {
      if (missing(value)) {
        private[[".alpha"]]
      } else {
        alpha <- as.vector(value)
        stopifnot(
          is.numeric(alpha),
          length(alpha) == 1L,
          !is.na(alpha),
          is.finite(alpha)
        )
        private[[".alpha"]] <- alpha
      }
    },

    #' @field degrees get or set the \code{degrees} field
    degrees = function(value) {
      if (missing(value)) {
        private[[".degrees"]]
      } else {
        degrees <- as.vector(value)
        stopifnot(
          is.logical(degrees),
          length(degrees) == 1L,
          !is.na(degrees)
        )
        private[[".degrees"]] <- degrees
      }
    }
  ),

  public = list(

    initialize = function(center, rmajor, rminor, alpha, degrees = TRUE) {
      center <- as.vector(center)
      stopifnot(
        is.numeric(center),
        length(center) == 2L,
        !any(is.na(center)),
        all(is.finite(center))
      )
      rmajor <- as.vector(rmajor)
      stopifnot(
        is.numeric(rmajor),
        length(rmajor) == 1L,
        !is.na(rmajor),
        is.finite(rmajor),
        rmajor > 0
      )
      rminor <- as.vector(rminor)
      stopifnot(
        is.numeric(rminor),
        length(rminor) == 1L,
        !is.na(rminor),
        is.finite(rminor),
        rminor > 0,
        rminor <= rmajor
      )
      alpha <- as.vector(alpha)
      stopifnot(
        is.numeric(alpha),
        length(alpha) == 1L,
        !is.na(alpha),
        is.finite(alpha)
      )
      degrees <- as.vector(degrees)
      stopifnot(
        is.logical(degrees),
        length(degrees) == 1L,
        !is.na(degrees)
      )
      private[[".center"]] <- center
      private[[".rmajor"]] <- rmajor
      private[[".rminor"]] <- rminor
      private[[".alpha"]] <- alpha
      private[[".degrees"]] <- degrees
    },


    print = function(...) {
      private[[".center"]] -> center
      private[[".rmajor"]] -> rmajor
      private[[".rminor"]] -> rminor
      private[[".alpha"]] -> alpha
      private[[".degrees"]] -> degrees
      cat("Ellipse:\n")
      cat("       center: ", toString(center), "\n", sep = "")
      cat(" major radius: ", toString(rmajor), "\n", sep = "")
      cat(" minor radius: ", toString(rminor), "\n", sep = "")
      cat("        angle: ",
          sprintf("%s %s", alpha,
                  ifelse(degrees,
                         ifelse(alpha %in% c(0,1,-1), "degree", "degrees"),
                         ifelse(alpha %in% c(0,1,-1), "radian", "radians"))
          ), "\n", sep = "")
    },


    isEqual = function(ell){
      if(is(ell, "Circle")) ell <- .circleAsEllipse(ell)
      private[[".center"]] -> center0
      private[[".rmajor"]] -> rmajor0
      private[[".rminor"]] -> rminor0
      private[[".alpha"]] -> alpha0
      private[[".degrees"]] -> degrees
      if(!degrees) alpha0 <- (alpha0 * 180/pi)
      alpha1 <- ell$alpha
      if(!ell$degrees) alpha1 <- (alpha1 * 180/pi)
      isTRUE(all.equal(
        c(center0, rmajor0, rminor0, alpha0 %% 180),
        c(ell$center, ell$rmajor, ell$rminor, alpha1 %% 180)
      ))
    },


    equation = function(){
      private[[".center"]] -> center
      private[[".rmajor"]]^2 -> a2
      private[[".rminor"]]^2 -> b2
      private[[".alpha"]] -> alpha
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      x <- center[1L]; y <- center[2L]
      sine <- sin(alpha); cosine <- cos(alpha)
      sine2 <- sine*sine; cosine2 <- 1-sine2
      A <- a2*sine2 + b2*cosine2
      B <- 2*(b2-a2)*sine*cosine
      C <- a2*cosine2 + b2*sine2
      D <- -2*A*x - B*y
      E <- -B*x - 2*C*y
      F <- A*x*x + B*x*y + C*y*y - a2*b2
      c(A = A, B = B, C = C, D = D, E = E, F = F)
    },


    includes = function(M){
      ABCDEF <- as.list(self$equation())
      x <- M[1L]; y <- M[2L]
      zero <- with(ABCDEF, A*x*x + B*x*y + C*y*y + D*x + E*y + F)
      isTRUE(all.equal(0, zero, check.attributes=FALSE))
    },


    contains = function(M){
      x <- M[1L]; y <- M[2L]
      ABCDEF <- as.list(self$equation())
      with(ABCDEF, A*x*x + B*x*y + C*y*y + D*x + E*y + F) <= 0
    },


    matrix = function(){
      ABCDEF <- as.list(self$equation())
      X <- with(ABCDEF, cbind(
        c(A, B/2, D/2),
        c(B/2, C, E/2),
        c(D/2, E/2, F)))
      K <- -det(X) / with(ABCDEF, A*C - B*B/4)
      X[-3L,-3L] / K
    },


    path = function(npoints = 100L, closed = FALSE, outer = FALSE){
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      if(outer) {
        .ellipsePointsOuter(
          closed,
          private[[".center"]],
          private[[".rmajor"]],
          private[[".rminor"]],
          alpha,
          as.integer(npoints)
        )
      } else {
        if(closed) {
          t_ <- seq(0, 2*pi, length.out = npoints)
        } else {
          t_ <- seq(0, 2*pi, length.out = npoints+1L)[-1L]
        }
        .ellipsePoints(
          t_,
          private[[".center"]],
          private[[".rmajor"]],
          private[[".rminor"]],
          alpha
        )
      }
    },


    diameter = function(t, conjugate = FALSE){
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      ts <- if(conjugate){
        c(t, t+pi, t+pi/2, t-pi/2)
      }else{
        c(t, t+pi)
      }
      pts <- .ellipsePoints(
        ts,
        private[[".center"]],
        private[[".rmajor"]],
        private[[".rminor"]],
        alpha
      )
      if(conjugate){
        list(
          Line$new(pts[1L,], pts[2L,], FALSE, FALSE),
          Line$new(pts[3L,], pts[4L,], FALSE, FALSE)
        )
      }else{
        Line$new(pts[1L,], pts[2L,], FALSE, FALSE)
      }
    },


    perimeter = function() {
      a <- private[[".rmajor"]]
      b <- private[[".rminor"]]
      4 * a * Re(elliptic_E(pi/2, 1-b^2/a^2, minerror = 1e-12))
    },


    pointFromAngle = function(theta, degrees = TRUE){
      theta <- as.vector(theta)
      stopifnot(
        is.numeric(theta),
        length(theta) >= 1L,
        !any(is.na(theta)),
        all(is.finite(theta))
      )
      O <- private[[".center"]]
      a <- private[[".rmajor"]]
      b <- private[[".rminor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      if(degrees) theta <- theta * pi/180
      #t <- sort(atan2(a * tan(theta), b) %% (2*pi))
      sgn <- ifelse(theta %% (2*pi) <= sqrt(.Machine$double.eps), 1, -1)
      t <- atan2(a/b, 1/tan(theta %% (2*pi))) +
        theta + sgn*sqrt(.Machine$double.eps) -
        (theta + sgn*sqrt(.Machine$double.eps)) %% pi
      # t <- ifelse(theta <= pi,
      #             atan2(a/b,1/tan(theta)),
      #             ifelse(theta <= 2*pi,
      #                    atan2(a/b,1/tan(theta)) + pi,
      #                    atan2(a/b,1/tan(theta)) + 2*pi))
      out <- .ellipsePoints(t, O, a, b, alpha)
      if(length(theta) == 1L) out <- c(out)
      out
    },

    pointFromEccentricAngle = function(t){
      t <- as.vector(t)
      stopifnot(
        is.numeric(t),
        length(t) >= 1L,
        !any(is.na(t)),
        all(is.finite(t))
      )
      O <- private[[".center"]]
      a <- private[[".rmajor"]]
      b <- private[[".rminor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      out <- .ellipsePoints(t, O, a, b, alpha)
      if(length(t) == 1L) out <- c(out)
      out
    },


    semiMajorAxis = function(){
      O <- private[[".center"]]
      a <- private[[".rmajor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      O_A <- a * c(cos(alpha), sin(alpha))
      Line$new(O, O + O_A, FALSE, FALSE)
    },


    semiMinorAxis = function(){
      O <- private[[".center"]]
      b <- private[[".rminor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      O_B <- b * c(-sin(alpha), cos(alpha))
      Line$new(O, O + O_B, FALSE, FALSE)
    },



    foci = function(){
      O <- private[[".center"]]
      a <- private[[".rmajor"]]; b <- private[[".rminor"]]
      e <- sqrt(1 - b*b/a/a)
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      O_A <- a * c(cos(alpha), sin(alpha))
      list(F1 = O + e*O_A, F2 = O - e*O_A)
    },


    tangent = function(t){
      t <- as.vector(t)
      stopifnot(
        is.numeric(t),
        length(t) == 1L,
        !is.na(t),
        is.finite(t)
      )
      O <- private[[".center"]]
      a <- private[[".rmajor"]]; b <- private[[".rminor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      x <- a*cos(t); y <- b*sin(t)
      cosalpha <- cos(alpha); sinalpha <- sin(alpha)
      T <- c(
        O[1L] + cosalpha*x - sinalpha*y,
        O[2L] + sinalpha*x + cosalpha*y
      )
      x <- -a*sin(t); y <- b*cos(t)
      v <- c(
        cosalpha*x - sinalpha*y,
        sinalpha*x + cosalpha*y
      )
      Line$new(T, T+v)
    },


    normal = function(t){
      t <- as.vector(t)
      stopifnot(
        is.numeric(t),
        length(t) == 1L,
        !is.na(t),
        is.finite(t)
      )
      O <- private[[".center"]]
      a <- private[[".rmajor"]]; b <- private[[".rminor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      cosalpha <- cos(alpha); sinalpha <- sin(alpha)
      x <- -a*sin(t); y <- b*cos(t)
      v <- c(
        sinalpha*x + cosalpha*y,
        -cosalpha*x + sinalpha*y
      )
      v / .vlength(v)
    },


    theta2t = function(theta, degrees = TRUE){
      theta <- as.vector(theta)
      stopifnot(
        is.numeric(theta),
        length(theta) == 1L,
        !is.na(theta),
        is.finite(theta)
      )
      a <- private[[".rmajor"]]; b <- private[[".rminor"]]
      if(degrees) theta <- theta * pi/180
      sgn <- ifelse(theta %% (2*pi) <= sqrt(.Machine$double.eps), 1, -1)
      atan2(a/b, 1/tan(theta %% (2*pi))) +
        theta + sgn*sqrt(.Machine$double.eps) -
        (theta + sgn*sqrt(.Machine$double.eps)) %% pi
    },


    regressionLines = function(){
      O <- private[[".center"]]
      a <- private[[".rmajor"]]; b <- private[[".rminor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      cosalpha <- cos(alpha); sinalpha <- sin(alpha)
      A <- -b*sinalpha; B <- -a*cosalpha
      t1 <- .solveTrigonometricEquation(A, B)
      A <- b*cosalpha; B <- -a*sinalpha
      t2 <- .solveTrigonometricEquation(A, B)
      pts <- .ellipsePoints(c(t1,t2), O, a, b, alpha)
      list(
        YonX = Line$new(pts[1L,], pts[2L,], FALSE, FALSE),
        XonY = Line$new(pts[3L,], pts[4L,], FALSE, FALSE)
      )
    },


    boundingbox = function(){
      O <- private[[".center"]]
      a <- private[[".rmajor"]]; b <- private[[".rminor"]]
      alpha <- private[[".alpha"]]
      if(private[[".degrees"]]) alpha <- alpha * pi/180
      cosalpha <- cos(alpha); sinalpha <- sin(alpha)
      A <- -b*sinalpha; B <- -a*cosalpha
      t1 <- .solveTrigonometricEquation(A, B)
      A <- b*cosalpha; B <- -a*sinalpha
      t2 <- .solveTrigonometricEquation(A, B)
      pts <- .ellipsePoints(c(t1,t2), O, a, b, alpha)
      list(
        x = sort(c(pts[1L,1L], pts[2L,1L])),
        y = sort(c(pts[3L,2L], pts[4L,2L]))
      )
    },


    randomPoints = function(n, where = "in"){
      where <- match.arg(where, c("in", "on"))
      S <- self$matrix()
      if(where == "in"){
        sims <- runif_in_ellipsoid(n, S, 1)
        sweep(sims, 2L, private[[".center"]], "+")
      }else{
        sims <- runif_on_ellipsoid(n, S, 1)
        sweep(sims, 2L, private[[".center"]], "+")
      }
    }
  )
)

#' Ellipse from center and matrix
#'

EllipseFromCenterAndMatrix <- function(center, S){
  stopifnot(isSymmetric(S))
  e <- eigen(S, symmetric = TRUE)
  if(any(e$values <= 0)) stop("`S` is not positive.")
  .EllipseFromCenterAndEigen(center, e)
}



GaussianEllipse <- function(mean, Sigma, p){
  stopifnot(
    isSymmetric(Sigma),
    p > 0, p < 1
  )
  e <- eigen(Sigma, symmetric = TRUE)
  if(any(e$values <= 0)) stop("`Sigma` is not positive.")
  r <- -2 * log1p(-p)
  e <- list(
    values = rev(1/e$values)/r,
    vectors = e$vectors %*% cbind(c(0,1),c(-1,0))
  )
  .EllipseFromCenterAndEigen(mean, e)
}


EllipseEquationFromFivePoints <- function(P1, P2, P3, P4, P5){
  P <- rbind(P1, P2, P3, P4, P5)
  if(anyDuplicated(P)) stop("The five points are not distinct.")
  x <- P[,1L]; y <- P[,2L]
  M <- cbind(x*x, x*y, y*y, x, y, 1)
  A <- det(M[,-1L])
  B <- -det(M[,-2L])
  C <- det(M[,-3L])
  if(B*B-4*A*C >= 0) stop("The five points do not lie on an ellipse.")
  D <- -det(M[,-4L])
  E <- det(M[,-5L])
  F <- -det(M[,-6L])
  c(A = A, B = B, C = C, D = D, E = E, F = F)
}



EllipseFromEquation <- function(A, B, C, D, E, F){
  stopifnot(A*C > 0)
  if(B*B-4*A*C >= 0) stop("These parameters do not define an ellipse.")
  #if(D*D + E*E <= 4*(A+C)*F) stop("These parameters do not define an ellipse.")
  Q <- rbind(c(2*A, B, D), c(B, 2*C, E), c(D, E, 2*F))
  if(det(Q) == 0) stop("These parameters do not define an ellipse.")
  M0 <- matrix(c(F, D/2, E/2, D/2, A, B/2, E/2, B/2, C), 3L, 3L)
  M <- matrix(c(A, B/2, B/2, C), 2L, 2L)
  lambda <- eigen(M, symmetric = TRUE)$values
  #if(abs(lambda[1L] - A) >= abs(lambda[1L] - C)) lambda <- rev(lambda)
  detM0 <- det(M0); detM <- det(M)
  a <- sqrt(-detM0 / (detM*lambda[1L]))
  b <- sqrt(-detM0 / (detM*lambda[2L]))
  x <- B*E - 2*C*D
  y <- B*D - 2*A*E
  phi <- if(is.nan(B/(A-C))){
    0
  }else{
    if(abs(C) > abs(A)) atan(B/(A-C))/2 else (pi/2 - atan(-B/(A-C))/2)
  }
  Ellipse$new(c(x,y)/(4*A*C - B*B), max(a,b), min(a,b), (phi*180/pi) %% 180)
}


EllipseFromFivePoints <- function(P1, P2, P3, P4, P5){
  cf <- EllipseEquationFromFivePoints(P1, P2, P3, P4, P5)
  EllipseFromEquation(cf[1L], cf[2L], cf[3L], cf[4L], cf[5L], cf[6L])
}




EllipseFromThreeBoundaryPoints <- function(P1, P2, P3){
  if(.collinear(P1, P2, P3)){
    stop("The three points are collinear.")
  }
  points <- rbind(P1, P2, P3)
  means <- colMeans(points)
  cpoints <- sweep(points, 2L, means)
  S <- 1.5 * solve(crossprod(cpoints))
  EllipseFromCenterAndMatrix(means, S)
}


EllipseFromFociAndOnePoint <- function(F1, F2, P){
  k <- .distance(P, F1) + .distance(P, F2)
  a <- k/2
  center <- (F1 + F2) / 2
  d <- .distance(center, F1)
  b <- sqrt(a*a - d*d)
  alpha <- atan(abs(F1[2]-F2[2])/abs(F1[1]-F2[1]))
  Ellipse$new(center, a, b, alpha, degrees = FALSE)
}


fitEllipse <- function(points){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(points) != 2L){
    stop("The `points` matrix must have two columns.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  fit <- fitConic(points, conicType = "e")
  if(fit[["exitCode"]] != 1){
    stop("The ellipse fitting has failed.", call. = TRUE)
  }
  cfs <- fit[["parA"]]
  fittedEllipse <- EllipseFromEquation(
    cfs[1L], cfs[2L], cfs[3L], cfs[4L], cfs[5L], cfs[6L]
  )
  attr(fittedEllipse, "RSS") <- fit[["RSS"]]
  fittedEllipse
}


maxAreaInscribedEllipse <- function(points, verbose = FALSE) {
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(points) != 2L){
    stop("The `points` matrix must have two columns.", call. = TRUE)
  }
  if(nrow(points) < 3L){
    stop("The `points` matrix must have at least three rows.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  # linear inequalities
  V <- makeV(points)
  H <- scdd(V)[["output"]]
  A <- - H[, -c(1L, 2L)]
  b <- H[, 2L]
  # problem variables
  Bvar <- Variable(2L, 2L, symmetric = TRUE)
  dvar <- Variable(2L)
  # objective
  objective <- Minimize(-log_det(Bvar))
  #constraints
  constraints <- list()
  for(i in 1L:nrow(A)) {
    constraints <- append(
      constraints, list(norm2(Bvar %*% A[i, ]) + sum(A[i, ]*dvar) <= b[i])
    )
  }
  # solve the problem
  program <- Problem(objective, constraints)
  solution <- psolve(program, solver = "SCS", verbose = verbose)
  status <- solution[["status"]]
  if(status != "optimal") {
    warning("Non-optimal solution.")
  }
  # get solutions
  B <- solution$getValue(Bvar)
  d <- c(solution$getValue(dvar))
  # get ellipse
  aff <- Affine$new(B, d)
  unitcircle <- CircleOA(c(0, 0), c(1, 0))
  ell <- aff$transformEllipse(unitcircle)
  attr(ell, "status") <- status
  ell
}

#circles



Circle <- R6Class(

  "Circle",

  private = list(
    .center = c(NA_real_, NA_real_),
    .radius = NA_real_
  ),

  active = list(
    #' @field center get or set the center
    center = function(value) {
      if (missing(value)) {
        private[[".center"]]
      } else {
        center <- as.vector(value)
        stopifnot(
          is.numeric(center),
          length(center) == 2L,
          !any(is.na(center))
        )
        private[[".center"]] <- center
      }
    },

    #' @field radius get or set the radius
    radius = function(value) {
      if (missing(value)) {
        private[[".radius"]]
      } else {
        radius <- as.vector(value)
        stopifnot(
          is.numeric(radius),
          length(radius) == 1L,
          !is.na(radius)
        )
        if(radius < 0) {
          warning("The radius is negative!")
        }
        private[[".radius"]] <- radius
      }
    }
  ),

  public = list(
    #' @description Create a new \code{Circle} object.
    #' @param center the center
    #' @param radius the radius
    #' @return A new \code{Circle} object.
    #' @examples circ <- Circle$new(c(1,1), 1)
    #' circ
    #' circ$center
    #' circ$center <- c(0,0)
    #' circ
    initialize = function(center, radius) {
      center <- as.vector(center)
      stopifnot(
        is.numeric(center),
        length(center) == 2L,
        !any(is.na(center))
      )
      radius <- as.vector(radius)
      stopifnot(
        is.numeric(radius),
        length(radius) == 1L,
        !is.na(radius)
      )
      if(radius < 0) {
        warning("The radius is negative!")
      }
      private[[".center"]] <- center
      private[[".radius"]] <- radius
    },

    #' @description Show instance of a circle object.
    #' @param ... ignored
    #' @examples Circle$new(c(0,0), 2)
    print = function(...) {
      cat("Circle:\n")
      cat(" center: ", toString(private[[".center"]]), "\n", sep = "")
      cat(" radius: ", toString(private[[".radius"]]), "\n", sep = "")
    },

    #' @description Get a point on the reference circle from its polar angle.
    #' @param alpha a number, the angle
    #' @param degrees logical, whether \code{alpha} is given in degrees
    #' @return The point on the circle with polar angle \code{alpha}.
    pointFromAngle = function(alpha, degrees = TRUE) {
      if(degrees) alpha <- alpha * pi/180
      private[[".center"]] + private[[".radius"]] * c(cos(alpha), sin(alpha))
    },

    #' @description Diameter of the reference circle for a given polar angle.
    #' @param alpha an angle in radians, there is one diameter for each value of
    #' \code{alpha} modulo \code{pi}
    #' @return A segment (\code{Line} object).
    #' @examples circ <- Circle$new(c(1,1), 5)
    #' diams <- lapply(c(0, pi/3, 2*pi/3), circ$diameter)
    #' plot(NULL, type="n", asp=1, xlim = c(-4,6), ylim = c(-5,7),
    #'      xlab = NA, ylab = NA)
    #' draw(circ, lwd = 2, col = "yellow")
    #' invisible(lapply(diams, draw, col = "blue"))
    diameter = function(alpha){
      t <- as.vector(alpha)
      stopifnot(
        is.numeric(t),
        length(t) == 1L,
        !is.na(t),
        is.finite(t)
      )
      O <- private[[".center"]]
      v <- private[[".radius"]] * c(cos(t), sin(t))
      Line$new(O+v, O-v, FALSE, FALSE)
    },

    #' @description Tangent of the reference circle at a given polar angle.
    #' @param alpha an angle in radians, there is one tangent for each value of
    #' \code{alpha} modulo \code{2*pi}
    #' @examples circ <- Circle$new(c(1,1), 5)
    #' tangents <- lapply(c(0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3), circ$tangent)
    #' plot(NULL, type="n", asp=1, xlim = c(-4,6), ylim = c(-5,7),
    #'      xlab = NA, ylab = NA)
    #' draw(circ, lwd = 2, col = "yellow")
    #' invisible(lapply(tangents, draw, col = "blue"))
    tangent = function(alpha){
      t <- as.vector(alpha)
      stopifnot(
        is.numeric(t),
        length(t) == 1L,
        !is.na(t),
        is.finite(t)
      )
      r <- private[[".radius"]]
      cost <- cos(t); sint <- sin(t)
      T <- private[[".center"]] + r*c(cost, sint)
      Line$new(T, T + c(-sint,cost))
    },

    #' @description Return the two tangents of the reference circle passing
    #' through an external point.
    #' @param P a point external to the reference circle
    #' @return A list of two \code{Line} objects, the two tangents; the
    #' tangency points are in the \code{B} field of the lines.
    tangentsThroughExternalPoint = function(P){
      P <- as.vector(P)
      stopifnot(
        is.numeric(P),
        length(P) == 2L,
        all(is.finite(P)),
        !any(is.na(P))
      )
      O <- private[[".center"]]
      if(.distance(O,P) <= private[[".radius"]]){
        stop("`P` is not external to the circle.")
      }
      M <- (O+P)/2
      circ <- Circle$new(M, .distance(O, M))
      Is <- intersectionCircleCircle(self, circ)
      list(T1 = Line$new(P, Is[[1L]]), T2 = Line$new(P, Is[[2L]]))
    },

    #' @description Check whether the reference circle equals another circle.
    #' @param circ a \code{Circle} object
    isEqual = function(circ){
      c0 <- private[[".center"]]; r0 <- private[[".radius"]]
      c1 <- circ$center; r1 <- circ$radius
      isTRUE(all.equal(c(c0[1L],c0[2L],r0), c(c1[1L],c1[2L],r1)))
    },

    #' @description Check whether the reference circle differs from another circle.
    #' @param circ a \code{Circle} object
    isDifferent = function(circ){
      !self$isEqual(circ)
    },

    #' @description Check whether the reference circle is orthogonal to a
    #' given circle.
    #' @param circ a \code{Circle} object
    isOrthogonal = function(circ){
      stopifnot(is(circ, "Circle"))
      d2 <- c(crossprod(private[[".center"]]-circ$center))
      R <- private[[".radius"]]
      isTRUE(all.equal(d2, R*R + circ$radius*circ$radius))
    },

    #' @description Angle between the reference circle and a given circle,
    #' if they intersect.
    #' @param circ a \code{Circle} object
    angle = function(circ){
      stopifnot(is(circ, "Circle"))
      center1 <- private[[".center"]]
      center2 <- circ$center
      r1 <- private[[".radius"]]
      r2 <- circ$radius
      d2 <- c(crossprod(center1 - center2))
      epsilon <- sqrt(.Machine$double.eps)
      if(d2 > (r1+r2)^2 + epsilon || d2 < (r1-r2)^2 - epsilon){
        message("The two circles do not intersect.")
        return(NULL)
      }
      b1 <- 1/r1; b2 <- 1/r2
      a1 <- b1*center1; a2 <- b2*center2
      bprime1 <- r1 * (c(crossprod(a1))-1)
      bprime2 <- r2 * (c(crossprod(a2))-1)
      cosTheta <- b1*bprime2/2 + b2*bprime1/2 - .dot(a1,a2)
      acos(cosTheta)
    },

    #' @description Check whether a point belongs to the reference circle.
    #' @param M a point
    includes = function(M){
      isTRUE(all.equal(private[[".radius"]]^2,
                       c(crossprod(M-private[[".center"]]))))
    },

    #' @description Orthogonal circle passing through two points on the reference circle.
    #' @param alpha1,alpha2 two angles defining two points on the reference circle
    #' @param arc logical, whether to return only the arc at the interior of the
    #' reference circle
    #' @return A \code{Circle} object if \code{arc=FALSE}, an \code{Arc} object
    #' if \code{arc=TRUE}, or a \code{Line} object: the diameter
    #' of the reference circle defined by the two points in case when the two
    #' angles differ by \code{pi}.
    #' @examples # hyperbolic triangle
    #' circ <- Circle$new(c(5,5), 3)
    #' arc1 <- circ$orthogonalThroughTwoPointsOnCircle(0, 2*pi/3, arc = TRUE)
    #' arc2 <- circ$orthogonalThroughTwoPointsOnCircle(2*pi/3, 4*pi/3, arc = TRUE)
    #' arc3 <- circ$orthogonalThroughTwoPointsOnCircle(4*pi/3, 0, arc = TRUE)
    #' opar <- par(mar = c(0,0,0,0))
    #' plot(0, 0, type = "n", asp = 1, xlim = c(2,8), ylim = c(2,8))
    #' draw(circ)
    #' draw(arc1, col = "red", lwd = 2)
    #' draw(arc2, col = "green", lwd = 2)
    #' draw(arc3, col = "blue", lwd = 2)
    #' par(opar)
    orthogonalThroughTwoPointsOnCircle = function(alpha1, alpha2, arc = FALSE) {
      I <- private[[".center"]]; r <- private[[".radius"]]
      dalpha <- alpha1 - alpha2
      if(dalpha %% pi == 0){
        eialpha1 <- c(cos(alpha1), sin(alpha1))
        A <- I + r*eialpha1; B <- I - r*eialpha1
        return(Line$new(A, B, !arc, !arc))
      }
      r0 <- r * abs(tan(dalpha/2))
      IO <- r / cos(dalpha/2)
      center <- I + IO * c(cos((alpha1+alpha2)/2), sin((alpha1+alpha2)/2))
      # Oy <- IO * sin((alpha1+alpha2)/2)
      if(arc){
        dalpha <- (alpha2-alpha1)%%(2*pi)# - alpha1%%(2*pi)
        delta <- ifelse(dalpha >= pi, pi, 0)
        beta1 <- -pi/2 + delta
        beta2 <- beta1 - pi + dalpha
        theta1 <- beta1+alpha1 #%% (2*pi)
        theta2 <- beta2+alpha1 #%% (2*pi)
        return(
          Arc$new(center, r0, min(theta1,theta2), max(theta1,theta2), FALSE)
        )
      }
      # Circle$new(I+c(Ox,Oy), r0)
      Circle$new(center, r0)
    },

    #' @description Orthogonal circle passing through two points within the reference circle.
    #' @param P1,P2 two distinct points in the interior of the reference circle
    #' @param arc logical, whether to return the arc joining the two points
    #' instead of the circle
    #' @return A \code{Circle} object or an \code{Arc} object,
    #' or a \code{Line} object if the two points are on a diameter.
    #' @examples circ <- Circle$new(c(0,0),3)
    #' P1 <- c(1,1); P2 <- c(1, 2)
    #' ocirc <- circ$orthogonalThroughTwoPointsWithinCircle(P1, P2)
    #' arc <- circ$orthogonalThroughTwoPointsWithinCircle(P1, P2, arc = TRUE)
    #' plot(0, 0, type = "n", asp = 1, xlab = NA, ylab = NA,
    #'      xlim = c(-3, 4), ylim = c(-3, 4))
    #' draw(circ, lwd = 2)
    #' draw(ocirc, lty = "dashed", lwd = 2)
    #' draw(arc, lwd = 3, col = "blue")
    orthogonalThroughTwoPointsWithinCircle = function(P1, P2, arc = FALSE) {
      if(isTRUE(all.equal(P1, P2, check.attributes = FALSE)))
        stop("`P1` and `P2` must be distinct.")
      I <- private[[".center"]]; r <- private[[".radius"]]; r2 <- r*r
      if(.distance(P1,I) >= r2){
        stop("`P1` is not in the interior of the reference circle.")
      }
      if(.distance(P2,I) >= r2){
        stop("`P2` is not in the interior of the reference circle.")
      }
      if(.collinear(I, P1, P2)){
        return(Line$new(P1, P2, !arc, !arc))
      }
      iota <- Inversion$new(I, r2)
      P1prime <- iota$invert(P1); P2prime <- iota$invert(P2)
      line1 <- Line$new(P1,P1prime); line2 <- Line$new(P2,P2prime)
      perp1 <- suppressMessages(line1$perpendicular((P1+P1prime)/2))
      perp2 <- suppressMessages(line2$perpendicular((P2+P2prime)/2))
      O <- .LineLineIntersection(perp1$A, perp1$B, perp2$A, perp2$B)
      if(arc){
        theta1 <- atan2(P1[2L]-O[2L], P1[1L]-O[1L]) %% (2*pi)
        theta2 <- atan2(P2[2L]-O[2L], P2[1L]-O[1L]) %% (2*pi)
        Arc$new(O, sqrt(c(crossprod(O-P1))),
                min(theta1,theta2), max(theta1,theta2), FALSE)
      }else{
        Circle$new(O, sqrt(c(crossprod(O-P1))))
      }
    },

    #' @description Power of a point with respect to the reference circle.
    #' @param M point
    #' @return A number.
    power = function(M) {
      private[[".radius"]] -> radius
      c(crossprod(M - private[[".center"]])) - radius*radius
    },

    #' @description Radical center of two circles.
    #' @param circ2 a \code{Circle} object
    #' @seealso \code{\link{radicalCenter}} for the radical center of three circles.
    radicalCenter = function(circ2){
      C1 <- private[[".center"]]; C2 <- circ2$center
      k <- private[[".radius"]]^2 - circ2$radius^2;
      C1_C2 <- C2 - C1
      C1C2sqr <- c(crossprod(C1_C2))
      K <- if(C1C2sqr == 0){
        c(Inf, Inf)
      }else{
        (C1+C2)/2 + k/2 * C1_C2/C1C2sqr
      }
      K#/C1[1L] # quid if C1[1] = 0 ?
    },

    #' @description Radical axis of two circles.
    #' @param circ2 a \code{Circle} object
    #' @return A \code{Line} object.
    radicalAxis = function(circ2){
      C1 <- private[[".center"]]; C2 <- circ2$center
      if(isTRUE(all.equal(C1,C2))){
        stop("The two circles must have distinct centers.")
      }
      C1_C2 <- C2 - C1
      v <- c(-C1_C2[2L], C1_C2[1L])
      R <- self$radicalCenter(circ2)
      Line$new(R, R+v, TRUE, TRUE)
      # l <- Line$new(C1, C2, TRUE, TRUE)
      # R <- self$radicalCenter(circ2)
      # l$perpendicular(R, TRUE, TRUE)
    },

    #' @description Rotate the reference circle.
    #' @param alpha angle of rotation
    #' @param O center of rotation
    #' @param degrees logical, whether \code{alpha} is given in degrees
    #' @return A \code{Circle} object.
    rotate = function(alpha, O, degrees = TRUE){
      alpha <- as.vector(alpha)
      stopifnot(
        is.numeric(alpha),
        length(alpha) == 1L,
        !is.na(alpha),
        is.finite(alpha)
      )
      O <- as.vector(O)
      stopifnot(
        is.numeric(O),
        length(O) == 2L,
        !any(is.na(O)),
        all(is.finite(O))
      )
      if(degrees){
        alpha <- alpha * pi/180
      }
      cosalpha <- cos(alpha); sinalpha <- sin(alpha)
      At <- private[[".center"]] - O
      RAt <- c(cosalpha*At[1L]-sinalpha*At[2L], sinalpha*At[1L]+cosalpha*At[2L])
      Circle$new(RAt + O, private[[".radius"]])
    },

    #' @description Translate the reference circle.
    #' @param v the vector of translation
    #' @return A \code{Circle} object.
    translate = function(v){
      v <- as.vector(v)
      stopifnot(
        is.numeric(v),
        length(v) == 2L,
        !any(is.na(v)),
        all(is.finite(v))
      )
      Circle$new(private[[".center"]] + v, private[[".radius"]])
    },

    #' @description Invert the reference circle.
    #' @param inversion an \code{Inversion} object
    #' @return A \code{Circle} object or a \code{Line} object.
    invert = function(inversion){
      inversion$invertCircle(self)
    },

    #' @description Convert the reference circle to an \code{Ellipse} object.
    asEllipse = function(){
      r <- abs(private[[".radius"]])
      Ellipse$new(private[[".center"]], r, r, 0)
    },

    #' @description Random points on or in the reference circle.
    #' @param n an integer, the desired number of points
    #' @param where \code{"in"} to generate inside the circle,
    #' \code{"on"} to generate on the circle
    #' @return The generated points in a two columns matrix with \code{n} rows.
    randomPoints = function(n, where = "in"){
      where <- match.arg(where, c("in", "on"))
      if(where == "in"){
        sims <- runif_in_sphere(n, 2, private[[".radius"]])
        sweep(sims, 2L, private[[".center"]], "+")
      }else{
        sims <- runif_on_sphere(n, 2, private[[".radius"]])
        sweep(sims, 2L, private[[".center"]], "+")
      }
    }
  )
)

#' Radical center

#'
#' @param circ1,circ2,circ3 \code{Circle} objects
#'

radicalCenter <- function(circ1, circ2, circ3){
  l1 <- circ1$radicalAxis(circ2)
  l2 <- circ1$radicalAxis(circ3)
  .LineLineIntersection(l1$A, l1$B, l2$A, l2$B)
}

#' Mid-circle(s)

#'
#' @param circ1,circ2 \code{Circle} objects
#'

midCircles <- function(circ1, circ2){
  stopifnot(
    is(circ1, "Circle"),
    is(circ2, "Circle")
  )

  r1 <- circ1$radius; r2 <- circ2$radius
  O1 <- circ1$center; O2 <- circ2$center

  epsilon <- sqrt(.Machine$double.eps)

  if(r1 == r2){
    if(isTRUE(all.equal(O1,O2))){ # O1=O2
      out <- list(
        C1 = circ1,
        C2 = "circles are equal; every diameter is a mid-circle"
      )
    }else{
      d2 <- c(crossprod(O1 - O2))
      sumRadii2 <- (r1+r2)^2
      I <- (O1 + O2) / 2
      O1_O2 <- O2 - O1
      v <- c(O1_O2[2L], -O1_O2[1L])
      line <- Line$new(I+v, I-v)
      if(d2 < sumRadii2){ # they intersect at two points
        out <- list(
          C1 = Circle$new(I, sqrt(abs(c(crossprod(I-O2)) - r2*r2))),
          C2 = line
        )
      }else{ # they are tangent or they do not intersect
        out <- line
      }
    }
  }else{ # r1 != r2
    d2 <- c(crossprod(O1 - O2))
    sumRadii2 <- (r1+r2)^2
    rho <- r1/r2
    if(d2 > sumRadii2 + epsilon){ # they are outside each other
      I <- O1 - rho/(1-rho)*(O2-O1)
      k <- rho * abs(c(crossprod(I-O2))-r2*r2)
      out <- Circle$new(I, sqrt(k))
    }else if(d2 < (r1-r2)^2 - epsilon){ # one contains the other
      I <- O1 + rho/(1+rho)*(O2-O1)
      k <- rho * abs(c(crossprod(I-O2))-r2*r2)
      out <- Circle$new(I, sqrt(k))
    }else if(sumRadii2 - d2 < epsilon){ # they are externally tangent
      I <- O1 - rho/(1-rho)*(O2-O1)
      k <- rho * abs(c(crossprod(I-O2))-r2*r2)
      out <- Circle$new(I, sqrt(k))
    }else if(d2 - (r1-r2)^2 < epsilon){ # they are internally tangent
      I <- O1 + rho/(1+rho)*(O2-O1)
      k <- rho * abs(c(crossprod(I-O2))-r2*r2)
      out <- Circle$new(I, sqrt(k))
    }else{ # they intersect at two points
      I1 <- O1 - rho/(1-rho)*(O2-O1)
      k1 <- rho * abs(c(crossprod(I1-O2))-r2*r2)
      I2 <- O1 + rho/(1+rho)*(O2-O1)
      k2 <- rho * abs(c(crossprod(I2-O2))-r2*r2)
      out <- list(
        C1 = Circle$new(I1, sqrt(k1)),
        C2 = Circle$new(I2, sqrt(k2))
      )
    }
  }
  out
}


#' @param c0 exterior circle, a \code{Circle} object
#' @param n number of circles, not including the inner circle; at least \code{3}
#' @param phi \code{-1 < phi < 1} controls the radii of the circles
#' @param shift any number; it produces a kind of rotation around the inner
#' circle; values between \code{0} and \code{n} cover all possibilities
#' @param ellipse logical; the centers of the circles of the Steiner chain lie
#' on an ellipse, and this ellipse is returned as an attribute if you set this
#' argument to \code{TRUE}
#'


SteinerChain <- function(c0, n, phi, shift, ellipse = FALSE){
  n <- as.vector(n)
  phi <- as.vector(phi)
  shift <- as.vector(shift)
  stopifnot(
    is(c0, "Circle"),
    .isInteger(n),
    length(n) == 1L,
    n > 2,
    is.numeric(phi),
    length(phi) == 1L,
    !is.na(phi),
    -1 < phi && phi < 1,
    is.numeric(shift),
    length(shift) == 1L,
    !is.na(shift),
    is.finite(shift)
  )
  circles0 <- .SteinerChain_phi0(c0 = c0, n = n, shift = shift)
  if(phi == 0) return(circles0)
  R <- c0$radius; O <- c0$center
  invphi <- 1/phi
  I <- c(R*invphi, 0) + O
  r2 <- R*R * (invphi*invphi-1)
  iota <- Inversion$new(I, r2)
  out <- lapply(circles0, function(circle) iota$invertCircle(circle))
  if(ellipse){
    O2 <- out[[n+1L]]$center
    r <- out[[n+1L]]$radius
    c <- (O2[1L] - O[1L])/2
    a <- (r + R)/2
    b <- sqrt(a*a - c*c)
    attr(out, "ellipse") <- Ellipse$new((O+O2)/2, a, b, 0)
  }
  return(out)
}




#' @param O the center of the circle
#' @param A a point of the circle
#'

CircleOA <- function(O, A){
  Circle$new(O, .distance(O,A))
}


#'
#' @param A,B the endpoints of the diameter
#'

CircleAB <- function(A, B){
  Circle$new((A+B)/2, .distance(A,B)/2)
}

#' Unit circle

unitCircle <- Circle$new(c(0,0), 1)

#other

.isInteger <- function(x){
  is.numeric(x) && all(is.finite(x)) && !any(is.na(x)) && all(trunc(x) == x)
}

.isPoint <- function(M){
  is.numeric(M) && length(M) == 2L && !any(is.na(M)) && all(is.finite(M))
}

.toCplx <- function(M){
  if(isTRUE(all.equal(M, Inf, check.attributes = FALSE))){
    Inf
  }else{
    complex(real = M[1L], imaginary = M[2L])
  }
}

.fromCplx <- function(z){
  c(Re(z), Im(z))
}

.Mod2 <- function(z){
  Re(z)*Re(z) + Im(z)*Im(z)
}

.distance <- function(A, B){
  sqrt(c(crossprod(A-B)))
}

.dot <- function(u, w = NULL){
  c(crossprod(u, w))
}

.vlength <- function(v){
  sqrt(c(crossprod(v)))
}

.isAlmostZero <- function(x){
  isTRUE(all.equal(x + 1, 1, tol = 1e-6))
}

.LineLineIntersection <- function (P1, P2, Q1, Q2) {
  dx1 <- P1[1L] - P2[1L]
  dx2 <- Q1[1L] - Q2[1L]
  dy1 <- P1[2L] - P2[2L]
  dy2 <- Q1[2L] - Q2[2L]
  D <- det(rbind(c(dx1, dy1), c(dx2, dy2)))
  if (D == 0) {
    return(c(Inf, Inf))
  }
  D1 <- det(rbind(P1, P2))
  D2 <- det(rbind(Q1, Q2))
  c(
    det(rbind(c(D1, dx1), c(D2, dx2))),
    det(rbind(c(D1, dy1), c(D2, dy2)))
  ) / D
}

.collinear <- function(A, B, C, tol = 0) {
  notdistinct <-
    isTRUE(all.equal(A, B, check.attributes = FALSE)) ||
    isTRUE(all.equal(A, C, check.attributes = FALSE)) ||
    isTRUE(all.equal(B, C, check.attributes = FALSE))
  if(notdistinct) return(TRUE)
  AB <- B-A; AC <- C-A
  z <- (AB[1] - 1i*AB[2]) * (AC[1] + 1i*AC[2])
  re <- Re(z); im <- Im(z)
  1 / (1 + im*im/re/re) >= 1 - tol
}

.CircleLineIntersection00 <- function(A1, A2, r) {
  x1 <- A1[1L]; y1 <- A1[2L]
  x2 <- A2[1L]; y2 <- A2[2L]
  dx <- x2 - x1; dy <- y2 - y1
  dr2 <- dx*dx + dy*dy
  D <- det(cbind(A1,A2))
  Delta <- r*r*dr2 - D*D
  if(Delta < 0){
    return(NULL)
  }
  if(Delta < sqrt(.Machine$double.eps)){
    return(D/dr2 * c(dy, -dx))
  }
  sgn <- ifelse(dy < 0, -1, 1)
  Ddy <- D*dy
  sqrtDelta <- sqrt(Delta)
  I1 <- c(
    Ddy + sgn*dx * sqrtDelta,
    -D*dx + abs(dy)*sqrtDelta
  ) / dr2
  I2 <- c(
    Ddy - sgn*dx * sqrtDelta,
    -D*dx - abs(dy)*sqrtDelta
  ) / dr2
  list(I1 = I1, I2 = I2)
}

.SteinerChain_phi0 <- function(c0, n, shift){
  R <- c0$radius; O <- c0$center
  sine <- sin(pi/n)
  Cradius <- R / (1+sine)
  Cside <- Cradius*sine
  circles0 <- vector("list", n+1)
  for(i in 1:n){
    beta <- (i-1+shift)*2*pi/n
    pti <- Cradius*c(cos(beta), sin(beta)) + O
    circ1 <- Circle$new(pti, Cside)
    circles0[[i]] <- circ1
  }
  circles0[[n+1]] <- Circle$new(O, R-2*Cside)
  circles0
}

.inversion2conjugateMobius <- function(iota){
  C <- .toCplx(iota$pole)
  k <- iota$power
  Mobius$new(rbind(c(C, k-C*Conj(C)), c(1, -Conj(C))))
}

.MobiusMappingThreePoints2ZeroOneInf <- function(z1, z2, z3){
  if(z1 == Inf){
    K <- z2 - z3
    return(Mobius$new(rbind(c(0,K),c(1,-z3))))
  }
  if(z2 == Inf){
    return(Mobius$new(rbind(c(1,-z1),c(1,-z3))))
  }
  if(z3 == Inf){
    K <- 1 / (z2 - z1)
    return(Mobius$new(rbind(c(K, -K*z1),c(0,1))))
  }
  K <- (z2 - z3) / (z2 - z1)
  Mobius$new(rbind(c(K, -K*z1),c(1,-z3)))
}

.ellipsePoints <- function(t, O, a, b, alpha){
  x <- a*cos(t); y <- b*sin(t)
  cosalpha <- cos(alpha); sinalpha <- sin(alpha)
  cbind(
    x = O[1L] + cosalpha*x - sinalpha*y,
    y = O[2L] + sinalpha*x + cosalpha*y
  )
}

.ellipsePointsOuter <- function(closed, O, a, b, alpha, n) {
  x0 <- O[1L]; y0 <- O[2L]
  if(!closed) {
    n <- n + 1L
  }
  theta <- c(seq(0, 2 * pi, length.out = n), 0)
  sintheta <- sin(theta); costheta <- cos(theta)
  cosalpha <- cos(alpha); sinalpha <- sin(alpha)
  slopes <- (-a * sintheta * sinalpha + b * costheta * cosalpha) /
    (-a * sintheta * cosalpha - b * costheta * sinalpha)
  crds <- cbind(
    a * costheta * cosalpha - b * sintheta * sinalpha + x0,
    a * costheta * sinalpha + b * sintheta * cosalpha + y0
  )
  intercepts <- crds[, 2L] - slopes * crds[, 1L]
  i <- 1L:(n-1L)
  x <- (intercepts[i] - intercepts[i+1L]) / (slopes[i+1L] - slopes[i])
  y <- slopes[i]*x + intercepts[i]
  out <- cbind(x = x, y = y)
  if(closed) {
    out <- rbind(out, c(x = x[1L], y = y[1L]))
  }
  out
}

.circlePoints <- function(t, O, r){
  cbind(
    x = O[1] + r*cos(t),
    y = O[2] + r*sin(t)
  )
}

.solveTrigonometricEquation <- function(a, b, D = 0){
  # solve a*cos(x) + b*sin(x) = D
  if(D == 0){
    return((atan2(b, a) + c(pi, -pi)/2) %% (2*pi))
  }
  d <- sqrt(a*a+b*b)
  if(abs(D)/d > 1) return(NULL)
  if(D == d) return(atan2(b,a) %% (2*pi))
  if(D == -d) return((atan2(b,a)+pi) %% (2*pi))
  e <- acos(D/d)
  (atan2(b,a) + c(e, -e)) %% (2*pi)
  # https://socratic.org/questions/how-do-you-use-linear-combinations-to-solve-trigonometric-equations
  # D = sqrt(a*a+b*b) * cos(x - atan2(b,a))
  # cos(x - atan2(b,a)) = D / sqrt(a*a+b*b)
  # x - atan2(b,a) = acos(D / sqrt(a*a+b*b)) or -acos(D / sqrt(a*a+b*b))
}

.circleAsEllipse <- function(circ){
  Ellipse$new(circ$center, circ$radius, circ$radius, 0)
}

.EllipseFromCenterAndEigen <- function(center, e){
  v <- e$vectors[,2L]
  alpha <- (atan2(v[2L],v[1L]) * 180/pi) %% 180
  a <- 1/sqrt(e$values[2L])
  b <- 1/sqrt(e$values[1L])
  Ellipse$new(center, a, b, alpha)
}

# .Jordan2x2 <- function(M){
#   detM <- det(M)
#   trM <- M[1L,1L] + M[2L,2L]
#   if(abs(trM*trM - 4*detM) < sqrt(.Machine$double.eps)){
#     lambda <- trM/2
#     if(isTRUE(all.equal(M, diag(c(lambda,lambda))))){
#       list(P = diag(2L), J = M)
#     }else{
#       N <- M - diag(c(lambda,lambda))
#       if(isTRUE(all.equal(N[,2L], c(0,0)))){
#         v2 <- c(1,0)
#         v1 <- N[,1L]
#       }else{
#         v2 <- c(0,1)
#         v1 <- N[,2L]
#       }
#       list(P = cbind(v1,v2), J = rbind(c(1,lambda),c(0,1)))
#     }
#   }else{
#     eig <- eigen(M)
#     list(P = eig$vectors, J = diag(eig$values))
#   }
# }

`%**%` <- function(M, k){
  Reduce(`%*%`, replicate(k, M, simplify = FALSE))
}

.htrigonometricEquation <- function(a, b, D) {
  # solution of a*cosh(x) + b*sinh(x) = D
  a2 <- a * a
  b2 <- b * b
  if(a2 > b2) {
    acosh(D / (a * sqrt(1 - b2/a2))) - atanh(b/a)
  } else if(a2 < b2) {
    asinh(D / (b * sqrt(1 - a2/b2))) - atanh(a/b)
  } else if(a == b) {
    log(D/a)
  } else {
    log(-D/a)
  }
}

.good_t <- function(H, xmin, xmax, ymin, ymax) {
  OAB <- H$OAB()
  O <- OAB[["O"]]
  A <- OAB[["A"]]
  B <- OAB[["B"]]
  L1 <- H$L1
  L2 <- H$L2
  theta1 <- L1$directionAndOffset()$direction
  theta2 <- L2$directionAndOffset()$direction
  sgn <- if(theta2 > pi) -1 else 1
  t1 <- .htrigonometricEquation(sgn*A[1L], B[1L], xmin - O[1L]) # prendre -g1 car branche ouest
  if(is.nan(t1)) t1 <- 0                                       # -> voir ça selon la direction de l1 ou l2
  t2 <- .htrigonometricEquation(A[1L], B[1L], xmax - O[1L])
  if(is.nan(t2)) t2 <- 0
  sgn <- if(theta1 > pi) 1 else -1
  t3 <- .htrigonometricEquation(sgn*A[2L], B[2L], ymin - O[2L])
  if(is.nan(t3)) t3 <- 0
  t4 <- .htrigonometricEquation(A[2L], B[2L], ymax - O[2L])
  if(is.nan(t4)) t4 <- 0
  print(c(t1, t2, t3, t4))
  #min(max(t1, t2), max(t3, t4)) # take absolute values?
  min(max(abs(t1), abs(t2)), max(abs(t3), abs(t4)))
}
