#' Kolmogorov-Smirnov Goodness of Fit Test for Distribution Generalized Inverse-Lindley
#'
#' @param x is numeric variable
#' @param gamaawal is numeric variable and initial values
#' @param alphacv is level of significance and use 0,05
#' @param epsilon is error and use it epsilon = 0.000001
#'
#' @export
#' @examples gamaawal = c(alpha,lambda)
#' @examples KStestGIL(x, gamaawal, alphacv = 0.05, eps=0.000001)

# Kolmogorov-Smirnov Goodness-of-fit test
KStestGIL <- function(x, gamaawal, alphacv = 0.05, eps=0.000001){
  # Newton Raphson GIL
  {
    n <- length(x)
    diff <- 1
    gama <- gamaawal # Gama is vector of alpha and lambda parameters
    alpha <- gama[1]
    lambda <- gama[2]
    iteration <- 1
    iterationresults <- NULL
    while(diff>eps)
    {
      gamaawal <- gama
      dla <- (n/alpha)-sum(((x^(-alpha))*log(x))/(1+(x^(-alpha))))-sum(log(x))+lambda*sum((x^(-alpha))*log(x))
      dlb <- (2*n/lambda)-(n/(1+lambda))-sum(x^(-alpha))
      d2la2 <- (-n/alpha^2)+sum(((x^(-alpha))*(log(x))^2)/(1+(x^(-alpha)))^2)-lambda*sum((x^(-alpha))*(log(x))^2)
      d2lab <- sum((x^(-alpha))*log(x))
      d2lba <- sum((x^(-alpha))*log(x))
      d2lb2 <- (-n*(lambda^2)-4*n*lambda-2*n)/((lambda^2)*(lambda+1)^2)

      dll <- matrix (c(d2la2, d2lab, d2lba, d2lb2), ncol = 2, byrow = T)
      dl <- matrix (c(dla, dlb), ncol = 1, byrow = T)
      gama <- gama - solve(dll)%*%dl
      alpha <- gama[1]
      lambda <- gama[2]
      diff <- sqrt(sum((gama - gamaawal)^2))
      eachiteration <- c(iteration, gama, diff)
      iterationresults <- rbind(iterationresults, eachiteration)
      iteration <- iteration+1
    }
    colnames(iterationresults) <- c("Iterasi", "alpha", "lambda", "diff")
    # output <- list(param=gama, iterationresults=iterationresults)
  }
  # Kolmogorov Smirnov
  m=sort(unique(unlist(x)),decreasing = FALSE)
  #FS is Cumulative Density Function using unique data
  FS = ((1+lambda+(lambda*m^(-alpha)))/(1+lambda))*exp((-lambda)*m^(-alpha))
  frek <- as.matrix(table(x))
  FT <- cumsum(frek)/length(x)
  D <- max(abs(FT-FS))
  if (alphacv == 0.10){
    critical_values = 1.22/(sqrt(n))
  } else if (alphacv == 0.05){
    critical_values = 1.36/(sqrt(n))
  } else if (alphacv == 0.01){
    critical_values = 1.63/(sqrt(n))
  }
  output=list(D=D,critical_values=critical_values)
  return(output)
} #end of function
