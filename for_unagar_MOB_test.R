mob_test <- function(formula, PATH=PATH, trim = .1, minsplit = 5, breakties = FALSE, CvM = FALSE){
  #install.packages("modeltools")
  #install.packages("strucchange")
  #install.packages("sandwich")
  data <- read.csv(PATH)
  library(modeltools)
  library(sandwich)
  library(strucchange)
  formula <- as.formula(formula)
  # Trying MOB model
  #model = mob(MV ~ LSTAT + RM | ZN + INDUS + CHAS + NOX + AGE + DIS + RAD + TAX + CRIM, data = boston, control = ctrl, model = linearModel)
  # Some function for formula (formula is used as the fm into the other functions)
  mobpp <- function(formula, data, model) {
    ff <- attr(ParseFormula(formula), "formula")
    ff$input[[3]] <- ff$input[[2]]
    ff$input[[2]] <- ff$response[[2]]
    dpp(model, as.formula(ff$input), other = list(part = as.formula(ff$blocks)), 
        data = data)
    }
  mf <- mobpp(formula, data = data, model = linearModel)
  # Extracting the data(independant varriable)  
  subset <- mf@get("part")
  
  fm <- fit(linearModel, mf)
  ## extract estimating functions  
  process <- as.matrix(estfun(fm))
  k <- NCOL(process)
  
  partvar <- mf@get("part")
  m <- NCOL(partvar)
  pval <- rep.int(0, m)
  stat <- rep.int(0, m)
  ifac <- rep.int(FALSE, m)
  
  n <- NROW(partvar)
  ## scale process
  process <- process/sqrt(n)
  J12 <- root.matrix(crossprod(process))
  process <- t(chol2inv(chol(J12)) %*% t(process))  
  
  from <- if(trim > 1) trim else ceiling(n * trim)
  from <- max(from, minsplit)
  to <- n - from
  lambda <- ((n-from)*to)/(from*(n-to))
  
  beta <- get("sc.beta.sup")
  logp.supLM <- function(x, k, lambda)
    {
    if(k > 40) {
      ## use Estrella (2003) asymptotic approximation
      logp_estrella2003 <- function(x, k, lambda)
        -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * (1 - k/x) + 2/x))
      ## FIXME: Estrella only works well for large enough x
      ## hence require x > 1.5 * k for Estrella approximation and
      ## use an ad hoc interpolation for larger p-values
      p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
      } else {
        ## use Hansen (1997) approximation
        m <- ncol(beta)-1
        if(lambda<1) tau <- lambda
        else tau <- 1/(1+sqrt(lambda))
        beta <- beta[(((k-1)*25 +1):(k*25)),]
        dummy <- beta[,(1:m)]%*%x^(0:(m-1))
        dummy <- dummy*(dummy>0)
        pp <- pchisq(dummy, beta[,(m+1)], lower.tail = FALSE, log.p = TRUE)
        if(tau==0.5)
          p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
        else if(tau <= 0.01)
          p <- pp[25]
        else if(tau >= 0.49)
          p <- log((exp(log(0.5-tau) + pp[1]) + exp(log(tau-0.49) + pchisq(x,k,lower.tail = FALSE, log.p = TRUE)))*100)
        else
          {
            taua <- (0.51-tau)*50
            tau1 <- floor(taua)
            p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1+1]))
          }
        }
    return(as.vector(p))
    }
  for(i in 1:m) {
    pvi <- partvar[,i]
    if(is.factor(pvi)) {
      proci <- process[ORDER(pvi), , drop = FALSE]
      ifac[i] <- TRUE
      # re-apply factor() added to drop unused levels
      pvi <- factor(pvi[ORDER(pvi)])
      # compute segment weights
      segweights <- as.vector(table(pvi))/n ## tapply(ww, pvi, sum)/n      
      # compute statistic only if at least two levels are left
      if(length(segweights) < 2) {
        stat[i] <- 0
        pval[i] <- NA
        } else {      
          stat[i] <- sum(sapply(1:k, function(j) (tapply(proci[,j], pvi, sum)^2)/segweights))
          pval[i] <- pchisq(stat[i], k*(length(levels(pvi))-1), log.p = TRUE, lower.tail = FALSE)
          }
      } else {
        oi <- if(breakties) {
          mm <- sort(unique(pvi))
          mm <- ifelse(length(mm) > 1, min(diff(mm))/10, 1)
          ORDER(pvi + runif(length(pvi), min = -mm, max = +mm))
          } else {
            ORDER(pvi)
            }    
        proci <- process[oi, , drop = FALSE]
        proci <- apply(proci, 2, cumsum)
        stat[i] <- if(CvM) sum((proci)^2)/n 
        else if(from < to) {
          xx <- rowSums(proci^2)
          xx <- xx[from:to]
          tt <- (from:to)/n
          max(xx/(tt * (1-tt)))	  
          } else {
            0
            }
        pval[i] <- if(CvM) log(approx(c(0, critval), c(1, 1-as.numeric(names(critval))), stat[i], rule=2)$y)
        else if(from < to) logp.supLM(stat[i], k, lambda) else NA
      }
  }

  ## select variable with minimal p-value
  
  
  # bonferroni adjustment
  #pval1 <- pmin(1, sum(!is.na(pval)) * pval)
  #pval2 <- 1 - (1-pval)^sum(!is.na(pval))
  #pval <- ifelse(!is.na(pval) & (pval > 0.01), pval2, pval1)
  
  pval[is.na(pval)] <- 1
  best <- which.min(pval)
  if(length(best) < 1) best <- NA
  rval <- list(pval = exp(pval), stat = stat, best = best)
  names(rval$pval) <- names(partvar)
  names(rval$stat) <- names(partvar)
  
  if (!all(is.na(rval$best)) && min(exp(pval)) <= .01){
    names(rval$best) <- names(partvar)[rval$best]
  }
  else {
    names(rval$best) <- "NA"
    rval$best <- -1
  }
  return (rval$best)
}
  