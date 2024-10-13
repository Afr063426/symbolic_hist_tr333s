#' author: Joshua Cervantes Artavia
#' date: 13.10.2024
#' description: This code have as utility to define a Mallows metric for rpart without the use
#' of C code, it is slower in comparation to C code
#' 
#' 
library(RSDA)
library(tidyverse)
histogramas_ejemplo <- iris  %>%
    group_by(Species) %>%
    reframe(Petal.Width = sym.histogram(x = Petal.Width, breaks = pretty(Petal.Width, 5) ))
prueba <- histogramas_ejemplo$Petal.Width[1]

plot(histogramas_ejemplo[1,2])
#' Uniform intervals, it has as objective to uniform the interval 
#' @param x it is a symbolic histogram
#' @param quantile uniform length of the quantiles
inv_quantile <- function(x, prob = 0.1){
    
    breaks <- x$breaks

    cdf <- cumsum(x$breaks)
    idx <- which(cdf >= prob)[1]
    w_l <- cdf[idx]
    w_r <- cdf[idx + 1]
    z_l <- breaks[idx]
    z_r <- breaks[idx + 1]

    z_l + (prob - w_l)/(w_r - w_l)*(z_r - z_l)

}

inv_quantile <- Vectorize(inv_quantile)


#' @param x symbolic histogram object that is going to be uniformed 
#' @param length_quantiles distance beteween the weights
#' @return symbolic histogram uniformed
uniform_histogram <- function(x, length_quantiles = 0.1, props = NULL){

#    x <- histogramas_ejemplo$Petal.Width[1]
    if(is.null(props)){
        w <- seq(0, 1, length_quantiles)
        props <- rep(length_quantiles, length(w)-1)
    }else{
        w <- cumsum(props)
    }
    

    breaks <- inv_quantile(x, w)

    out <- list(
        breaks = breaks, 
        props = props
    )

    vctrs::new_vctr(list(out), class = "symbolic_histogram")
    
}


#Necesito ver como vectorizo esto

#uniform_histogram <- Vectorize(uniform_histogram, vectorize.args = "x")


#' @param y response variable
#' @param offset offset if it is found in the formual, in this case is not used
#' @param parms parameters suplied by the user
#' @param wt weight vector for each observation 
#' 
mallows_initialization <- function(y, offset, parms, wt){
    if(is.matrix(y) && ncol(y) > 1) stop("Matrix response not allowed")

    if (!missing(parms) && length(parms) > 0) warning("parameter argument ignored")


    sfun <- function(yval, dev, wt, ylevel, digits ) {
        paste(" mean=", format(signif(yval, digits)),
        ", RMSE=" , format(signif(dev/wt, digits)),
        sep = '')
    }
    environment(sfun) <- .GlobalEnv

    list(y = c(y), parms = NULL, numresp = 1, numy = 1, summary = sfun)
    
}

#' @param y value of the observation
#' @param wt weight of the observation 
#' @param parms parameteres given by the user
mallows_evaluation <- function(y, wt, parms){
    
    #Diferencias con medida dada por brito

    # wmean <- sum(y*wt)/sum(wt)
    # rss <- sum(wt*(y-wmean)^2)
    #list(label = wmean, deviance = rss)
}



mallows_split <- function(y, wt, x, parms, continuous){
    # Center y
    n <- length(y)
    y <- y- sum(y*wt)/sum(wt)
    if (continuous) {
    # continuous x variable
    temp <- cumsum(y*wt)[-n]
    left.wt <- cumsum(wt)[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp/left.wt
    rmean <- -temp/right.wt
    goodness <- (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2)
    list(goodness = goodness, direction = sign(lmean))
    } else {
    # Categorical X variable
    ux <- sort(unique(x))
    wtsum <- tapply(wt, x, sum)
    ysum <- tapply(y*wt, x, sum)
    means <- ysum/wtsum
    # For anova splits, we can order the categories by their means
    # then use the same code as for a non-categorical
    ord <- order(means)
    n <- length(ord)
    temp <- cumsum(ysum[ord])[-n]
    left.wt <- cumsum(wtsum[ord])[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp/left.wt
    rmean <- -temp/right.wt
    list(goodness= (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2),
    direction = ux[ord])
    }
}


# Ejemplo
# > ulist <- list(eval = mallows_evaluation, split = mallows_split, init = mallows_initialization)
# > fit1 <- rpart(murder ~ population + illiteracy + income + life.exp +
# hs.grad + frost + region, data = mystate,
# method = ulist, minsplit = 10)
# > fit2 <- rpart(murder ~ population + illiteracy + income + life.exp +
# hs.grad + frost + region, data = mystate,
# method = 'anova', minsplit = 10, xval = 0)