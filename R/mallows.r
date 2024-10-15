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


mean.symbolic_histogram <- function(x){
    breaks <- x$breaks 

    breaks_mean <- Reduce(`+`, breaks)/length(breaks[[1]])
    props <- x$props[1]
    out <- list(
        breaks = breaks_mean, 
        props = props
    )
    vctrs::new_vctr(list(out), class = "symbolic_histogram")   
}

mse_histogram <- function(y, pred){
    if(length(y) > 1){
        breaks <- y$breaks
        centers <- lapply(breaks, function(x) na.omit((x+lag(x))/2))
        centers[[1]]
        ranges <- lapply(breaks, function(x) na.omit(x-lag(x)))
        
        breaks_mean <- pred$breaks

        centers_mean <- (breaks_mean + lag(breaks_mean))/2
        ranges_mean <- (breaks_mean - lag(breaks_mean))
        centers_mean <- centers_mean[!is.na(centers_mean)]
        ranges_mean <- ranges_mean[!is.na(ranges_mean)]

        props <- y$props
        
        sum_centers <- sapply(centers, function(x) sum(props[[1]][1]*(x - centers_mean)^2))
        sum_ranges <- 1/3*sapply(ranges, function(x) sum(props[[1]][1]*(x - ranges_mean)^2))
        return(sum_centers + sum_ranges)
    }else{
        return(0)
    }
    
}

histogramas_ejemplo <- histogramas_ejemplo %>%
    group_by(Species) %>%
    mutate(histograma = uniform_histogram(Petal.Width)) 
histogramas_ejemplo$prueba <- c(1,2,3)
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
        paste(" hist=", format(signif(yval, digits)),
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
    y <- histogramas_ejemplo$histograma
    
    #Diferencias con medida dada por brito
    breaks <- y$breaks 

    breaks_mean <- Reduce(`+`, breaks)/length(breaks[[1]])
    centers <- lapply(breaks, function(x) na.omit((x+lag(x))/2))
    ranges <- lapply(breaks, function(x) na.omit(x-lag(x)))

    centers_mean <- (breaks_mean + lag(breaks_mean))/2
    ranges_mean <- (breaks_mean - lag(breaks_mean))
    centers_mean <- centers_mean[!is.na(centers_mean)]
    ranges_mean <- ranges_mean[!is.na(ranges_mean)]

    props <- y$props
    
    sum_centers <- sapply(centers, function(x) sum(props[[1]][1]*(x - centers_mean)^2))*wt
    sum_ranges <- 1/3*sapply(ranges, function(x) sum(props[[1]][1]*(x - ranges_mean)^2))*wt
    deviance <- (sum(sum_centers) + sum(sum_ranges))


    # izquierdo <- mean_centers - mean_ranges
    # derecho <-  mean_centers + mean_ranges
    out <- list(
        breaks = breaks_mean, 
        props = props[[1]]
    )

    
    # rss <- sum(wt*(y-wmean)^2)
    list(label = vctrs::new_vctr(list(out), class = "symbolic_histogram"), deviance = deviance)
}




mallows_split <- function(y, wt = NULL, x, parms = NULL, continuous = TRUE){
    # Center y
    
    if (continuous) {
    mse <- Inf
    x <- histogramas_ejemplo$prueba
        for(i in 1:(length(x) - 1)){
            
            left <- which(x[i] >= x)

            right <- which(x[i] < x)
            direction <- ifelse(x[i] >= x, -1, 1)

            left_mean <- mean(y[left])
            right_mean <- mean(y[right])
            left_mse <- mse_histogram(y[left], left_mean)
            right_mse <- mse_histogram(y[right], right_mean)

            mse_aux <- sum((sum(x <= x[i])*left_mse + sum(x > x[i])*right_mse))/length(x)
            print(mse_aux)
            if(mse > mse_aux){
                df_aux <- data.frame(mse = c(left_mse, right_mse), indice = c(left, right)) %>%
                    arrange(indice)
                mse <- mse_aux
                goodness <- -df_aux$mse
                split <- direction
            }

        }
    
    
    }
    list(goodness = goodness, direction = split)
}



# ------------------------------------------------
# ---------------- PRUEBA -----------------------

ulist <- list(eval = mallows_evaluation, split = mallows_split, init = mallows_initialization)
histogramas_ejemplo_aux <- histogramas_ejemplo %>%
    ungroup %>%
    mutate(texto = "a", collapse = ":") %>% dplyr::select(texto, prueba)

prueba <- rpart(texto ~ prueba , method = ulist, data = histogramas_ejemplo_aux )
formula <- histograma ~ prueba

# Ejemplo
# > ulist <- list(eval = mallows_evaluation, split = mallows_split, init = mallows_initialization)
# > fit1 <- rpart(murder ~ population + illiteracy + income + life.exp +
# hs.grad + frost + region, data = mystate,
# method = ulist, minsplit = 10)
# > fit2 <- rpart(murder ~ population + illiteracy + income + life.exp +
# hs.grad + frost + region, data = mystate,
# method = 'anova', minsplit = 10, xval = 0)


