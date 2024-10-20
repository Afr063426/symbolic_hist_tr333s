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
    reframe(mediana = median(Petal.Width), 
        Petal.Width = sym.histogram(x = Petal.Width, breaks = pretty(Petal.Width, 5) )
    )
prueba <- histogramas_ejemplo$Petal.Width[1]


#' Uniform intervals, it has as objective to uniform the interval 
#' @param x it is a symbolic histogram
#' @param quantile uniform length of the quantiles
inv_quantile <- function(x, prob = 0.1){
    
    breaks <- x$breaks
    props <- x$props 

    if(prob != 1){
        
        w_i <- c(0, cumsum(props))
        idx <- which(prob < w_i)[1]
        idx <- idx - 1
        w_l <- w_i[idx]
        w_r <- w_i[idx + 1]

        z_l <- breaks[idx]
        z_r <- breaks[idx + 1]

        z_l + (prob - w_l)/(w_r - w_l)*(z_r - z_l)
    }else{
        max(breaks)
    }
    

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

    breaks_mean <- Reduce(`+`, breaks)/length(x)
    props <- x$props[1]
    out <- list(
        breaks = breaks_mean, 
        props = props
    )
    vctrs::new_vctr(list(out), class = "symbolic_histogram")   
}

median.symbolic_histogram <- function(x){
    
    breaks <- x$breaks 

    props <- x$props


    w_i <- c(0, cumsum(x$props))
    idx <- which(w_i > 0.5)[1]
    idx <- idx-1

    
    q <- breaks[idx] + 0.5/(w_i[idx+1] - w_i[idx]) * (breaks[idx+1] - breaks[idx])

    q

}
median.symbolic_histogram <- Vectorize(median.symbolic_histogram)

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




mallows_split <- function(x, y){
    # Center y
    
    # if (continuous) {
    mse <- Inf
    x_aux <- sort(unique(x))
    if(length(x_aux) > 1){

        for(i in 1:(length(x_aux) - 1)){
            
            left <- which(x_aux[i] >= x)

            right <- which(x_aux[i] < x)
            #direction <- ifelse(x[i] >= x, -1, 1)

            left_mean <- mean(y[left])
            right_mean <- mean(y[right])
            
            left_mse <- mse_histogram(y[left], left_mean)
            right_mse <- mse_histogram(y[right], right_mean)

            mse_aux <- sum((sum(x <= x_aux[i])*left_mse + sum(x > x_aux[i])*right_mse))/length(x)
            print(mse_aux)
            if(mse > mse_aux){
                # df_aux <- data.frame(mse = c(left_mse, right_mse), indice = c(left, right)) %>%
                #     arrange(indice)
                mse <- mse_aux
                #goodness <- -df_aux$mse
                #split <- direction
                split <- x_aux[i]
            }

        }
    }else{
        mse <- Inf
        split <- NULL    
    }
    
    
    #}
    return(c(sse =  mse, split = split))
    #list(goodness = mse, direction = split)
}


model.matrix.symbolic_histogram <- function(formula, data){
  
  # extract target
  y <- data[, as.character(formula)[2]]

  data[as.character(formula)[2]] <- 0
  
  # get the design matrix
  X <- model.matrix(formula, data)

  #X[as.character(formula)[2]] <- y

  X
  
}

#' reg_tree
#' Fits a simple regression tree with SSE splitting criterion. The estimator function
#' is the mean.
#' 
#' @param formula an object of class formula
#' @param data a data.frame or matrix
#' @param minsize a numeric value indicating the minimum size of observations
#'                in a leaf
#'
#' @return \itemize{
#' \item tree - the tree object containing all splitting rules and observations
#' \item fit - our fitted values, i.e. X %*% theta
#' \item formula - the underlying formula
#' \item data - the underlying data
#' }
#' @export
#'
#' @examples # Complete runthrough see: www.github.com/andrebleier/cheapml
reg_tree <- function(formula, data, minsize) {
  
  # coerce to data.frame
  data <- as.data.frame(data)
  
  # handle formula
  formula <- terms.formula(formula)
  

  # get the design matrix
  X <- model.matrix.symbolic_histogram(formula, data)
  
  # extract target
  y <- data[, as.character(formula)[2]]
  
  # initialize while loop
  do_splits <- TRUE
  
  # create output data.frame with splitting rules and observations
  tree_info <- data.frame(NODE = 1, NOBS = nrow(data), FILTER = NA, TERMINAL = "SPLIT",
                          stringsAsFactors = FALSE)
  
  # keep splitting until there are only leafs left
  while(do_splits) {
    
    # which parents have to be splitted
    to_calculate <- which(tree_info$TERMINAL == "SPLIT")
    
    for (j in to_calculate) {
      
      # handle root node
      if (!is.na(tree_info[j, "FILTER"])) {
        # subset data according to the filter
        this_data <- subset(data, eval(parse(text = tree_info[j, "FILTER"])))
        # get the design matrix
        X <- model.matrix.symbolic_histogram(formula, this_data)
      } else {
        this_data <- data
      }
      
      # estimate splitting criteria
      splitting <- apply(X,  MARGIN = 2, FUN = mallows_split, y = this_data[, all.vars(formula)[1]])

      # get the min SSE
      tmp_splitter <- which.min(sapply(splitting, `[`, 1))
      
      # define maxnode
      mn <- max(tree_info$NODE)
      
      # paste filter rules
      tmp_filter <- c(paste(gsub(".sse", "",names(tmp_splitter)), ">=", 
                            splitting[[tmp_splitter]][2]),
                      paste(gsub(".sse", "",names(tmp_splitter)), "<", 
                            splitting[[tmp_splitter]][2]))
      
      # Error handling! check if the splitting rule has already been invoked
      split_here  <- !sapply(tmp_filter,
                             FUN = function(x,y) any(grepl(x, x = y)),
                             y = tree_info$FILTER)
      
      # append the splitting rules
      if (!is.na(tree_info[j, "FILTER"])) {
        tmp_filter  <- paste(tree_info[j, "FILTER"], 
                             tmp_filter, sep = " & ")
      } 
      
      # get the number of observations in current node
      tmp_nobs <- sapply(tmp_filter, 
                         FUN = function(i, x) {
                           nrow(subset(x = x, subset = eval(parse(text = i))))
                         },
                         x = this_data)  
      
      # insufficient minsize for split
      if (any(tmp_nobs <= minsize)) {
        split_here <- rep(FALSE, 2)
      }
      
      # create children data frame
      children <- data.frame(NODE = c(mn+1, mn+2),
                             NOBS = tmp_nobs,
                             FILTER = tmp_filter,
                             TERMINAL = rep("SPLIT", 2),
                             row.names = NULL)[split_here,]
      
      # overwrite state of current node
      tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "PARENT")
       
      # bind everything
      tree_info <- rbind(tree_info, children)
      
      # check if there are any open splits left
      do_splits <- !all(tree_info$TERMINAL != "SPLIT")
    } # end for
  } # end while
  
  # calculate fitted values
  leafs <- tree_info[tree_info$TERMINAL == "LEAF", ]
  fitted <- c()
  for (i in seq_len(nrow(leafs))) {
    # extract index
    ind <- as.numeric(rownames(subset(data, eval(parse(text = leafs[i, "FILTER"])))))
    # estimator is the mean y value of the leaf
    fitted[ind] <- mean(y[ind])
  }
  
  # return everything
  return(list(tree = tree_info, fit = fitted, formula = formula, data = data))
}


regression_prueba <- reg_tree(formula = histograma ~ mediana, data = histogramas_ejemplo, minsize = 1)

# ------------------------------------------------
# ---------------- PRUEBA -----------------------


ulist <- list(eval = mallows_evaluation, split = mallows_split, init = mallows_initialization)
histogramas_ejemplo_aux <- histogramas_ejemplo %>%
    ungroup %>%
    mutate(texto = "a", collapse = ":") %>% dplyr::select(texto, prueba)
library(rpart)
matriz_prueba <- model.matrix(mediana ~ prueba, histogramas_ejemplo)
matriz_prueba <- model.frame.default(mediana ~ Species, histogramas_ejemplo)
matriz_prueba$mediana <- histogramas_ejemplo$histograma
class(matriz_prueba$mediana)
as.numeric.symbolic_histogram <- function(x){
    x
}
as.numeric.symbolic_histogram(matriz_prueba$mediana)
prueba <- rpart(histograma ~ prueba, method = ulist, model = matriz_prueba)
formula <- histograma ~ prueba

# Ejemplo
# > ulist <- list(eval = mallows_evaluation, split = mallows_split, init = mallows_initialization)
# > fit1 <- rpart(murder ~ population + illiteracy + income + life.exp +
# hs.grad + frost + region, data = mystate,
# method = ulist, minsplit = 10)
# > fit2 <- rpart(murder ~ population + illiteracy + income + life.exp +
# hs.grad + frost + region, data = mystate,
# method = 'anova', minsplit = 10, xval = 0)


# Se crea la tabla de pacientes de Diday y Billard
pulso_cardiaco = c(
    vctrs::new_vctr(list(list(breaks = c(44,60,68), props = c(0.8, 0.2))), class = "symbolic_histogram"), 
    vctrs::new_vctr(list(list(breaks = c(60,70,72), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(56,80,90), props = c(0.6, 0.4))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(70,75,112), props = c(0.4, 0.6))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(54,56,72), props = c(0.2, 0.8))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(70,80,100), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(63,73,75), props = c(0.4, 0.6))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(72,79,100), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(76,80,98), props = c(0.2, 0.8))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(86,94,96), props = c(0.8, 0.2))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(86,89,89), props = c(0.6, 0.4))), class = "symbolic_histogram")
    )

presion_sistolica <- c(
    vctrs::new_vctr(list(list(breaks = c(90,95,100), props = c(0.2, 0.8))), class = "symbolic_histogram"), 
    vctrs::new_vctr(list(list(breaks = c(90,110,110), props = c(0.4, 0.6))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(140,160,180), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(110,120,142), props = c(0.2, 0.8))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(90,98,100), props = c(0.6, 0.4))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(130,150,160), props = c(0.4, 0.6))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(140,145,150), props = c(0.2, 0.8))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(130,140,160), props = c(0.4, 0.8))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(110,160,190), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(138,142,180), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(110,135,150), props = c(0.2, 0.8))), class = "symbolic_histogram")
    )

presion_diastolica <- c(
    vctrs::new_vctr(list(list(breaks = c(50,60,70), props = c(0.4, 0.6))), class = "symbolic_histogram"), 
    vctrs::new_vctr(list(list(breaks = c(70,80,90), props = c(0.2, 0.8))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(90,92,100), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(80,85,108), props = c(0.6, 0.4))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(50,63,70), props = c(0.4, 0.6))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(80,90,100), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(60,80,100), props = c(0.2, 0.8))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(76,85,90), props = c(0.5, 0.5))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(70,100,110), props = c(0.4, 0.6))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(90,100,110), props = c(0.4, 0.6))), class = "symbolic_histogram"),
    vctrs::new_vctr(list(list(breaks = c(78,88,100), props = c(0.2, 0.8))), class = "symbolic_histogram")
    )
df_cardiological_hist <- data.frame(
    individuo = c(1:11)
)


df_cardiological_hist$pulso_cardiaco <- pulso_cardiaco
df_cardiological_hist$presion_sistolica <- presion_sistolica
df_cardiological_hist$presion_diastolica <- presion_diastolica


# Se procede a preparar los datos 
df_cardiological_hist <- df_cardiological_hist %>%
    mutate(individuo = as.character(individuo)) %>%
    group_by(individuo) %>%
    mutate(
        pulso_cardiaco = uniform_histogram(pulso_cardiaco), 
        presion_sistolica = median(presion_sistolica), 
        presion_diastolica = median(presion_diastolica)
    ) %>%
    ungroup()


# Se procede a probar la funcion de regresion 
prueba_reg_tree <- reg_tree(formula = pulso_cardiaco ~ presion_diastolica + presion_sistolica, data = df_cardiological_hist %>% dplyr::select(-individuo), minsize = 1)
