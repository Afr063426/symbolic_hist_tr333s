#' author: Joshua Cervantes Artavia
#' date: 13.10.2024
#' description: This code have as utility to define a Mallows metric for rpart without the use
#' of C code, it is slower in comparation to C code
#' 
#' 
library(RSDA)
library(tidyverse)
library(ggalt)


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
        
        ranges <- lapply(breaks, function(x) na.omit(x-lag(x)))
        
        breaks_mean <- pred$breaks

        centers_mean <- (breaks_mean + lag(breaks_mean))/2
        ranges_mean <- (breaks_mean - lag(breaks_mean))
        centers_mean <- centers_mean[!is.na(centers_mean)]
        ranges_mean <- ranges_mean[!is.na(ranges_mean)]

        props <- pred$props
        
        sum_centers <- sapply(centers, function(x) as.numeric(props[[1]][1])*sum((as.numeric(x) - as.numeric(centers_mean))^2))
        sum_ranges <- 1/3*sapply(ranges, function(x) as.numeric(props[[1]][1])*sum((as.numeric(x) - as.numeric(ranges_mean))^2))
        return(sum_centers + sum_ranges)
    }else{
        return(0)
    }
    
}

mse_histogram_vec <- function(y, pred){
    breaks <- y$breaks
    centers <- (breaks + lag(breaks))/2
    
    ranges <- (breaks - lag(breaks))

    centers <- centers[!is.na(centers)]
    ranges <- ranges[!is.na(ranges)]      


    breaks_mean <- pred$breaks

    centers_mean <- (breaks_mean + lag(breaks_mean))/2
    ranges_mean <- (breaks_mean - lag(breaks_mean))
    centers_mean <- centers_mean[!is.na(centers_mean)]
    ranges_mean <- ranges_mean[!is.na(ranges_mean)]

    props <- y$props[1]
    
    sum_centers <- props * sum((centers_mean - centers)^2)
    sum_ranges <- 1/3 * props * sum((ranges_mean - ranges)^2)
    return(sum_centers + sum_ranges)

}
mse_histogram_vec <- Vectorize(mse_histogram_vec)



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

            mse_aux <- sum((sum(left_mse) + sum(right_mse)))/length(x)
            #print(mse_aux)
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
    #print(mse)
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
hist_tree <- function(formula, data, minsize) {
  
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
      tmp_filter <- c(paste(gsub(".sse", "",names(tmp_splitter)), ">", 
                            splitting[[tmp_splitter]][2]),
                      paste(gsub(".sse", "",names(tmp_splitter)), "<=", 
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





cv_tree <- function(formula, data, minsize = 2, size_split = 0.8, cv = 10){
  r_squared <- 0
  omega <- 0
  rmse <- 0
  variable <- as.character(formula)[2]
  sym_variable <- sym(variable)
  promedio <-  data[, as.character(formula)[2]] %>% summarise(promedio = mean(!!sym_variable))
  promedio <- promedio$promedio
  promedio_billard <- (promedio$breaks + lag(promedio$breaks))/2
  promedio_billard <- sum(promedio$props[[1]][1]*promedio_billard[!is.na(promedio_billard)])
  promedio_billard <- list(breaks = rep(promedio_billard, length(promedio$breaks)), props = rep(promedio$props[1], length(promedio$props)))
  folds <- caret::createFolds(1:nrow(data), k = 10)
  for(i in 1:cv){

    muestra <- folds[[i]]

    df_train <- data[-muestra, ]

    df_test <- data[muestra, ]

    tree <- hist_tree(formula =  formula, data = data %>% dplyr::select(-individuo),  minsize = minsize)

    df_tree <- tree$tree

    for(j in 2:nrow(df_tree)){
      
      if(j == 2){
        df <- subset(df_train, eval(parse(text = df_tree$FILTER[j])))
        df$node <- j

        df_t <- subset(df_test, eval(parse(text = df_tree$FILTER[j])))
        df_t$node <- j
      }else{
        df_aux <- subset(df_train, eval(parse(text = df_tree$FILTER[j])))
        df_aux$node <- j
        df_t_aux <- subset(df_test, eval(parse(text = df_tree$FILTER[j])))
        df_t_aux$node <- j
        df <- rbind(df, df_aux)
        df_t <- rbind(df, df_t_aux)
      }
    }

    df <- df %>%
      arrange(desc(node)) %>%
      group_by(individuo) %>%
      slice(1) %>%
      ungroup()
    
    df_estimador <- df %>%
      group_by(node) %>%
      summarise(estimador =  mean(!!sym_variable))
    
    df_t <- df_t %>%
      arrange(desc(node)) %>%
      group_by(individuo) %>%
      slice(1) %>%
      ungroup() %>%
      left_join(df_estimador, by = "node")

    y_test <- df_t[, as.character(formula)[2]]
    colnames(y_test) <- "y_test"
    y_pred <- df_t[, "estimador"]
    
    r_squared <- r_squared + 1 - mean(mse_histogram_vec(pred = y_pred$estimador, y = y_test$y_test))/mean(mse_histogram(y_test$y_test, promedio))
    rmse <- rmse + sqrt(mean(mse_histogram_vec(pred = y_pred$estimador, y = y_test$y_test)))

    
    omega <- omega + 1 - mean(mse_histogram_vec(pred = y_pred$estimador, y = y_test$y_test))/mean(mse_histogram(y_test$y_test, promedio_billard))
  }

  r_squared <- r_squared/cv
  rmse <- rmse/cv
  omega <- omega/cv

  data.frame("minsize" = minsize, "r_squared" = r_squared, "rmse" = rmse, "omega" = omega)
}
