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
    print(mse)
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
install.packages("partykit")

library(partykit)

sym_tree <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...){
  y_hat <- mean(y)

  loss <- sum(mse_histogram(y, y_hat))

  out <- list(coefficients = y_hat, objfun = loss)
  class(out) <- c("sym_tree", class(out))
  out
}

logLik.sym_tree <- function(x){
  -x$objfun
}


coef.sym_tree <- function(x){
  -x$coefficients
}


class(df_OzoneFull) <- c("symbolic_tbl", class(df_OzoneFull))

model.frame <- function(drop.unused.levels, formula, data){
  # extract target
  y <- data[, as.character(formula)[2]]

  data[as.character(formula)[2]] <- 0
  
  # get the design matrix
  X <- stats::model.frame(drop.unused.levels, formula, data)

  

  X
}


mob <- function(formula, data, subset, na.action, weights, offset, cluster,
  fit, control = mob_control(), ...)
{
  ## check fitting function
  fitargs <- names(formals(fit))
  if(!all(c("y", "x", "start", "weights", "offset") %in% fitargs)) {
    stop("no suitable fitting function specified")
  }

  ## augment fitting function (if necessary)
  if(!all(c("estfun", "object") %in% fitargs)) {
    afit <- function(y,
      x = NULL, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...,
      estfun = FALSE, object = FALSE)
    {
      obj <- if("cluster" %in% fitargs) {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, cluster = cluster, ...)
      } else {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, ...)
      }
      list(
        coefficients = coef(obj),
        objfun = -as.numeric(logLik(obj)),
        estfun =  NULL,
        object =  NULL
      )
    }
  } else {
    if("cluster" %in% fitargs) {
      afit <- fit
    } else {
      afit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ..., estfun = FALSE, object = FALSE) {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, ..., estfun = estfun, object = object)
      }
    }
  }

  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula FIXME: y ~ . or y ~ x | .
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::Formula(formula(Formula::as.Formula(formula(formula), ~ 0), rhs = 2L:1L))
    xreg <- FALSE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1L:2L))
      warning("Formula must not have more than two RHS parts")
    }
    xreg <- TRUE
  }
  mf$formula <- formula

  ## evaluate model.frame
#  mf[[1L]] <- quote(model.frame)
#  mf <- eval(mf, parent.frame())

  ## extract terms, response, regressor matrix (if any), partitioning variables
  mt <- terms(formula, data = data)
  mtY <- terms(formula, data = data, rhs = if(xreg) 1L else 0L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  # Y <- switch(control$ytype,
  #   "vector" = Formula::model.part(formula, mf, lhs = 1L)[[1L]],
  #   "matrix" = model.matrix(~ 0 + ., Formula::model.part(formula, mf, lhs = 1L)),
  #   "data.frame" = Formula::model.part(formula, mf, lhs = 1L)
  # )
  # X <- if(!xreg) NULL else switch(control$xtype,
  #   "matrix" = model.matrix(mtY, mf),
  #   "data.frame" = Formula::model.part(formula, mf, rhs = 1L)
  # )
  # if(!is.null(X) && ncol(X) < 1L) {
  #   X <- NULL
  #   xreg <- FALSE
  # }
  # if(xreg) {
  #   attr(X, "formula") <- formula(formula, rhs = 1L)
  #   attr(X, "terms") <- mtY
  #   attr(X, "offset") <- cl$offset
  #   attr(X, "xlevels") <- .getXlevels(mtY, mf)
  # }
  # Z <- Formula::model.part(formula, mf, rhs = 2L)
  nyx <- length(mf) - length(Z) - as.numeric("(weights)" %in% names(mf)) - as.numeric("(offset)" %in% names(mf)) - as.numeric("(cluster)" %in% names(mf))
  varindex <- match(names(Z), names(mf))
  Y <- df_OzoneFull$Ozone.Conc.ppb
  X <- df_OzoneFull %>%
    dplyr::select(-Ozone.Conc.ppb, -individuo)

  Z <- X
  n <- nrow(Z)
  ## weights and offset
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1L
  if(length(weights) == 1L) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  offset <- if(xreg) model.offset(mf) else NULL
  cluster <- mf[["(cluster)"]]

  ## process pruning options (done here because of "n")
  if(!is.null(control$prune)) {
    if(is.character(control$prune)) {
      control$prune <- tolower(control$prune)
      control$prune <- match.arg(control$prune, c("aic", "bic", "none"))
      control$prune <- switch(control$prune,
        "aic" = {
	  function(objfun, df, nobs) (2 * objfun[1L] + 2 * df[1L]) < (2 * objfun[2L] + 2 * df[2L])
	}, "bic" = {
	  function(objfun, df, nobs) (2 * objfun[1L] + log(n) * df[1L]) < (2 * objfun[2L] + log(n) * df[2L])
	}, "none" = {
	  NULL
	})      
    }
    if(!is.function(control$prune)) {
      warning("Unknown specification of 'prune'")
      control$prune <- NULL
    }
  }
  
  
  ## grow the actual tree
  nodes <- partykit:::mob_partynode(Y = Y, X = X, Z = Z, weights = weights, offset = offset, cluster = cluster,
    fit = afit, control = control, varindex = varindex)

  ## compute terminal node number for each observation
  fitted <- fitted_node(nodes, data = mf)
  fitted <- data.frame(
      "(fitted)" = fitted,
      ## "(response)" = Y, ## probably not really needed
      check.names = FALSE,
      row.names = rownames(mf))
  if(!identical(weights, rep.int(1L, n))) fitted[["(weights)"]] <- weights
  if(!is.null(offset)) fitted[["(offset)"]] <- offset
  if(!is.null(cluster)) fitted[["(cluster)"]] <- cluster

  ## return party object
  rval <- party(nodes, 
    data = if(control$model) mf else mf[0,],
    fitted = fitted,
    terms = mt,
    info = list(
      call = cl,
      formula = oformula,
      Formula = formula,
      terms = list(response = mtY, partitioning = mtZ),
      fit = afit,
      control = control,
      dots = list(...),
      nreg = max(0L, as.integer(xreg) * (nyx - NCOL(Y))))
  )
  class(rval) <- c("modelparty", class(rval))
  return(rval)
}


model.frame(Ozone.Conc.ppb ~ Temperature.C + Solar.Radiation.WattM2 + Wind.Speed.mSec, df_OzoneFull %>% dplyr::select(-individuo))


prueba <- mob(formula = Ozone.Conc.ppb ~ 1 | Temperature.C + Solar.Radiation.WattM2 + Wind.Speed.mSec, 
  data = df_OzoneFull %>% dplyr::select(-individuo),  
  
  fit = sym_tree
  )






mob_partynode <- function(Y, X, Z, weights = NULL, offset = NULL, cluster = NULL,
  fit, control = mob_control(), varindex = 1L:NCOL(Z), ...)
{
  ## are there regressors?
  if(missing(X)) X <- NULL
  xreg <- !is.null(X)
  n <- nrow(Z)
  if(is.null(weights)) weights <- 1L
  if(length(weights) < n) weights <- rep(weights, length.out = n)

  ## control parameters (used repeatedly)
  minsize <- control$minsize
  if(!is.null(minsize) && !is.integer(minsize)) minsize <- as.integer(minsize)
  verbose <- control$verbose
  rnam <- c("estfun", "object")
  terminal <- lapply(rnam, function(x) x %in% control$terminal)
  inner    <- lapply(rnam, function(x) x %in% control$inner)
  names(terminal) <- names(inner) <- rnam

  ## convenience functions
  w2n <- function(w) if(control$caseweights) sum(w) else sum(w > 0)
  suby <- function(y, index) {
    if(control$ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  subx <- if(xreg) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels")   <- attr(x, "xlevels")
      attr(sx, "formula")   <- attr(x, "formula")
      attr(sx, "terms")     <- attr(x, "terms")
      attr(sx, "offset")    <- attr(x, "offset")
      sx
    }
  } else {
    function(x, index) NULL
  }
  subz <- function(z, index) z[index, , drop = FALSE]
  ## from strucchange
  root.matrix <- function(X) {
    if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
      X.eigen <- eigen(X, symmetric = TRUE)
      if(any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
      sqomega <- sqrt(diag(X.eigen$values))
      V <- X.eigen$vectors
      return(V %*% sqomega %*% t(V))
    }
  }

## core mob_grow_* functions

  ## variable selection: given model scores, conduct
  ## all M-fluctuation tests for orderins in z
  mob_grow_fluctests <- function(estfun, z, weights, obj = NULL, cluster = NULL){  
    ## set up return values
    m <- NCOL(z)
    pval <- rep.int(NA_real_, m)
    stat <- rep.int(0, m)
    ifac <- rep.int(FALSE, m)
    
    ## variables to test
    mtest <- if(m <= control$mtry) 1L:m else sort(sample(1L:m, control$mtry))

    ## estimating functions (dropping zero weight observations)
    process <- as.matrix(estfun)
    ww0 <- (weights > 0)
    process <- process[ww0, , drop = FALSE]
    z <- z[ww0, , drop = FALSE]
    k <- NCOL(process)
    n <- NROW(process)
    nobs <- if(control$caseweights && any(weights != 1L)) sum(weights) else n

    ## scale process
    process <- process/sqrt(nobs)
    vcov <- control$vcov
    if(is.null(obj)) vcov <- "opg"
    if(vcov != "opg") {
      bread <- vcov(obj) * nobs
    }
    if(vcov != "info") {
      ## correct scaling of estfun for variance estimate:
      ## - caseweights=FALSE: weights are integral part of the estfun -> squared in estimate
      ## - caseweights=TRUE: weights are just a factor in variance estimate -> require division by sqrt(weights)
      meat <- if(is.null(cluster)) {
        crossprod(if(control$caseweights) process/sqrt(weights) else process)
      } else {
        ## nclus <- length(unique(cluster)) ## nclus / (nclus - 1L) * 
        crossprod(as.matrix(apply(if(control$caseweights) process/sqrt(weights) else process, 2L, tapply, as.numeric(cluster), sum)))
      }
    }
    J12 <- root.matrix(switch(vcov,
      "opg" = chol2inv(chol(meat)),
      "info" = bread,
      "sandwich" = bread %*% meat %*% bread
    ))
    process <- t(J12 %*% t(process)) ## NOTE: loses column names

    ## select parameters to test
    if(!is.null(control$parm)) {
      if(is.character(control$parm)) colnames(process) <- colnames(estfun)
      process <- process[, control$parm, drop = FALSE]
    }
    k <- NCOL(process)

    ## get critical values for supLM statistic
    from <- if(control$trim > 1) control$trim else ceiling(nobs * control$trim)
    from <- max(from, minsize)
    to <- nobs - from
    lambda <- ((nobs - from) * to)/(from * (nobs - to))

    beta <- mob_beta_suplm
    logp.supLM <- function(x, k, lambda)
    {
      if(k > 40L) {
        ## use Estrella (2003) asymptotic approximation
        logp_estrella2003 <- function(x, k, lambda)
  	  -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * (1 - k/x) + 2/x))
        ## FIXME: Estrella only works well for large enough x
        ## hence require x > 1.5 * k for Estrella approximation and
        ## use an ad hoc interpolation for larger p-values
        p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
      } else {
        ## use Hansen (1997) approximation
        nb <- ncol(beta) - 1L
        tau <- if(lambda < 1) lambda else 1/(1 + sqrt(lambda))
        beta <- beta[(((k - 1) * 25 + 1):(k * 25)),]
        dummy <- beta[,(1L:nb)] %*% x^(0:(nb-1))
        dummy <- dummy * (dummy > 0)
        pp <- pchisq(dummy, beta[,(nb+1)], lower.tail = FALSE, log.p = TRUE)
        if(tau == 0.5) {
          p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
        } else if(tau <= 0.01) {
          p <- pp[25L]
        } else if(tau >= 0.49) {
          p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
	  ## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
	  if(!is.finite(p)) p <- mean(c(pp[1L], pchisq(x, k, lower.tail = FALSE, log.p = TRUE)))
        } else {
  	  taua <- (0.51 - tau) * 50
    	  tau1 <- floor(taua)
  	  p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1 + 1L]))
	  ## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
	  if(!is.finite(p)) p <- mean(pp[tau1 + 0L:1L])
        }
      }
      return(as.vector(p))
    }

    ## compute statistic and p-value for each ordering
    for(i in mtest) {
      zi <- z[,i]
      if(length(unique(zi)) < 2L) next
      if(is.factor(zi)) {
        oi <- order(zi)
        proci <- process[oi, , drop = FALSE]
        ifac[i] <- TRUE
	iord <- is.ordered(zi) & (control$ordinal != "chisq")

        ## order partitioning variable
        zi <- zi[oi]
        # re-apply factor() added to drop unused levels
        zi <- factor(zi, levels = unique(zi))
        # compute segment weights
        segweights <- if(control$caseweights) tapply(weights[oi], zi, sum) else table(zi)
	segweights <- as.vector(segweights)/nobs

        # compute statistic only if at least two levels are left
        if(length(segweights) < 2L) {
          stat[i] <- 0
	  pval[i] <- NA_real_
        } else if(iord) {
          proci <- apply(proci, 2L, cumsum)
	  tt0 <- head(cumsum(table(zi)), -1L)
          tt <- head(cumsum(segweights), -1L)
          if(control$ordinal == "max") {
  	    stat[i] <- max(abs(proci[tt0, ] / sqrt(tt * (1-tt))))
	    pval[i] <- log(as.numeric(1 - mvtnorm::pmvnorm(
	      lower = -stat[i], upper = stat[i],
	      mean = rep(0, length(tt)),
	      sigma = outer(tt, tt, function(x, y)
	        sqrt(pmin(x, y) * (1 - pmax(x, y)) / ((pmax(x, y) * (1 - pmin(x, y))))))
	      )^k))
	  } else {
	    proci <- rowSums(proci^2)
  	    stat[i] <- max(proci[tt0] / (tt * (1-tt)))
	    pval[i] <- log(strucchange::ordL2BB(segweights, nproc = k, nrep = control$nrep)$computePval(stat[i], nproc = k))
	  }
	} else {      
          stat[i] <- sum(sapply(1L:k, function(j) (tapply(proci[,j], zi, sum)^2)/segweights))
          pval[i] <- pchisq(stat[i], k*(length(levels(zi))-1), log.p = TRUE, lower.tail = FALSE)
        }
      } else {
        oi <- if(control$breakties) {
          mm <- sort(unique(zi))
	  mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 1)
  	  order(zi + runif(length(zi), min = -mm, max = +mm))
        } else {
          order(zi)
        }
        proci <- process[oi, , drop = FALSE]
        proci <- apply(proci, 2L, cumsum)
	tt0 <- if(control$caseweights && any(weights != 1L)) cumsum(weights[oi]) else 1:n
	from_to <- tt0 >= from & tt0 <= to
        stat[i] <- if(sum(from_to) > 0L) {
	  xx <- rowSums(proci^2)
	  xx <- xx[from_to]
	  tt <- tt0[from_to]/nobs
	  max(xx/(tt * (1 - tt)))	  
	} else {
	  0
	}
        pval[i] <- if(sum(from_to) > 0L) logp.supLM(stat[i], k, lambda) else NA
      }
    }

    ## select variable with minimal p-value
    best <- which.min(pval)
    if(length(best) < 1L) best <- NA
    rval <- list(pval = exp(pval), stat = stat, best = best)
    names(rval$pval) <- names(z)
    names(rval$stat) <- names(z)
    if(!all(is.na(rval$best)))
      names(rval$best) <- names(z)[rval$best]
    return(rval)
  }

  ### split in variable zselect, either ordered (numeric or ordinal) or nominal
  mob_grow_findsplit <- function(y, x, zselect, weights, offset, cluster, ...){
    ## process minsize (to minimal number of observations)
    if(minsize > 0.5 & minsize < 1) minsize <- 1 - minsize
    if(minsize < 0.5) minsize <- ceiling(w2n(weights) * minsize)
  
    if(is.numeric(zselect)) {
    ## for numerical variables
      uz <- sort(unique(zselect))
      if (length(uz) == 0L) stop("Cannot find admissible split point in partitioning variable")
      
      ## if starting values are not reused then the applyfun() is used for determining the split
      if(control$restart) {
        get_dev <- function(i) {
          zs <- zselect <= uz[i]
          if(w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
            return(Inf)
          } else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, zs), start = NULL,
	      weights = weights[zs], offset = offset[zs], cluster = cluster[zs], ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, !zs), start = NULL,
	      weights = weights[!zs], offset = offset[!zs], cluster = cluster[!zs], ...)
	    return(fit_left$objfun + fit_right$objfun)
          }
	}
        dev <- unlist(control$applyfun(1L:length(uz), get_dev))
      } else {
      ## alternatively use for() loop to go through all splits sequentially
      ## and reuse previous parameters as starting values
        dev <- vector(mode = "numeric", length = length(uz))

        start_left <- NULL
        start_right <- NULL

        for(i in 1L:length(uz)) {
          zs <- zselect <= uz[i]
  	  if(control$restart ||
	     !identical(names(start_left), names(start_right)) ||
	     !identical(length(start_left), length(start_right)))
	  {
	    start_left <- NULL
	    start_right <- NULL
 	  }
          if(w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
            dev[i] <- Inf
          } else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, zs), start = start_left,
	      weights = weights[zs], offset = offset[zs], cluster = cluster[zs], ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, !zs), start = start_right,
	      weights = weights[!zs], offset = offset[!zs], cluster = cluster[!zs], ...)
  	    start_left <- fit_left$coefficients
	    start_right <- fit_right$coefficients
	    dev[i] <- fit_left$objfun + fit_right$objfun
          }
        }
      }

      ## maybe none of the possible splits is admissible
      if(all(!is.finite(dev))) {
        split <- list(
          breaks = NULL,
	  index = NULL
        )
      } else {
        split <- list(
          breaks = if(control$numsplit == "center") {
	    as.double(mean(uz[which.min(dev) + 0L:1L]))
	  } else {
	    as.double(uz[which.min(dev)])
	  },
          index = NULL
        )
      }

    } else {

      if(!is.ordered(zselect) & control$catsplit == "multiway") {
        return(list(breaks = NULL, index = seq_along(levels(zselect))))      
      }

      ## for categorical variables      
      olevels <- levels(zselect) ## full set of original levels
      zselect <- factor(zselect) ## omit levels that do not occur in the data
      al <- mob_grow_getlevels(zselect)

      get_dev <- function(i) {
        w <- al[i,]
        zs <- zselect %in% levels(zselect)[w]
        if(w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
          return(Inf)
        } else {
	  if(nrow(al) == 1L) 1 else {
            fit_left <- fit(y = suby(y, zs), x = subx(x, zs), start = NULL,
	      weights = weights[zs], offset = offset[zs], cluster = cluster[zs], ...)
            fit_right <- fit(y = suby(y, !zs), x = subx(x, !zs), start = NULL,
	      weights = weights[!zs], offset = offset[!zs], cluster = cluster[zs], ...)
    	    fit_left$objfun + fit_right$objfun
	  }
        }
      }
      dev <- unlist(control$applyfun(1L:nrow(al), get_dev))

      if(all(!is.finite(dev))) {
        split <- list(
          breaks = NULL,
	  index = NULL
        )
      } else {
        if(is.ordered(zselect)) {
	  ## map back to set of full original levels
          split <- list(
            breaks = match(levels(zselect)[which.min(dev)], olevels),
	    index = NULL
          )
        } else {
	  ## map back to set of full original levels
	  ix <- structure(rep.int(NA_integer_, length(olevels)), .Names = olevels)
	  ix[colnames(al)] <- !al[which.min(dev),]
	  ix <- as.integer(ix) + 1L
          split <- list(
            breaks = NULL,
	    index = ix
          )
        }
      }
    }
  
    return(split)
  }

  ## grow tree by combining fluctuation tests for variable selection
  ## and split selection recursively
  mob_grow <- function(id = 1L, y, x, z, weights, offset, cluster, ...){
    if(verbose) {
      if(id == 1L) cat("\n")
      cat(sprintf("-- Node %i %s\n", id, paste(rep("-", 32 - floor(log10(id)) + 1L), collapse = "")))
      cat(sprintf("Number of observations: %s\n", w2n(weights)))
      ## cat(sprintf("Depth: %i\n", depth))
    }

    ## fit model
    mod <- fit(y, x, weights = weights, offset = offset, cluster = cluster, ...,
      estfun = TRUE, object = terminal$object | control$vcov == "info")
    mod$test <- NULL
    mod$nobs <- w2n(weights)
    mod$p.value <- NULL

    ## set default for minsize if not specified
    if(is.null(minsize)) minsize <<- as.integer(ceiling(10L * length(mod$coefficients)/NCOL(y)))

    ## if too few observations or maximum depth: no split = return terminal node
    TERMINAL <- FALSE
    if(w2n(weights) < 2 * minsize) {
      if(verbose) cat(sprintf("Too few observations, stop splitting (minsize = %s)\n\n", minsize))
      TERMINAL <- TRUE
    }
    if(depth >= control$maxdepth) {
      if(verbose) cat(sprintf("Maximum depth reached, stop splitting (maxdepth = %s)\n\n", control$maxdepth))
      TERMINAL <- TRUE
    }
    if(TERMINAL) {
      return(partynode(id = id, info = mod))
    }

    ## conduct all parameter instability tests
    test <- if(is.null(mod$estfun)) NULL else try(mob_grow_fluctests(mod$estfun, z, weights, mod$object, cluster))

    if(!is.null(test) && !inherits(test, "try-error")) {
      if(control$bonferroni) {
        pval1 <- pmin(1, sum(!is.na(test$pval)) * test$pval)
        pval2 <- 1 - (1 - test$pval)^sum(!is.na(test$pval))
        test$pval <- ifelse(!is.na(test$pval) & (test$pval > 0.001), pval2, pval1)
      }

      best <- test$best
      TERMINAL <- is.na(best) || test$pval[best] > control$alpha
      mod$p.value <- as.numeric(test$pval[best])

      if (verbose) {
        cat("\nParameter instability tests:\n")
        print(rbind(statistic = test$stat, p.value = test$pval))
        cat(sprintf("\nBest splitting variable: %s", names(test$stat)[best]))
        cat(sprintf("\nPerform split? %s", ifelse(TERMINAL, "no\n\n", "yes\n")))
      }
    } else {
      if(verbose && inherits(test, "try-error")) cat("Parameter instability tests failed\n\n")
      TERMINAL <- TRUE
      test <- list(stat = NA, pval = NA)
    }
    
    ## update model information
    mod$test <- rbind("statistic" = test$stat, "p.value" = test$pval)

    if(TERMINAL) {
      return(partynode(id = id, info = mod))
    } else {
      zselect <- z[[best]]
      sp <- mob_grow_findsplit(y, x, zselect, weights, offset, cluster, ...)
    
      ## split successful?
      if(is.null(sp$breaks) & is.null(sp$index)) {
        if(verbose) cat(sprintf("No admissable split found in %s\n\n", sQuote(names(test$stat)[best])))
        return(partynode(id = id, info = mod))
      } else {
        sp <- partysplit(as.integer(best), breaks = sp$breaks, index = sp$index)
        if(verbose) cat(sprintf("Selected split: %s\n\n",
	  paste(character_split(sp, data = z)$levels, collapse = " | ")))
      }
    }
  
    ## actually split the data
    kidids <- kidids_split(sp, data = z)
    
    ## set-up all daugther nodes
    depth <<- depth + 1L
    kids <- vector(mode = "list", length = max(kidids))
    for(kidid in 1L:max(kidids)) {
      ## select obs for current next node
      nxt <- kidids == kidid

      ## get next node id
      if(kidid > 1L) {
        myid <- max(nodeids(kids[[kidid - 1L]]))
      } else {
        myid <- id
      }

      ## start recursion on this daugther node
      kids[[kidid]] <- mob_grow(id = myid + 1L,
        suby(y, nxt), subx(x, nxt), subz(z, nxt), weights[nxt], offset[nxt], cluster[nxt], ...)
    }
    depth <<- depth - 1L

    ## shift split varid from z to mf
    sp$varid <- as.integer(varindex[sp$varid])
    
    ## return nodes
    return(partynode(id = id, split = sp, kids = kids, info = mod))
  }

  mob_prune <- function(node){
    ## turn node to list
    nd <- as.list(node)

    ## if no pruning selected
    if(is.null(control$prune)) return(nd)

    ## node information (IDs, kids, ...)
    id <- seq_along(nd)
    kids <- lapply(nd, "[[", "kids")
    tmnl <- sapply(kids, is.null)

    ## check nodes that only have terminal kids
    check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    while(any(check)) {

      ## pruning node information
      pnode <- which(check)
      objfun <- sapply(nd, function(x) x$info$objfun)
      pok <- sapply(pnode, function(i) control$prune(
        objfun = c(objfun[i], sum(objfun[kids[[i]]])),
	df = c(length(nd[[1]]$info$coefficients), length(kids[[i]]) * length(nd[[1]]$info$coefficients) + as.integer(control$dfsplit)),
        nobs = c(nd[[i]]$info$nobs, n)
      ))

      ## do any nodes need pruning?
      pnode <- pnode[pok]
      if(length(pnode) < 1L) break

      ## prune nodes and relabel IDs  
      pkids <- sort(unlist(sapply(pnode, function(i) nd[[i]]$kids)))
      for(i in id) {
        nd[[i]]$kids <- if(nd[[i]]$id %in% pnode || is.null(kids[[i]])) {
          NULL
        } else {    
          nd[[i]]$kids - sapply(kids[[i]], function(x) sum(pkids < x))
        }
      }
      nd[pkids] <- NULL
      id <- seq_along(nd)
      for(i in id) nd[[i]]$id <- i
      
      ## node information
      kids <- lapply(nd, "[[", "kids")
      tmnl <- sapply(kids, is.null)
      check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    }
   
    ## return pruned list 
    return(nd)
  }

  ## grow tree
  depth <- 1L
  nodes <- mob_grow(id = 1L, Y, X, Z, weights, offset, cluster)

  ## prune tree
  if(verbose && !is.null(control$prune)) cat("-- Post-pruning ---------------------------\n")
  nodes <- mob_prune(nodes)  
  for(i in seq_along(nodes)) {
    if(is.null(nodes[[i]]$kids)) {
      nodes[[i]]$split <- NULL  
      if(!terminal$estfun) nodes[[i]]$info$estfun <- NULL
      if(!terminal$object) nodes[[i]]$info$object <- NULL      
    } else {
      if(!inner$estfun) nodes[[i]]$info$estfun <- NULL
      if(!inner$object) nodes[[i]]$info$object <- NULL      
    }
  }
  
  ## return as partynode
  as.partynode(nodes)
}