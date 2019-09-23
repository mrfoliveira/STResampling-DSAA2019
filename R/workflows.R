#' Summarize metrics
#'
#' Used for internal validation
#' 
#' @param int.res A list with results obtained by running estimates
#' @param metrics a list of metrics that should be summarized
#'
#' @return a named vector with median, IQR, mean, standard-deviation,
#' and number of non-NA values of each metric
summarize_metrics <- function(int.res, metrics){
  
  int.res <- int.res[, metrics, drop=F]
  
  medians <- apply(int.res, 2, stats::median, na.rm=TRUE)
  iqrs <- apply(int.res, 2, stats::IQR, na.rm=TRUE)
  means <- apply(int.res, 2, mean, na.rm=TRUE)
  sds <- apply(int.res, 2, stats::sd, na.rm=TRUE)
  nonNAs <- apply(int.res, 2, function(x) length(which(!is.na(x))))
  
  nms <- c(paste0("MED", metrics), paste0("IQR", metrics),
           paste0("MEAN", metrics), paste0("SD", metrics),
           paste0("nonNAs", metrics))
  
  res <- c(medians, iqrs, means, sds, nonNAs)
  names(res) <- nms
  res
}

#' Handling NAs in train and test
#' 
#' Discard columns/rows with too many NAs and 
#' then impute with central value.
#'
#' @param train training data set
#' @param test testing data set
#' @param nORp minimum percentage of NA values in a row/column for
#' that row/column to be discarded from training set
#'
#' @return list with an entry for the training and test sets
#' (in this order), both now with no NA values
#' 
#' @export
centralImputNAs <- function(train, test, nORp){
  # discard columns with too many NAs if there would be still predictors left
  discCols <- which(sapply(train, 
                           function(y) length(which(is.na(y)))/length(y)) > nORp)
  
  if(length(discCols) > 0) {
    if( (ncol(train) - length(discCols) - 3) > 0 )
      train <- train[, -discCols]
  }
  
  # fill in empty values in test
  if(anyNA(test)) test <- DMwR2::centralImputation(test)
  # check if any columns that are being used in train are empty in test
  if(anyNA(test)){
    emptyCols <- which(apply(test, 2, function(x) 
      length(which(is.na(x)))==length(x)))
    emptyColNms <- names(test)[emptyCols]
    
    # fill in test columns with central value of train column
    if(any(emptyColNms %in% colnames(train))){
      for(e in emptyColNms){
        if(e %in% colnames(train)) 
          test[,e] <- DMwR2::centralValue(train[,e])
      }
    }
  }
  
  # discard rows with too many NAs in train
  suppressWarnings( idxs <- DMwR2::manyNAs(train, nORp = nORp) )
  if(length(idxs)) train <- train[-idxs, ]
  # fill in empty value in train
  if(anyNA(train)) train <- DMwR2::centralImputation(train)
  
  list(train=train, test=test)
}


#' A simple learning and prediction workflow
#' 
#' A simple learning and prediction workflow that may deal
#' with NAs and use re-sampling techniques to balance an
#' imbalanced regression problem.
#' 
#' @param train a data frame for training
#' @param test a data frame for testing
#' @param time the name of the column in \code{train} and
#' \code{test} containing time-stamps
#' @param site_id the name of the column in \code{train} and
#' \code{test} containing location IDs
#' @param form a formula describing the model to learn
#' @param model the name of the algorithm to use
#' @param resample re-sampling technique to be used. Default is NULL.
#' @param resample.pars parameters to be passed to re-sample function.
#' Default is NULL.
#' @param handleNAs string indicating how to deal with NAs.
#' If "centralImput", training observations with at least 80\%
#' of non-NA columns, will have their NAs substituted by the mean
#' value and testing observatiosn will have their NAs filled in with
#' mean value regardless. Default is NULL.
#' @param min_train a minimum number of observations that must be
#' left to train a model. If there are not enough observations, 
#' predictions will be \code{NA}. Default is 2.
#' @param nORp a maximum number or fraction of columns/rows with missing
#' values above which a row/column will be removed from train before 
#' learning the model. Only works if \code{handleNAs} was
#' set to centralImputation. Default is 0.2.
#' @param ... other parameters to feed to \code{model}
#' 
#' @return a data frame containing time-stamps, location IDs,
#' true values and predicted values
#' 
#' @export
simple_workflow <- function(train, test, form, model="lm", 
                            resample = NULL, resample.pars = NULL,
                            handleNAs=NULL, 
                            min_train=2, nORp = 0.2,
                            time="time", site_id="site", ...){
  dotargs <- list(...)
  
  # get true values
  trues <- responseValues(form, test)
  
  # save original state
  orig_train <- train
  
  col.inds <- which(colnames(train) %in% c(time, site_id))
  # correct default mtry if model is ranger and there is no argument given
  if(model=="ranger" & !("mtry" %in% dotargs) & is.numeric(trues))
    dotargs$mtry <- max(floor(ncol(train[,-col.inds])/3), 1)
  
  # pre-process NAs
  if(!is.null(handleNAs)){
    if(handleNAs=="centralImput"){
      data <- centralImputNAs(train, test, nORp)
      train <- data$train
      test <- data$test
    }
  }
  
  succ <- TRUE
  if(!is.null(resample)){
    nrs <- nrow(train)
    
    if(resample %in% c("under", "over", "smote")){
      
      FUNS <- c(under="RandUnderRegress", over="RandOverRegress", smote="SmoteRegress")
      
      # ublfun <- get(FUNS[resample], asNamespace("UBL"))
      
      try( train <- do.call(FUNS[resample], #ublfun, 
                            c(list(form=form, dat=train[, -col.inds]), 
                                         resample.pars)) )
      
    }else{
      assertthat::assert_that(resample %in% c("stunder", "stover", "stsmote"))
      FUNS <- c(stunder="randUnderRegress_ST", stover="randOverRegress_ST", stsmote="SmoteRegress_ST")
      
      try( train <- do.call(FUNS[resample], c(list(form=form, dat=train, 
                                              site_id = site_id, time=time), 
                                              resample.pars)) )
    }
    
    succ <- (nrow(train)!=nrs)
  }
  
  if(nrow(train)>=min_train & succ){
    # check if columns need to be removed from train
    tr.col.inds <- which(colnames(train) %in% c(time, site_id))
    ts.col.inds <- which(colnames(test) %in% c(time, site_id))
    # removing offending column indices from train
    if(length(tr.col.inds)) train <- train[,-tr.col.inds]
    # train model
    m <- do.call(model, c(list(form, train), dotargs))
    # make predictions
    if(model=="ranger"){
      preds <- stats::predict(m, test[,-ts.col.inds])$predictions
    } else{
      preds <- stats::predict(m, test[,-ts.col.inds])
      if (is.numeric(train[[as.character(form[[2]])]]) && !is.null(dim(preds)))
        preds <- preds[, 1]
    } 
    # prepare result object
    res <- list( results = data.frame(time=test[[time]], site_id=test[[site_id]],
                      trues=trues, preds=preds) )
  }else{
    warning("nrow(train)<min_train", call. = FALSE)
    res <- list( results = data.frame(time=test[[time]], site_id=test[[site_id]],
                      trues=trues, preds=as.numeric(NA)) )
  }
  res <- c(res,
           list(orig_nrow_tr = nrow(orig_train),
                final_nrow_tr = nrow(train),
                final_trainCols = colnames(train),
                initial_tgt_sum = summary(orig_train[[as.character(form[[2]])]]),
                final_tgt_sum = summary(train[[as.character(form[[2]])]])) )
  colnames(res$results)[1:2] <- c(time, site_id)
  res
}


#' A learning and prediction workflow with internal validation
#' 
#' A learning and prediction workflow that may deal
#' with NAs and use internal validation to parametrize a
#' re-sampling technique to balance an imbalanced regression problem.
#' 
#' @inheritParams simple_workflow
#' @param internal.est character string identifying the internal estimator
#'  function to use
#' @param internal.est.pars named list of internal estimator
#' parameters (e.g., tr.perc or nfolds)
#' @param internal.evaluator character string indicating internal 
#' evaluation function
#' @param internal.eval.pars named list of parameters to feed to internal 
#' evaluation function
#' @param stat parameter indicating summary statistic that should be 
#' used to determine the best internal evaluation metric: 
#' "MED" (for median) or "MEAN" (for mean)
#' @param metrics vector of names of two metrics to be used to determine the best 
#' parametrization (the second metric is only used in case of ties)
#' @param metrics.max vector of Booleans indicating whether each metric in 
#' parameter metrics should be maximized (TRUE) or minimized (FALsE)
#' for best results
#' @param resample.grid a data.frame with columns indicating
#' resample.pars to test using internal.est
#' @param .full_intRes a Boolean indicating whether the full results
#' object for internal validation should be returned as well. Defaults to FALSE
#' 
#' @return a data frame containing time-stamps, location IDs,
#' true values and predicted values
#' 
#' @export
internal_workflow <- function(train, test, form, model="lm", 
                            resample = NULL, resample.pars = NULL,
                            internal.est = NULL, internal.est.pars = NULL,
                            internal.evaluator = NULL, internal.eval.pars = NULL,
                            metrics = NULL, metrics.max = NULL, stat="MED",
                            resample.grid = NULL, 
                            handleNAs=NULL, 
                            min_train=2, nORp = 0.2,
                            time="time", site_id="site", 
                            .full_intRes = FALSE, ...){
  dotargs <- list(...)
  
  # get true values
  trues <- responseValues(form, test)
  
  # save original state
  orig_train <- train
  
  col.inds <- which(colnames(train) %in% c(time, site_id))
  # correct default mtry if model is ranger and there is no argument given
  if(model=="ranger" & !("mtry" %in% dotargs) & is.numeric(trues))
    dotargs$mtry <- max(floor(ncol(train[,-col.inds])/3), 1)
  
  
  
  if(is.null(resample.pars)) resample.pars <- list()
  df <- data.frame(matrix(ncol = length(metrics)*5, nrow = nrow(resample.grid)))
  colnames(df) <- c(paste0("MED", metrics), paste0("IQR", metrics),
                    paste0("MEAN", metrics), paste0("SD", metrics),
                    paste0("nonNAs", metrics))
  grid.res <- cbind(resample.grid, df)
  
  # the relevance function for internal evaluation should be calculated 
  # using the whole available training values when re-sampling
  int_resample.pars <- resample.pars
  rel <- int_resample.pars$rel
  if(is.null(rel)) rel <- "auto"
  
  if(rel == "auto"){
    cf <- if(!("cf" %in% names(internal.eval.pars))) 1.5 else internal.eval.pars$cf
    method <- if(!("method" %in% names(internal.eval.pars))) "extremes" else internal.eval.pars$method
    
    int_resample.pars$rel <- uba::phi.control(train[, as.character(form[[2]])], 
                                              method=method, coef=cf)
  }
  
  if(is.list(int_resample.pars$rel)){
      
    # the relevance function for internal evaluation should be calculated 
    # using the whole available training values when calculating internal evaluation metrics
    if(!("y_train" %in% names(internal.eval.pars))){
      internal.eval.pars$y_train <- train[, as.character(form[[2]])]
    }
    
    int.results <- list()
    # do internal validation 
    for(i in 1:nrow(resample.grid)){
      
      # update resample.pars to params to test internally
      for(j in 1:ncol(resample.grid))
        int_resample.pars[[colnames(resample.grid)[j]]] <- resample.grid[i,j]
      
      # get results with these parameters
      int.res <- NULL
      try( int.res <- estimates(data = train, form = form, 
                           estimator = internal.est,
                           est.pars = internal.est.pars, 
                           workflow = "simple_workflow", 
                           wf.pars = c(list(model = model,
                                            handleNAs = handleNAs, 
                                            min_train = min_train, 
                                            nORp = nORp,
                                            resample = resample, 
                                            resample.pars = int_resample.pars), 
                                            dotargs), 
                           evaluator = internal.evaluator, 
                           eval.pars = internal.eval.pars, 
                           seed = NULL) )
      
      # if re-sampling didn't fail, save results
      if(!is.null(int.res)){
        if(.full_intRes){
          int.results[[i]] <- int.res
        }else{
          int.results[[i]] <- list(evalRes=int.res$evalRes)
        } 
        int.results[[i]]$resample.pars <- int_resample.pars
        int.results[[i]]$internal.eval.pars <- c(list(internal.evaluator=internal.evaluator),
                                                 internal.eval.pars)
        
        int.res <- int.res$evalRes
        grid.res[i, c((ncol(resample.grid) + 1):ncol(grid.res))] <- summarize_metrics(int.res, metrics)
      }
    }
  
    if(!is.null(handleNAs)){
      # pre-process NAs
      if(handleNAs=="centralImput"){
        data <- centralImputNAs(train, test, nORp)
        train <- data$train
        test <- data$test
      }
    }
    
    nrs <- nrow(train)  
    if(!(all(is.na(grid.res[, paste0(stat, metrics)])))){
      # find best results from internal validation
      if(metrics.max[1] == TRUE){
        best1 <- max(grid.res[,paste0(stat, metrics[1])], na.rm=TRUE) 
      }else{
        best1 <- min(grid.res[,paste0(stat, metrics[1])], na.rm=TRUE) 
      }
        
      best <- which(grid.res[,paste0(stat, metrics[1])] == best1)
      if(length(best) >1){
        if(metrics.max[2] == TRUE)
          best2 <- max(grid.res[best, paste0(stat, metrics[2])], na.rm=TRUE)
        else
          best2 <- min(grid.res[best, paste0(stat, metrics[2])], na.rm=TRUE)
        best <- which(grid.res[,paste0(stat, metrics[1])] == best1 & 
                        grid.res[,paste0(stat, metrics[2])] == best2)[1]
      }
      
    # update resample.pars to best results
    for(i in 1:ncol(resample.grid))
      resample.pars[[colnames(resample.grid)[i]]] <- resample.grid[best,i]
    
      if(resample %in% c("under", "over", "smote")){
        
        FUNS <- c(under="RandUnderRegress", over="RandOverRegress", smote="SmoteRegress")
        
        #ublfun <- get(FUNS[resample], asNamespace("UBL"))
        
        try( train <- do.call(FUNS[resample], #ublfun, 
                              c(list(form=form, dat=train[, -col.inds]), 
                                   resample.pars)) )
        
      }else{
        assertthat::assert_that(resample %in% c("stunder", "stover", "stsmote"))
        FUNS <- c(stunder="randUnderRegress_ST", stover="randOverRegress_ST", stsmote="SmoteRegress_ST")
        
        try( train <- do.call(FUNS[resample], c(list(form=form, dat=train, 
                                                site_id = site_id, time=time), 
                                           resample.pars)) )
      }
    }
    
    if( nrow(train)!=nrs & nrow(train)>=min_train & 
        !all(is.na(grid.res[,paste0(stat, metrics)])) ){
      # check if columns need to be removed from train
      tr.col.inds <- which(colnames(train) %in% c(time, site_id))
      ts.col.inds <- which(colnames(test) %in% c(time, site_id))
      # removing offending column indices from train
      if(length(tr.col.inds)) train <- train[,-tr.col.inds]
        # train model
        m <- do.call(model, c(list(form, train), dotargs))
        # make predictions
        if(model=="ranger"){
          preds <- stats::predict(m, test[,-ts.col.inds])$predictions
        } else{
          preds <- stats::predict(m, test[,-ts.col.inds])
          if (is.numeric(train[[as.character(form[[2]])]]) && !is.null(dim(preds)))
            preds <- preds[, 1]
        } 
        
        # prepare result object
        res <- list( results = data.frame(time=test[[time]], site_id=test[[site_id]],
                          trues=trues, preds=preds),
                     resample_grid = grid.res,
                     resample_chosen = grid.res[best, 1:ncol(resample.grid), drop=F])
        if(!is.null(int.results))
          res$internal_results <- int.results
        
    }else{
      
      warning("resampling failed OR nrow(train)<min_train", call. = FALSE)
      
      res <- list( results = data.frame(time=test[[time]], site_id=test[[site_id]],
                                        trues=trues, preds=as.numeric(NA)),
                   resample_grid = grid.res,
                   resample_chosen = grid.res[integer(0),])
    
    }
  }else{
    warning("resampling failed", call. = FALSE)
    res <- list(results = data.frame(time=test[[time]], site_id=test[[site_id]],
                      trues=trues, preds=as.numeric(NA)),
                resample_grid = grid.res,
                resample_chosen = grid.res[integer(0),])
  }
  
  res <- c(res,
           list(orig_nrow_tr = nrow(orig_train),
           final_nrow_tr = nrow(train),
           final_trainCols = colnames(train),
           initial_tgt_sum = summary(orig_train[[as.character(form[[2]])]]),
           final_tgt_sum = summary(train[[as.character(form[[2]])]])) )
  colnames(res$results)[1:2] <- c(time, site_id)
  
  res
}
