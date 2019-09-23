#' Biased under-sampling for imbalanced regression spatio-temporal problems
#' 
#' Based on randUnderRegress (\code{R} package \code{UBL}). 
#' The function performs a random under-sampling strategy for
#' imbalanced regression problems with a bias based on spatio-temporal
#' contextual information. Essentially, a percentage of 
#' cases of the "class(es)" (bumps below a relevance threshold defined) 
#' selected by the user are randomly removed with a sampling bias
#' based on a spatio-temporal weight. Alternatively, the strategy 
#' can be applied to either balance all the existing "classes"" or to 
#' "smoothly invert" the frequency of the examples in each "class".
#'
#' @param form a model formula
#' @param dat the original training set (with the unbalanced distribution)
#' @param rel relevance determined automatically (default) with uba package 
#' or provided by the user
#' @param thr.rel relevance threshold above which a case is considered as 
#' belonging to the rare "class"
#' @param C.perc A vector containing the over-sampling percentage/s to apply to all/each 
#' "class" (bump) obtained with the relevance threshold. Replicas of the examples are are
#'  randomly added in each "class". If only one percentage is provided this value is reused 
#'  in all the "classes" that have values above the relevance threshold. A different percentage 
#'  can be provided to each "class". In this case, the percentages should be provided in 
#'  ascending order of target variable value. The over-sampling percentage(s), should be 
#'  numbers above 0, meaning that the important cases (cases above the threshold) are over-sampled 
#'  by the corresponding percentage. If the number 1 is provided then the number of extreme examples 
#'  will be doubled. 
#'  Alternatively, C.perc parameter may be set to "balance" or "extreme", cases where the 
#'  over-sampling percentages are automatically estimated to either balance or invert the
#'  frequencies of the examples in the "classes" (bumps).
#' @param repl allowed to perform sampling with replacement
#' @param type character string indicating the type of bias used. Default is "add".
#' More types to be added in future work
#' @inheritParams sample_wts
#'
#' @return The function returns a data frame with the new data set resulting 
#' from the application of the spatio-temporally biased under-sampling strategy.
#' 
#' @references Paula Branco, Rita P. Ribeiro, Luis Torgo (2016)., 
#' UBL: an R Package for Utility-Based Learning, 
#' CoRR abs/1604.08079 [cs.MS], URL: http://arxiv.org/abs/1604.08079
#' 
#' @seealso \code{\link[UBL]{RandUnderRegress}}, \code{\link{sample_wts}}
#' 
#' @export
randUnderRegress_ST <- function(form, dat, alpha=0.5, beta=0.9,
                                rel="auto", thr.rel=0.5, epsilon=1E-4,
  C.perc="balance", repl=FALSE, 
  type = "add",
  site_id="site_id", time="time", 
  sites_sf = NULL, 
  lon=NULL, lat=NULL, crs = NULL) {
  
  assertthat::assert_that(type %in% c("add", "mult", "addBetaPhi", "addPhiMult"), msg = "Check ST bias type")

  if( (is.list(C.perc) && any(unlist(C.perc)>1)) || (is.vector(C.perc) && is.numeric(C.perc) && any(C.perc>1)) ){
    stop("The under-sampling percentages provided in parameter C.perc
      can not be higher than 1!")
  }
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  
  if (is.matrix(rel)) { 
    pc <- UBL::phi.control(y, method = "range", control.pts = rel)
  } else if (is.list(rel)) { 
    pc <- rel
  } else if (rel == "auto") {
    pc <- UBL::phi.control(y, method = "extremes")
  } else {# handle other relevance functions and not using the threshold!
    stop("future work!")
  }
  
  temp <- y.relev <- UBL::phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. 
      Please, redefine your relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. 
      Please, redefine your relevance function!")
  }
  
  stprobs <- sample_wts(form = form, df = dat, phi.control = pc, 
                        alpha = alpha, beta=beta, rel.thr=thr.rel,
    site_id=site_id, time=time, sites_sf = sites_sf, 
    lon=lon, lat=lat, crs = crs, epsilon=epsilon)
  stprob <- stprobs[,paste0("stprob_", type), drop=TRUE]
  names(stprob) <- rownames(dat)
  
  
  #  temp[which(y.relev >= thr.rel)] <- -temp[which(y.relev >= thr.rel)]
  bumps <- c()
  for (i in 1:(length(y) - 1)) {
    #     if (temp[i] * temp[i + 1] < 0) {
    #       bumps <- c(bumps, i)
    #     }
    if ((temp[i] >= thr.rel && temp[i+1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i+1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  nbump <- length(bumps) + 1 # number of different "classes"
  
  # collect the indexes in each "class"
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for (i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] + 1
  }
  obs.ind[[nbump]] <- s.y[last:length(s.y)]
  
  
  imp <- sapply(obs.ind, function(x) mean(UBL::phi(x, pc)))
  
  und <- which(imp < thr.rel)
  ove <- which(imp >= thr.rel)
  
  newdata <- NULL
  for (j in 1:length(ove)) {
    newdata <- rbind(newdata, dat[names(obs.ind[[ove[j]]]), ])
    # start with the examples from the minority "classes"
  }
  
  # set the under-sampling percentages
  if (is.list(C.perc) || (is.vector(C.perc) && is.numeric(C.perc))) {
    if (length(und) > 1 & length(C.perc) == 1) { 
      # the same under-sampling percentage is applied to all the "classes"
      C.perc <- rep(C.perc[1], length(und))
    } else if (length(und) > length(C.perc) & length(C.perc) > 1) {
      stop("The number of under-sampling percentages must be equal 
        to the number of bumps below the threshold defined!")      
    }else if (length(und) < length(C.perc)) {
      stop("the number of under-sampling percentages must be at most
        the number of bumps below the threshold defined!")
    }
    } else if (C.perc == "balance") {
      B <- sum(sapply(obs.ind[ove], length))
      obj <- B/length(und)
      C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 5))
  } else if (C.perc == "extreme") {
    Bove <- sum(sapply(obs.ind[ove], length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und], length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 5))
  }
  
  for (j in 1:length(und)) {
    sel <- sample(names(obs.ind[[und[j]]]),
      C.perc[[j]] * length(obs.ind[[und[j]]]), replace = repl,
      prob = stprob[match(names(obs.ind[[und[j]]]),names(stprob))])
    newdata <- rbind(newdata, dat[sel, ])
  }
  
  newdata
  
}


#' Biased over-sampling for imbalanced regression spatio-temporal problems
#' 
#' Based on randOverRegress (\code{R} package \code{UBL}).
#' This function performs a random over-sampling strategy for imbalanced 
#' regression problems with a bias based on spatio-temporal
#' contextual information. Basically a percentage of cases of the "class(es)" 
#' (bumps above a relevance threshold defined) selected by the user are randomly 
#' over-sampled with a sampling bias based on a spatio-temporal weight.
#' Alternatively, it can either balance all the existing "classes" 
#' (the default) or it can "smoothly invert" the frequency of the examples in each class.
#'
#' @inheritParams randUnderRegress_ST
#' @inherit randUnderRegress_ST return
#' 
#' @references Paula Branco, Rita P. Ribeiro, Luis Torgo (2016)., 
#' UBL: an R Package for Utility-Based Learning, 
#' CoRR abs/1604.08079 [cs.MS], URL: http://arxiv.org/abs/1604.08079
#' 
#' @seealso \code{\link[UBL]{RandOverRegress}}, \code{\link{sample_wts}}
#' 
#' @return The function returns a data frame with the new data set resulting 
#' from the application of the spatio-temporally biased over-sampling strategy.
#' 
#' @export
randOverRegress_ST <- function(form, dat, alpha = 0.5, beta=0.9, 
                               rel="auto", thr.rel=0.5, epsilon=1E-4,
  C.perc="balance", repl=TRUE, 
  type = "add",
  site_id="site_id", time="time", 
  sites_sf = NULL, 
  lon=NULL, lat=NULL, crs = NULL) {
  
assertthat::assert_that(type %in% c("add", "mult", "addBetaPhi", "addPhiMult"), msg = "Check ST bias type")
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)
  
  s.y <- sort(y)
  if (is.matrix(rel)) { 
    pc <- UBL::phi.control(y, method = "range", control.pts = rel)
  } else if (is.list(rel)) { 
    pc <- rel
  } else if (rel == "auto") {
    pc <- UBL::phi.control(y, method = "extremes")
  } else {# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }
  
  temp <- y.relev <- UBL::phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. Please, redefine your relevance
      function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. Please, redefine your relevance
      function!")
  }
  
  stprobs <- sample_wts(form = form, df = dat, phi.control = pc,
                        alpha = alpha, beta = beta, rel.thr=thr.rel,
    site_id=site_id, time=time, sites_sf = sites_sf, 
    lon=lon, lat=lat, crs = crs, epsilon = epsilon)
  stprob <- stprobs[,paste0("stprob_", type), drop=TRUE]
  names(stprob) <- rownames(dat)
  
  #  temp[which(y.relev >= thr.rel)] <- -temp[which(y.relev >= thr.rel)]
  bumps <- c()
  for (i in 1:(length(y) - 1)) {
    #     if (temp[i] * temp[i + 1] < 0) {
    #       bumps <- c(bumps, i)
    #     }
    if ((temp[i] >= thr.rel && temp[i+1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i+1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  nbump <- length(bumps) + 1 # number of different "classes"
  
  # collect the indexes in each "class"
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for (i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] + 1
  }
  obs.ind[[nbump]] <- s.y[last:length(s.y)]
  
  imp <- sapply(obs.ind, function(x) mean(UBL::phi(x, pc)))
  
  ove <- which(imp >= thr.rel)
  und <- which(imp < thr.rel)
  
  # set the over-sampling percentages
  if (is.list(C.perc) || (is.vector(C.perc) && is.numeric(C.perc))) {
    if (length(ove) > 1 & length(C.perc) == 1) {
      # only one percentage to apply to all the "classes" 
      C.perc <- rep(C.perc[1], length(ove))
    } else if (length(ove) > length(C.perc) & length(C.perc) > 1) {
      stop("The number of over-sampling percentages must be equal to the
        number of bumps above the threshold defined!")      
    } else if (length(ove) < length(C.perc)) {
      stop("The number of over-sampling percentages must be at most the 
        number of bumps above the threshold defined!")
    }
  } else if (C.perc == "balance") {
      B <- sum(sapply(obs.ind[und], length))
      obj <- B/length(ove)
      C.perc <- as.list(round(obj/sapply(obs.ind[ove], length), 5))
  } else if (C.perc == "extreme") {
    Bund <- sum(sapply(obs.ind[und], length))/length(und)
    obj <- Bund^2/sapply(obs.ind[ove], length)
    C.perc <- as.list(round(obj/sapply(obs.ind[ove], length), 5))
  }
  
  
  newdata <- dat   
  
  for (j in 1:length(ove)) {
    sel <- sample(names(obs.ind[[ove[j]]]),
      C.perc[[j]] * length(obs.ind[[ove[j]]]), 
      replace = repl,
      prob = stprob[match(names(obs.ind[[ove[j]]]),names(stprob))])
    newdata <- rbind(newdata, dat[sel, ])
  }  
  
  newdata
  
  }

