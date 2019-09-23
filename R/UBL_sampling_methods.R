#' Random under-sampling for imbalanced regression problems
#' 
#' This function is identical to the function of the same name available on the 
#' R package \code{UBL}. \cr
#' The function performs a random under-sampling strategy for
#' imbalanced regression problems. Essentially, a percentage of 
#' cases of the "class(es)" (bumps below a relevance threshold defined) 
#' selected by the user are randomly removed. Alternatively, the strategy 
#' can be applied to either balance all the existing "classes"" or to 
#' "smoothly invert" the frequency of the examples in each "class".
#'
#' @inherit UBL::RandUnderRegress return
#' @inheritParams UBL::RandUnderRegress
#' @param C.perc A vector containing the under-sampling percentage/s to apply 
#' to all/each "class" (bump) obtained with the relevance threshold. 
#' Examples are randomly removed from the "class(es)". 
#' If only one percentage is provided this value is reused in all the "classes" 
#' that have values below the relevance threshold. A different percentage can be 
#' provided to each "class". In this case, the percentages should be provided in 
#' ascending order of target variable value. The under-sampling percentage(s), 
#' should be a number below 1, meaning that the normal cases (cases below the threshold) 
#' are under-sampled by the corresponding percentage. If the number 1 is provided then 
#' those examples are not changed. Alternatively, C.perc parameter may be set to "balance"
#'  or "extreme", cases where the under-sampling percentages are automatically estimated
#'  to either balance or invert the frequencies of the examples in the "classes" (bumps).
#' 
#' @details The only difference between this function and the original function is in the requirements 
#' imposed on the argument C.perc. \cr This function performs a random under-sampling strategy for dealing with 
#' imbalanced regression problems. The examples removed are randomly selected among the
#'  examples belonging to the normal "class(es)" (bump of relevance below the threshold defined). 
#'  The user can chose one or more bumps to be under-sampled.
#' 
#' @references Paula Branco, Rita P. Ribeiro, Luis Torgo (2016)., 
#' UBL: an R Package for Utility-Based Learning, 
#' CoRR abs/1604.08079 [cs.MS], URL: http://arxiv.org/abs/1604.08079
#' 
#' @seealso \code{\link[UBL]{RandUnderRegress}}, \code{\link{RandOverRegress}}
#' 
#' @export
RandUnderRegress <- function (form, dat, rel = "auto", thr.rel = 0.5, C.perc = "balance", 
                              repl = FALSE) 
{
  if ( (is.list(C.perc) & any(unlist(C.perc) > 1)) || 
       (is.vector(C.perc) && is.numeric(C.perc) && any(C.perc > 1) ) ){
    stop("The under-sampling percentages provided in parameter C.perc\n         can not be higher than 1!")
  }
  tgt <- which(names(dat) == as.character(form[[2]]))
  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  if (is.matrix(rel)) {
    pc <- UBL::phi.control(y, method = "range", control.pts = rel)
  }
  else if (is.list(rel)) {
    pc <- rel
  }
  else if (rel == "auto") {
    pc <- UBL::phi.control(y, method = "extremes")
  }
  else {
    stop("future work!")
  }
  temp <- y.relev <- UBL::phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. \n         Please, redefine your relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. \n         Please, redefine your relevance function!")
  }
  bumps <- c()
  for (i in 1:(length(y) - 1)) {
    if ((temp[i] >= thr.rel && temp[i + 1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i + 1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  nbump <- length(bumps) + 1
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
    newdata <- rbind(newdata, dat[names(obs.ind[[ove[j]]]), 
                                  ])
  }
  if (is.list(C.perc) || (is.vector(C.perc) && is.numeric(C.perc)) ) {
    if (length(und) > 1 & length(C.perc) == 1) {
      C.perc <- rep(C.perc[1], length(und))
    }
    else if (length(und) > length(C.perc) & length(C.perc) > 
             1) {
      stop("The number of under-sampling percentages must be equal \n           to the number of bumps below the threshold defined!")
    }
    else if (length(und) < length(C.perc)) {
      stop("the number of under-sampling percentages must be at most\n           the number of bumps below the threshold defined!")
    }
  }
  else if (C.perc == "balance") {
    B <- sum(sapply(obs.ind[ove], length))
    obj <- B/length(und)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 
                            5))
  }
  else if (C.perc == "extreme") {
    Bove <- sum(sapply(obs.ind[ove], length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und], length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 
                            5))
  }
  for (j in 1:length(und)) {
    sel <- sample(names(obs.ind[[und[j]]]), C.perc[[j]] * 
                    length(obs.ind[[und[j]]]), replace = repl)
    newdata <- rbind(newdata, dat[sel, ])
  }
  newdata
}

#' Random over-sampling for imbalanced regression problems
#' 
#' This function is identical to the function of the same name available on the 
#' R package \code{UBL}. The only difference is in the requirements
#' imposed on the argument C.perc. \cr
#' This function performs a random over-sampling strategy for imbalanced 
#' regression problems. Basically a percentage of cases of the "class(es)" 
#' (bumps above a relevance threshold defined) selected by the user are randomly 
#' over-sampled. Alternatively, it can either balance all the existing "classes" 
#' (the default) or it can "smoothly invert" the frequency of the examples in each class.
#'
#' @inherit UBL::RandOverRegress return
#' @inheritParams UBL::RandOverRegress
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
#' 
#' @details The only difference between this function and the original function is in the requirements 
#' imposed on the argument C.perc. \cr This function performs a random over-sampling strategy for dealing 
#' with imbalanced regression problems. The new examples included in the new data set
#'  are randomly selected replicas of the examples already present in the original data set.
#' 
#' @references Paula Branco, Rita P. Ribeiro, Luis Torgo (2016)., 
#' UBL: an R Package for Utility-Based Learning, 
#' CoRR abs/1604.08079 [cs.MS], URL: http://arxiv.org/abs/1604.08079
#' 
#' @seealso \code{\link[UBL]{RandOverRegress}}, \code{\link{RandUnderRegress}}
#' 
#' @export
RandOverRegress <- function (form, dat, rel = "auto", thr.rel = 0.5, C.perc = "balance", 
                             repl = TRUE) 
{
  #if (is.list(C.perc) & any(unlist(C.perc) < 1)) {
  #  stop("The over-sampling percentages provided in parameter C.perc\n         can not be lower than 1!")
  #}
  tgt <- which(names(dat) == as.character(form[[2]]))
  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  if (is.matrix(rel)) {
    pc <- UBL::phi.control(y, method = "range", control.pts = rel)
  }
  else if (is.list(rel)) {
    pc <- rel
  }
  else if (rel == "auto") {
    pc <- UBL::phi.control(y, method = "extremes")
  }
  else {
    stop("future work!")
  }
  temp <- y.relev <- UBL::phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. Please, redefine your relevance\n         function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. Please, redefine your relevance\n         function!")
  }
  bumps <- c()
  for (i in 1:(length(y) - 1)) {
    if ((temp[i] >= thr.rel && temp[i + 1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i + 1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  nbump <- length(bumps) + 1
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
  if (is.list(C.perc) || (is.vector(C.perc) && is.numeric(C.perc)) ) {
    if (length(ove) > 1 & length(C.perc) == 1) {
      C.perc <- rep(C.perc[1], length(ove))
    }
    else if (length(ove) > length(C.perc) & length(C.perc) > 
             1) {
      stop("The number of over-sampling percentages must be equal to the\n           number of bumps above the threshold defined!")
    }
    else if (length(ove) < length(C.perc)) {
      stop("The number of over-sampling percentages must be at most the \n           number of bumps above the threshold defined!")
    }
  }
  else if (C.perc == "balance") {
    B <- sum(sapply(obs.ind[und], length))
    obj <- B/length(ove)
    C.perc <- as.list(round(obj/sapply(obs.ind[ove], length), 
                            5))
  }
  else if (C.perc == "extreme") {
    Bund <- sum(sapply(obs.ind[und], length))/length(und)
    obj <- Bund^2/sapply(obs.ind[ove], length)
    C.perc <- as.list(round(obj/sapply(obs.ind[ove], length), 
                            5))
  }
  newdata <- dat
  for (j in 1:length(ove)) {
    sel <- sample(names(obs.ind[[ove[j]]]), C.perc[[j]] * 
                    length(obs.ind[[ove[j]]]), replace = repl)
    newdata <- rbind(newdata, dat[sel, ])
  }
  newdata
}
