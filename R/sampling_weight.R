#' Calculate utility-based relevance
#' 
#' Calculate relevance of values given a parametrization
#' of the relevance function. 
#' Most relevant: phi -> 1; less relevant: phi -> 0.
#' @param y vector of values to calculate relevance of
#' @param phi.control list of parameters as returned
#' by function \code{UBL::phi.control}
#' @seealso \code{\link[UBL]{phi}}, \code{\link[UBL]{phi.control}}
get_phi <- function(y, phi.control){
  # require(UBL)
  UBL::phi(y, phi.control)
}

#' Calculate temporally-biased re-sampling weights
#'
#' Calculate weights for re-sampling with a temporal bias.
#' Most recent observations have weights that tend to 1,
#' while the oldest observations have weights that tend to 0
#' (meaning they are less likely to be kept).
#' Most recent observations: w -> 1; oldest: w -> 0.
#' 
#' @param times a vector of time-stamps
#' @param phi a vector of the relevance values of 
#' \code{df}'s target variable
#' @param rel.thr a relevance threshold above which an 
#' observation is considered relevant
#' 
#' @return A vector of temporally-biased re-sampling weights, scaled
#' to fit within range [0,1].
#' @author Mariana Oliveira
get_time_wts <- function(times, phi, rel.thr){
  # check types
  assertthat::assert_that(lubridate::is.Date(times) | lubridate::is.POSIXct(times) | 
                            lubridate::is.POSIXlt(times) | lubridate::is.POSIXt(times), 
                          msg = "times must be of type Date or POSIX")
  
  # overall normal and relevant inds
  norm_inds <- which(phi < rel.thr)
  relev_inds <- which(phi >= rel.thr)
  
  time_wts <- vector(mode="numeric", length=length(times))
  time_wts <- rep(NA, length(times))
  
  # scale time so most recent = 1
  time_wts[norm_inds] <- as.numeric( lubridate::seconds( lubridate::interval(times[norm_inds], min(times[norm_inds])))) / 
    as.numeric( lubridate::seconds( lubridate::interval(max(times[norm_inds]), min(times[norm_inds])) ))

  # scale time so most recent = 1
  time_wts[relev_inds] <- as.numeric( lubridate::seconds( lubridate::interval(times[relev_inds], min(times[relev_inds])))) / 
    as.numeric( lubridate::seconds( lubridate::interval(max(times[relev_inds]), min(times[relev_inds])) ))
  
  time_wts
}

#' Calculate spatially-biased re-sampling weights
#'
#' Calculate weights for re-sampling with a spatial bias.
#' Observations have a distance that tends to 1 as 
#' they are farther away from the closest relevant case (besides itself)
#' at time slice \code{t} (meaning they are more likely to be kept).
#' Farthest away from relevant cases at time slice t: d -> 1.
#' 
#' @param df a data frame
#' @param phi a vector of the relevance values of 
#' \code{df}'s target variable
#' @param rel.thr a relevance threshold above which an 
#' observation is considered relevant
#' @param time the column name of the time-stamp
#' @param sites_sf An sf obejct containing station and IDs and 
#' geometry points of the locations. As an alternative, provide
#' \code{lon}, \code{lat}, and \code{crs}
#' @inheritParams df2site_sf
#'
#' @return A vector of spatially-biased re-sampling weights, scaled
#' to fit within range [0,1].
get_space_wts <- function(df, phi, rel.thr, sites_sf=NULL,
                          lon=NULL, lat=NULL, crs=NULL, site_id, time){
  
  # get sites into right format
  if(is.null(sites_sf)){
    assertthat::assert_that(!is.null(lon), !is.null(lat), !is.null(crs), 
                msg = "Please provide locations object of type sf or 
                CRS code and names of longitude and latitude columns")
    sites_sf <- df2site_sf(df, site_id, lon, lat, crs)
  } 
  
  # create distance matrix
  dists <- get_spatial_dist_mat(sites_sf, site_id)
  max_dist <- max(dists)
  
  timz <- df[[time]]
  space_wts <- vector(mode="numeric", length=nrow(df))
  space_wts <- rep(NA, length(space_wts))
  for(i in 1:length(unique(timz))){
    # get time slice
    t <- unique(timz)[i]
    inds_t <- which(df[[time]]==t)
    
    # get indices of relevant cases at time slice t
    relev_inds <- inds_t[which(phi[inds_t] >= rel.thr)]
    
    # get indices of normal cases
    norm_inds <- setdiff(inds_t, relev_inds)
    if(!length(relev_inds)){
      # if there are no relevant cases, all have max distance 
      # (will be normalized to d=1)
      space_wts[inds_t] <- max_dist # 1
    }else{
      # otherwise, for each case 
      # find minimum distance to relevant case (at time slice t)
      relev_sites <- df[relev_inds, site_id]
      for(i in inds_t){
        s <- df[i, site_id]
        #if(length(setdiff(relev_sites, s))==0) browser()
        
        # if i is the only relevant case, it has maximum distance to other relevant cases
        if((length(unique(relev_sites))==1) && (s %in% relev_sites)){
          d <- max_dist # 1 
        # get minimum distance (to a relevant case)
        }else{
          # check row for site s
          row <- which(rownames(dists)==paste0("SITE_",s))
          # check columns of sites that were relevant at this time slice (except itself)
          cols <- which(colnames(dists) %in% paste0("SITE_", setdiff(relev_sites, s)))
          
          d <- min(dists[row, cols])
        } 
        
        # this is the raw space weight
        space_wts[i] <- d
      }
      
    }
    
    if(t==timz[1]) assertthat::assert_that(all(df[which(!is.na(space_wts)),time]==t))
  }
  
  # overall normal and relevant inds
  norm_inds <- which(phi < rel.thr)
  relev_inds <- which(phi >= rel.thr)
  
  # each group of weights is normalized to scale [0,1]
  space_wts[norm_inds] <- norm_scale(space_wts[norm_inds])
  space_wts[relev_inds] <- norm_scale(space_wts[relev_inds])
  
  space_wts
}

#' Get spatio-temporal re-sampling weights
#'
#' A function that calculates different weights for
#' re-sampling that is temporally and/or spatially biased.
#' 
#' @details \code{phi} gives the target variable's relevance 
#' (higher relevance: phi -> 1; lower relevance: phi -> 0);
#' \code{time_wts} gives the observation's temporally biased
#' re-sampling weight (most recent observations: w -> 1; 
#' oldest: w -> 0.); \code{space_wts} gives the observation's
#' spatially biased re-sampling weight (farthest away from other 
#' relevant cases at time slice: d -> 1.).
#' High \code{time_wts} or \code{space_wts} means the observation is
#' more likely to be kept.
#' 
#' @param form a formula describing the learning task
#' @param df a data frame
#' @param alpha weighting parameter for temporal and spatial
#' re-sampling probabilities. Default 0.5
#' @param beta weighting parameter for spatiotemporal weight and phi for
#' re-sampling probabilities. Default 0.9
#' @param epsilon minimum weight to be added to all observations. 
#' Default 1E-4
#' @inheritParams get_phi
#' @inheritParams get_space_wts
#'
#' @return a data.frame with relevance \code{phi},
#' temporally biased weights \code{time_wts},
#' and spatially biased weights \code{space_wts} for
#' each row in \code{df}.
#'  
#' @seealso \code{\link{get_phi}}, \code{\link{get_time_wts}},
#'  \code{\link{get_space_wts}}.
#'  
#' @export
sample_wts <- function(form, df, phi.control, alpha = 0.5, beta = 0.9, 
                       rel.thr=0.9, 
                       epsilon=1E-4,
         site_id="site_id", time="time", sites_sf = NULL, 
         lon=NULL, lat=NULL, crs = NULL){
  
  # require(assertthat)
  
  assertthat::assert_that(alpha>=0, alpha<=1, msg = "alpha must be between 0 and 1")
    
  # check that there are no NAs in time and space tags
  assertthat::assert_that(!any(is.na(df[[time]])), 
                          !any(is.na(df[[site_id]])),
                          msg = "variables 'time' and 'site_id' cannot contain any NAs")
  
  # check that there are no NAs in target
  y <- stats::model.response(stats::model.frame(form, df, na.action = NULL))
  assertthat::assert_that(!any(is.na(y)), 
              msg = "target variable must not contain any NAs")
  
  # check that either sites_sf or lon/lat are provided
  if(is.null(sites_sf)){
    assertthat::assert_that(!is.null(lon), !is.null(lat), !is.null(crs), 
                msg = "please provide locations object of type sf or 
                CRS code and names of longitude and latitude columns")
    assertthat::assert_that(!any(is.na(df[[lat]])), !any(is.na(df[[lon]])), 
                msg = "variables 'lat' and 'lon' cannot contain any NAs")
  }
  
  # RELEVANCE
  phi <- get_phi(y, phi.control)
  
  # TIME
  timz <- df[[time]]
  time_wts <- get_time_wts(times = timz, phi = phi, rel.thr = rel.thr)
  
  # SPACE
  space_wts <- get_space_wts(df = df, phi = phi, rel.thr = rel.thr, site_id = site_id,
                             sites_sf = sites_sf, lon = lon, lat = lat, time = time, crs = crs)
  
  assertthat::assert_that(length(y)==length(phi),
              length(phi)==length(time_wts),
              length(time_wts)==length(space_wts))
  
  stprob <- data.frame(phi=phi, time_wts=time_wts, space_wts=space_wts)
  stprob$stprob_add <- (alpha*stprob$time_wts+(1-alpha)*stprob$space_wts) + epsilon
  
  stprob
}
