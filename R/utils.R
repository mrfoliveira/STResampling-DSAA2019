#' Calculate spatial distance matrix
#'
#' A function that calculates the geographical distance
#' matrix between the locations of an \code{sf} object. 
#' @param sites_sf an \code{sf} object with the geograhic
#' information of the locations (as returned by \code{df2site_sf})
#' @param site_id the column name of the location ID
#'
#' @return a matrix of distances. Row and column names
#' are a concatenation of "SITE_" and the location IDs.
#' 
#' @seealso \code{\link{df2site_sf}}
#' 
#' @importFrom lwgeom st_geod_distance
#' @importFrom sf st_distance
#' 
get_spatial_dist_mat <- function(sites_sf, site_id){
  
  assertthat::assert_that(any("sf" %in% class(sites_sf)))
  suppressPackageStartupMessages( requireNamespace("lwgeom", quietly = TRUE) )
  
  # unique location ids
  sites_sf <- unique(sites_sf)
  sids <- sites_sf[[site_id]]
  # calculate distance matrix using st_distance
  dists <- sf::st_distance(sites_sf, sites_sf)
  dists <- apply(dists, 2, as.numeric)
  colnames(dists) <- paste0("SITE_", sids)  
  rownames(dists) <- paste0("SITE_", sids)
  
  assertthat::assert_that(all(!is.na(dists)), all(dists>=0))
  
  dists
}

#' Create an sf object of available sites
#' 
#' Extracts the location information from a data frame
#' and transforms into a \code{sf} object.
#' @param df a data frame of the data set
#' @param site_id the name of the column containing location IDs
#' @param lon the name of the column containing the location's longitude
#' @param lat the name of the column containing the location's latitude
#' @param crs the code for the Coordinate Reference System
#'
#' @return a sf object, containing the geographic information for
#' each location in \code{df}
#' @seealso \code{\link[sf]{st_as_sf}}
#' 
df2site_sf <- function(df, site_id, lon, lat, crs){
  
  # not sf class
  if(!("sf" %in% class(df))){
    assertthat::assert_that(is.numeric(df[[lon]])) #, all(df[[lon]] > -180), all(df[[lon]] < 180),
    #msg = "variable 'lon' must be numeric between -180 and 180")
    assertthat::assert_that(is.numeric(df[[lat]]))#, all(df[[lat]] > -180), all(df[[lon]] < 180),
    #msg = "variable 'lon' must be numeric between -180 and 180")
    
    # create dataset for locations
    sites_df <- df[which(!duplicated(df[[site_id]])), c(site_id, lon, lat)]
    sites_sf <- sf::st_as_sf(sites_df, coords=c(lon, lat), crs=crs)
  }else{
    sites_sf <- df[which(!duplicated(df[[site_id]])), site_id]
  }
  sites_sf
}

#' Feature scaling
#' 
#' Normalize values to be within the range between [0,1].
#' @param x a vector of values
#' @return a scaled vector
#' 
#' @export
norm_scale <- function(x){
  if(min(x)!=max(x)) (x - min (x)) / ( max(x) - min(x) )
  else 1
} 

#' Get response values of a dataset from a formula
#' 
#' @param formula learning formula
#' @param data data set to get the target values from
#' @param na what action to perform if NAs are present. Default is na.fail
#' @return A vector of the target values.
responseValues <- function (formula, data, na = NULL) 
  stats::model.response(stats::model.frame(formula, data, na.action = na))

#' Shuffle values/rows
#' 
#' Shuffle the values or rows of a vector or data frame
#' @param x a vector or data frame
#' @return a vector or data frame
shuffle <- function(x){
  if(is.null(dim(x))) x[sample(NROW(x))]
  else x[sample(NROW(x)),]
}

#' Cut into folds
#' 
#' Assigns rows of a data frame into folds for cross-validation.
#' @param x a data.frame
#' @param nfolds number of folds
#' @return a vector with the fold assignment of each row
cv_folds <- function(x, nfolds) {
  cut(seq_len(NROW(x)), breaks = nfolds, labels = FALSE)
}

