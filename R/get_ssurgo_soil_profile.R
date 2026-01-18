#' Data source is USDA-NRCS Soil Data Access. See package soilDB for more details
#'
#' @title Retrieve soil profile data and convert it to an object of class \sQuote{soil_profile}
#' @description Generate a synthetic soil profile based on the information in SSURGO database
#' @name get_ssurgo_soil_profile
#' @param lonlat Longitude and latitude vector (e.g. c(-93, 42))
#' @param shift simple mechanism for creating an area of interest by displacing the point indicated in
#' lonlat by some amount of distance (e.g. 300 - in meters)
#' @param nmapunit number of mapunits to select (see \code{\link{ssurgo2sp}})
#' @param nsoil number of soils to select (see \code{\link{ssurgo2sp}}). If the
#' number of soils is negative or NA it will fetch all the soils in the mapunit
#' @param xout see \code{\link{ssurgo2sp}}
#' @param soil.bottom see \code{\link{ssurgo2sp}}
#' @param method interpolation method see \code{\link{ssurgo2sp}}
#' @param nlayers number for layer for the new soil profile
#' @param check whether to check for reasonable values using \code{\link{check_apsimx_soil_profile}}.
#' TRUE by default. If \sQuote{fix} is TRUE, it will be applied only after the fix attempt.
#' @param fix whether to fix compatibility between saturation and bulk density (default is FALSE).
#' @param verbose default FALSE. Whether to print messages.
#' @param xargs additional arguments passed to \code{\link{apsimx_soil_profile}} function.
#' @return this function will always return a list. Each element of the list will
#' be an object of class \sQuote{soil_profile}
#' @export
#' @examples
#' \dontrun{
#' require(soilDB)
#' require(sp)
#' require(sf)
#' require(spData)
#' require(ggplot2)
#' ## Soil inforation for a single point
#' sp <- get_ssurgo_soil_profile(lonlat = c(-93, 42))
#' ## The initial attempt throws warnings, so better to use 'fix'
#' sp <- get_ssurgo_soil_profile(lonlat = c(-93, 42), fix = TRUE)
#' plot(sp[[1]])
#' plot(sp[[1]], property = "water")
#' ## Add initial water
#' iwat <- initialwater_parms(Thickness = sp[[1]]$soil$Thickness,
#'                            InitialValues = sp[[1]]$soil$DUL * 0.8)
#' sp[[1]]$initialwater <- iwat
#' plot(sp[[1]], property = "initialwater")
#' }
#'
#'

get_ssurgo_soil_profile <- function(lonlat, shift = -1,
                                    nmapunit = 1, nsoil = 1,
                                    xout = NULL, soil.bottom = 200,
                                    method = c("constant", "linear"),
                                    nlayers = 10,
                                    check = TRUE,
                                    fix = FALSE,
                                    verbose = FALSE,
                                    xargs = NULL){

  if(!requireNamespace("soilDB", quietly = TRUE)){
    stop("The soilDB package is required for this function")
    return(NULL)
  }

  if(!requireNamespace("sf", quietly = TRUE)){
    stop("The sf package is required for this function")
    return(NULL)
  }

  if(!requireNamespace("spData", quietly = TRUE)){
    stop("The spData package is required for this function")
    return(NULL)
  }

  if(length(lonlat) != 2 || !is.numeric(lonlat))
    stop("lonlat should be a vector with length equal to 2")

  lon <- lonlat[1]
  lat <- lonlat[2]

  ## Determine if the location is in the US
  if(requireNamespace("maps", quietly = TRUE)){
    country <- maps::map.where(x = lon, y = lat)
    if(country != "USA" || is.na(country))
      stop("These coordinates do not correspond to a location in the USA. \n Did you specify the coordinates correctly?")
  }

  if(shift <= 0){
    spg <- sf::st_as_sf(data.frame(x = lon, y = lat),
                        coords = c("x", "y"),
                        crs = "EPSG:4326")
  }else{
    shift <- (shift / 111) * 0.001 ## This is now in degrees
    lonlat.mat <- rbind(lonlat, ##root
                        lonlat + c(shift * 0.75, 0), ## x = 1, y = 0
                        lonlat + c(shift * 0.75, shift), ## x = 1, y = 1
                        lonlat + c(0, shift), ## x = 0, y = 1
                        lonlat) ## back to root
    rownames(lonlat.mat) <- NULL
    ## the previous matrix is a rectangle
    spg <- sf::st_as_sf(
      data.frame(geometry = sprintf(
        "POLYGON ((%s))", paste0(paste(lonlat.mat[, 1], lonlat.mat[, 2]), collapse = ",")
      )),
      wkt = "geometry",
      crs = "EPSG:4326"
    )
  }

  res <- .suppressExpression(soilDB::SDA_spatialQuery(spg, what = 'mupolygon', geomIntersection = TRUE), verbose)

  mu.is <- soilDB::format_SQL_in_statement(res$mukey)
  sql <- sprintf("mukey IN %s", mu.is)

  fSDA_mapunit <- .suppressExpression(soilDB::SDA_query(
    paste0(
      "SELECT legend.areasymbol, musym, muname, muacres,
              farmlndcl, iacornsr,
              legend.lkey, mukey
         FROM mapunit
         INNER JOIN legend ON legend.lkey = mapunit.lkey
         WHERE mapunit.", sql
    )
  ), verbose)

  fSDA_component <- .suppressExpression(soilDB::SDA_query(
    paste0(
      "SELECT compname, comppct_r, slope_r, elev_r,
              drainagecl, taxsubgrp, taxpartsize, taxclname,
              mapunit.mukey,cokey
         FROM component
         INNER JOIN mapunit ON component.mukey = mapunit.mukey
             AND mapunit.", sql
    )
  ), verbose)

  fSDA_chorizon <- .suppressExpression(soilDB::SDA_query(
    paste0(
      "SELECT hzname, hzdept_r, hzdepb_r, hzdepb_r - hzdept_r AS hzthk_r,
              sandtotal_r, silttotal_r, claytotal_r, om_r, partdensity,
              ksat_r, awc_r, wthirdbar_r, wfifteenbar_r, wsatiated_r,
              ph1to1h2o_r, caco3_r,
              component.cokey, chkey
         FROM component
         LEFT JOIN chorizon ON chorizon.cokey = component.cokey
         WHERE component.", sql
    )
  ), verbose)

  mapunit <- fSDA_mapunit

  ### Component ###
  cmpnt <- fSDA_component
  names(cmpnt) <- gsub("_", ".", names(cmpnt), fixed = TRUE)
  cmpnt$geomdesc <- NA #TODO cmpnt$geompos

  ## Retrieve the state from the areasymbol
  if(shift <= 0 || length(unique(mapunit$areasymbol)) == 1){
    cmpnt$state <- unique(strtrim(mapunit$areasymbol, 2))
  }else{
    cmpnt$state <- NA
    warning("This area includes more than one state.
            I have not though about how to get the state in this case. Please submit an issue
            with a reproducible example to https://github.com/femiguez/apsimx/issues")
  }

  ### Chorizon ###
  chrzns <- fSDA_chorizon
  names(chrzns) <- gsub("_", ".", names(chrzns), fixed = TRUE)
  ### Things missing from horizons: hzthk.r, partdensity, wsatiated.r, wfifteenbar.r, wtenthbar.r, wthirdbar.r,
  if(sum(grepl("partdensity", names(chrzns))) == 0) chrzns$partdensity <- NA
  if(sum(grepl("hzthk", names(chrzns))) == 0) chrzns$hzthk.r <- NA
  if(sum(grepl("wsatiated", names(chrzns))) == 0) chrzns$wsatiated.r <- NA
  if(sum(grepl("wfifteenbar", names(chrzns))) == 0) chrzns$wfifteenbar.r <- NA
  if(sum(grepl("wthirdbar", names(chrzns))) == 0) chrzns$wthirdbar.r <- NA

  if(shift <= 0){
    spg.sf <- spg
    spg.sf[["MUKEY"]] <- res$mukey
    spg.sf[["AREASYMBOL"]] <- mapunit$areasymbol
    mapunit.shp <- spg.sf
  }else{
    mapunit.shp <- sf::st_as_sf(res)
  }

  sp0 <- ssurgo2sp(mapunit = mapunit, component = cmpnt,
                   chorizon = chrzns, mapunit.shp = mapunit.shp,
                   nmapunit = nmapunit, nsoil = nsoil, xout = xout,
                   soil.bottom = soil.bottom, method = method, nlayers = nlayers,
                   verbose = verbose)

  ans <- vector("list", length(sp0))

  for(i in seq_along(sp0)){
    metadata <- attributes(sp0[[i]])
    metadata$DataSource <- paste("SSURGO (https://sdmdataaccess.nrcs.usda.gov/) through R package soilDB, R package apsimx function ssurgo2sp. Timestamp",Sys.time())
    metadata$names <- NULL; metadata$class <- NULL; metadata$row.names <- NULL;

    if(fix){
      icheck <- check
      check <- FALSE
    }

    if(!is.null(xargs)){
      if(!is.null(xargs$crops)){
        crops <- xargs$crops
      }
    }else{
      crops <- c("Maize", "Soybean", "Wheat")
    }

    asp <- apsimx_soil_profile(nlayers = nlayers,
                               Thickness = sp0[[i]]$Thickness * 10,
                               BD = sp0[[i]]$BD,
                               AirDry = sp0[[i]]$AirDry,
                               LL15 = sp0[[i]]$LL15,
                               DUL = sp0[[i]]$DUL,
                               SAT = sp0[[i]]$SAT,
                               KS = sp0[[i]]$KS,
                               Carbon = sp0[[i]]$Carbon,
                               crop.LL = sp0[[i]]$LL15,
                               ParticleSizeClay = sp0[[i]]$ParticleSizeClay,
                               ParticleSizeSilt = sp0[[i]]$ParticleSizeSilt,
                               ParticleSizeSand = sp0[[i]]$ParticleSizeSand,
                               soil.bottom = soil.bottom,
                               metadata = metadata,
                               check = check,
                               crops = crops)

    if(fix){
      asp <- fix_apsimx_soil_profile(asp, verbose = verbose)
      if(icheck)
        check_apsimx_soil_profile(asp)
    }

    ans[[i]] <- asp
  }

  return(ans)
}
