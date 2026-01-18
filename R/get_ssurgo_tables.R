#' Data source is USDA-NRCS Soil Data Access. See package soilDB for more details
#'
#' * If a point is requested then an object of class \sQuote{sf} is returned (for mapunit.shp)
#' with the MUKEY and AREASYMBOL with GEOMETRY type: POINT.
#'
#' * If a the request is for a spatial polygon, then an object of class \sQuote{sf} is returned
#' with gid, mukey and area_ac with GEOMETRY type: POLYGON.
#'
#' @title Retrieve soil profile data and return a (list) with data frames (tables)
#' @description This function does partially what get_ssurgo_soil_profile does, but it
#' returns a list with tables for mapunit, component, chorizon and mapunit.shp (object of class sf)
#' @name get_ssurgo_tables
#' @param lonlat Longitude and latitude vector (e.g. c(-93, 42))
#' @param shift simple mechanism for creating an area of interest by displacing the point indicated in
#' lonlat by some amount of distance (e.g. 300 - in meters)
#' @param aoi area of interest, if supplied the lonlat and shift arguments will be ignored. Should be
#' of class \sQuote{sp::SpatialPolygons} or \sQuote{sf}.
#' @param verbose whether to print messages and warnings to the console default FALSE
#' @param ... additional arguments to be passed to \sQuote{soilDB} functions.
#' @return a list with elements: mapunit, component, chorizon and mapunit.shp
#' @export
#' @examples
#' \dontrun{
#' require(soilDB)
#' require(sp)
#' require(sf)
#' require(spData)
#' ## retrieve data from lon -93, lat = 42
#' stbls <- get_ssurgo_tables(lonlat = c(-93, 42))
#'
#' stbls2 <- get_ssurgo_tables(lonlat = c(-93, 42), shift = 200)
#'
#' }

get_ssurgo_tables <- function(lonlat, shift = -1, aoi, verbose = FALSE, ...){

  if(!requireNamespace("soilDB", quietly = TRUE)){
    warning("The soilDB package is required for this function")
    return(NULL)
  }

  if(!requireNamespace("sf", quietly = TRUE)){
    warning("The sf package is required for this function")
    return(NULL)
  }

  if(!requireNamespace("spData", quietly = TRUE)){
    warning("The spData package is required for this function")
    return(NULL)
  }

  if(!missing(lonlat) && !missing(aoi))
    stop("Either use 'lonlat' or 'aoi', but not both", call. = FALSE)

  if(missing(aoi)){
    lon <- lonlat[1]
    lat <- lonlat[2]

    if(lon < -180 || lon > 180) stop("longitude should be between -180 and 180")
    if(lat < -90 || lat > 90) stop("latitude should be between -90 and 90")

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
  }else{

    ## It appears that coercion is not needed
    ## if(inherits(aoi, "sf")) aoi <- sf::as_Spatial(aoi, "Spatial")

    if(!inherits(aoi, c("SpatialPolygons", "sf", "SpatRaster", "SpatVector")))
      stop("'aoi' should be of class 'SpatialPolygons' or 'sf'.", call. = TRUE)
    spg <- aoi
  }

  res <- .suppressExpression(soilDB::SDA_spatialQuery(spg, what = 'mupolygon', geomIntersection = TRUE), verbose)

  mu.is <- soilDB::format_SQL_in_statement(res$mukey)
  sql <- sprintf("mukey IN %s", mu.is)

  mapunit <- .suppressExpression(soilDB::get_mapunit_from_SDA(sql), verbose)

  cmpnt <- .suppressExpression(soilDB::SDA_query(
      paste0(
          "SELECT compname, comppct_r, slope_r, elev_r,
              drainagecl, taxsubgrp, taxpartsize, taxclname,
              STRING_AGG(CONCAT(geomposmntn, geomposhill, geomposflats, geompostrce), ',') AS geomdesc,
              mapunit.mukey, component.cokey
         FROM component
         INNER JOIN cogeomordesc ON component.cokey = cogeomordesc.cokey
         INNER JOIN cosurfmorphgc ON cogeomordesc.cogeomdkey = cosurfmorphgc.cogeomdkey
         INNER JOIN mapunit ON component.mukey = mapunit.mukey AND mapunit.", sql, "
         GROUP BY component.cokey, compname, comppct_r, slope_r, elev_r,
              drainagecl, taxsubgrp, taxpartsize, taxclname, mapunit.mukey"
      )
  ), verbose)

  chrzns <- .suppressExpression(soilDB::SDA_query(
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

  ### Component ###
  names(cmpnt) <- gsub("_", ".", names(cmpnt), fixed = TRUE)
  cmpnt$geomdesc <- cmpnt$geompos

  ### Getting interpretation information ----
  cointer <- .suppressExpression(
      soilDB::get_SDA_interpretation(
          "NCCPI - National Commodity Crop Productivity Index (Ver 3.0)",
          method = "None",
          WHERE = sql
      ),
      verbose
  )

  ### retrieve the state information ###
  ## This is unnecessary, but it works for now
  states <- spData::us_states ## This object contains states and their geometries
  states <- sf::st_transform(states, crs = 3857)
  pspg <- sf::st_transform(sf::st_as_sf(spg), crs = 3857)
  ii <- as.integer(sf::st_intersects(pspg, states))
  cmpnt$state <- states[["NAME"]][[ii]]

  ### Chorizon ###
  names(chrzns) <- gsub("_", ".", names(chrzns), fixed = TRUE)
  ### Things missing from horizons: hzthk.r, partdensity, wsatiated.r, wfifteenbar.r, wtenthbar.r, wthirdbar.r,
  if(sum(grepl("partdensity", names(chrzns))) == 0) chrzns$partdensity <- NA
  if(sum(grepl("hzthk", names(chrzns))) == 0) chrzns$hzthk.r <- NA
  if(sum(grepl("wsatiated", names(chrzns))) == 0) chrzns$wsatiated.r <- NA
  if(sum(grepl("wfifteenbar", names(chrzns))) == 0) chrzns$wfifteenbar.r <- NA
  if(sum(grepl("wthirdbar", names(chrzns))) == 0) chrzns$wthirdbar.r <- NA

  if(shift <= 0 && missing(aoi)){
    spg.sf <- sf::st_as_sf(spg)
    spg.sf[["MUKEY"]] <- res$mukey
    spg.sf[["AREASYMBOL"]] <- mapunit$areasymbol
    mapunit.shp <- spg.sf
  }else{
    mapunit.shp <- sf::st_as_sf(res)
  }

  ans <- list(mapunit = mapunit, component = cmpnt, cointerp = cointer,
              chorizon = chrzns, mapunit.shp = mapunit.shp)
  ans
}
