########################################################################
############ SDMs using MAXENT and ENMwizard ###########
# Rversion = 4.3.2 
# Gabriela Alves Ferreira (gabriela-alves77@hotmail.com)
# Adjusted for GitHub Repository Structure

# ==============================================================================
# 0. Load Packages
# ==============================================================================
library(devtools)
library(usethis)
library(ENMeval)
library(ENMwizard)
library(dismo)
library(rgdal)
library(raster)
library(rgeos)
library(spThin)
library(rJava)

# ==============================================================================
# 1. Occurrence Data Preparation
# ==============================================================================
## Read the consolidated and thinned occurrence dataset from the 'data' folder
spp <- read.csv("data/spps_occ_env_thin_f.csv", head = TRUE, sep = ",")

## Make it a list (required by ENMwizard)
spp.occ.list <- lapply(unique(spp$sp), function(i, x){
  x[x$sp==i, c("long", "lat")]
}, x=spp)

names(spp.occ.list) <- unique(spp$sp)
head(spp.occ.list)

## Create occurrence polygons to crop rasters prior to modeling
occ.polys <- set_calibarea_b(spp.occ.list)

## Create buffer around occurrence polygons using 1.5 degrees
occ.b <- buffer_b(occ.polys, width = 1.5)


# ==============================================================================
# 2. Environmental Layers (Current)
# ==============================================================================
## Path to the dataset folder (Ensure you have downloaded WorldClim data here)
filePath2 <- "data/env_layers/current/" 

## Import all .tif files
vars <- list.files(path = filePath2, pattern='.tif', all.files=TRUE)
current <- stack(paste0(filePath2, vars))

names(current) <- sub(".*[_$]", "bio_", names(current))

## Cut environmental variables for each species
pred.cut <- cut_calibarea_b(occ.b, current)

## Remove correlated variables (cutoff = 0.75)
pred.cut <- select_vars_b(pred.cut, cutoff=.75, names.only = FALSE)


# ==============================================================================
# 3. Model Tuning & Calibration (ENMeval)
# ==============================================================================
## Test combinations of Feature Classes and Regularization Multipliers
## Method based on the number of occs: jackknife (<=15), block (>15)
method <- ifelse(sapply(spp.occ.list, nrow) <= 15, "jackknife", "block")

ENMeval.res.lst <- ENMevaluate_b(spp.occ.list, pred.cut, 
                                 RMvalues = seq(0.5, 5, 0.5),
                                 fc = c("L", "P", "Q", "LP", "LQ", "PQ", "LPQ"),
                                 method="block", algorithm="maxent.jar", 
                                 numCores = 8, parallel = TRUE)

## Save tuning results
saveRDS(ENMeval.res.lst, "outputs/models/ENMeval.res.lst.rds")

## Model fitting (calibration) using optimal parameters
mxnt.mdls.preds.lst <- calib_mdl_b(ENMeval.o.l = ENMeval.res.lst, 
                                   a.calib.l = pred.cut, 
                                   occ.l = spp.occ.list, 
                                   mSel = c("EBPM"))

rm(ENMeval.res.lst)
saveRDS(mxnt.mdls.preds.lst, "outputs/models/mxnt_mdls_preds_list.rds")


# ==============================================================================
# 4. Future Projections Preparation
# ==============================================================================
## Read data for future projection (2050/2070) - UPDATE PATHS AS NEEDED
path_future <- "data/env_layers/future/"

r_CCSM4_45 <- brick(paste0(path_future, "wc2.1_2.5m_bioc_CCSM4_ssp245_2041-2060.tif"))
r_CCSM4_85 <- brick(paste0(path_future, "wc2.1_2.5m_bioc_CCSM4_ssp585_2041-2060.tif"))
r_MIROC_45 <- brick(paste0(path_future, "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060.tif"))
r_MIROC_85 <- brick(paste0(path_future, "wc2.1_2.5m_bioc_MIROC6_ssp585_2041-2060.tif"))
r_MPI_45 <- brick(paste0(path_future, "wc2.1_2.5m_bioc_MPI_ssp245_2041-2060.tif"))
r_MPI_85 <- brick(paste0(path_future, "wc2.1_2.5m_bioc_MPI_ssp585_2041-2060.tif"))

## Set CRS for projections
crs_wgs84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
proj4string(r_MPI_45) <- crs_wgs84
proj4string(r_MPI_85) <- crs_wgs84
proj4string(r_MIROC_45) <- crs_wgs84
proj4string(r_MIROC_85) <- crs_wgs84
proj4string(r_CCSM4_45) <- crs_wgs84
proj4string(r_CCSM4_85) <- crs_wgs84

## Fix raster names to match baseline
names(r_MPI_45) <- sub(".*[_$]", "bio_", names(r_MPI_45))
names(r_MPI_85) <- sub(".*[_$]", "bio_", names(r_MPI_85))
names(r_MIROC_45) <- sub(".*[_$]", "bio_", names(r_MIROC_45))
names(r_MIROC_85) <- sub(".*[_$]", "bio_", names(r_MIROC_85))
names(r_CCSM4_45) <- sub(".*[_$]", "bio_", names(r_CCSM4_45))
names(r_CCSM4_85) <- sub(".*[_$]", "bio_", names(r_CCSM4_85))

## Group predictors
predictors.l <- list(current = current,
                     r_MPI_SSP_45 = r_MPI_45,
                     r_MPI_SSP_85 = r_MPI_85,
                     r_MIROC_SSP_45 = r_MIROC_45,
                     r_MIROC_SSP_85 = r_MIROC_85,
                     r_CCSM4_SSP_45 = r_CCSM4_45,
                     r_CCSM4_SSP_85 = r_CCSM4_85)

# ==============================================================================
# 5. Model Projections & Consensus
# ==============================================================================
## Preparing projection area using the provided shapefile
para <- rgdal::readOGR("data/para_state.shp") 
proj4string(para) <- crs_wgs84

## Cut projection area
pa.env.proj.l <- cut_projarea_rst_mscn_b(para, predictors.l, occ.polys, mask = TRUE)

## Project models
mxnt.mdls.preds.cf <- proj_mdl_b(mxnt.mdls.preds.lst, a.proj.l = pa.env.proj.l)

## Grouping codes for Consensus
ssp <- c("SSP_45", "SSP_85")
groups <- list(ssp)
clim.scn.nms <- names(mxnt.mdls.preds.cf[[1]]$mxnt.preds)

## Compute consensual projections
mxnt.mdls.preds.cf_cons <- consensus_scn_b(mcmp.l= mxnt.mdls.preds.cf, 
                                           groups = list(ssp), ref = "current", 
                                           save = TRUE, numCores = 8)
saveRDS(mxnt.mdls.preds.cf_cons, "outputs/models/mxnt.mdls.preds.cf_cons.rds")

# ==============================================================================
# 6. Binary Thresholding & Suitable Area (TSA)
# ==============================================================================
## Applying thresholds
mods.thrshld.lst <- thrshld_b(mxnt.mdls.preds.cf_cons, thrshld.i = c(4,5), numCores = 7)
saveRDS(mods.thrshld.lst, "outputs/models/mods_thrshld_lst.rds")

## Compute species' total suitable area
get_tsa_b(mods.thrshld.lst)