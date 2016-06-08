import numpy as np
from osgeo import gdal
import os
import osr
import sys
from netCDF4 import Dataset
from datetime import date
import calendar

import accuracy_assessment as acc
import auxilliary_functions as aux

Mask_file = '/home/dmilodow/DataStore_DTM/IPSP/DATA/SVM_CLASSIFICATION_VERSION2/Acre_mask_with_river.tif'
RapidEye_early_file = '/home/dmilodow/DataStore_DTM/IPSP/DATA/SVM_CLASSIFICATION_VERSION2/RECLASSIFICATION/TEMPORAL_RECLASS/Acre_deforestation_date_early_opt.tif'   
RapidEye_late_file = '/home/dmilodow/DataStore_DTM/IPSP/DATA/SVM_CLASSIFICATION_VERSION2/RECLASSIFICATION/TEMPORAL_RECLASS/Acre_deforestation_date_late_opt.tif'
SAVEDIR = "../output"
FILENAME = "Acre_strat_rand_sample_pts"


mask = aux.get_nodata_mask(Mask_file)
FLe,geoTrans,coord = aux.load_GeoTIFF_band_and_georeferencing(RapidEye_early_file)
FLl,geoTrans,coord = aux.load_GeoTIFF_band_and_georeferencing(RapidEye_late_file)
ChangeMapAcre = aux.convert_RapidEye_forest_loss_to_thematic_map(FLe,FLl,mask)
aux.save_GeoTIFF_band(ChangeMapAcre,geoTrans,coord,"../output/AcreChangeMap.tif")

XMin=geoTrans[0]
YMax=geoTrans[3]
XResolution=5.
YResolution=5.

Classkeys=np.arange(0,3)
PredictedClassUA = np.asarray([0.95,0.95,0.80])

sample_points = acc.retrieve_stratified_random_sample(ChangeMapAcre,Classkeys,PredictedClassUA,XMin,YMax,XResolution,YResolution)

acc.write_sample_points_to_csv(sample_points,FILENAME,SAVEDIR)



Mask_file = '/home/dmilodow/DataStore_DTM/IPSP/DATA/SVM_CLASSIFICATION_VERSION2/Rondonia_mask.tif'
RapidEye_early_file = '/home/dmilodow/DataStore_DTM/IPSP/DATA/SVM_CLASSIFICATION_VERSION2/RECLASSIFICATION/TEMPORAL_RECLASS/RONDONIA_deforestation_date_early_opt.tif'   
RapidEye_late_file = '/home/dmilodow/DataStore_DTM/IPSP/DATA/SVM_CLASSIFICATION_VERSION2/RECLASSIFICATION/TEMPORAL_RECLASS/RONDONIA_deforestation_date_late_opt.tif'
SAVEDIR = "../output"
FILENAME = "Rondonia_strat_rand_sample_pts"

mask = aux.get_nodata_mask(Mask_file)
FLe,geoTrans,coord = aux.load_GeoTIFF_band_and_georeferencing(RapidEye_early_file)
FLl,geoTrans,coord = aux.load_GeoTIFF_band_and_georeferencing(RapidEye_late_file)
ChangeMapAcre = aux.convert_RapidEye_forest_loss_to_thematic_map(FLe,FLl,mask)
aux.save_GeoTIFF_band(ChangeMapAcre,geoTrans,coord,"../output/RondoniaChangeMap.tif")

XMin=geoTrans[0]
YMax=geoTrans[3]
XResolution=5.
YResolution=5.

Classkeys=np.arange(0,3)
PredictedClassUA = np.asarray([0.95,0.95,0.80])

sample_points = acc.retrieve_stratified_random_sample(ChangeMapAcre,Classkeys,PredictedClassUA,XMin,YMax,XResolution,YResolution)

acc.write_sample_points_to_csv(sample_points,FILENAME,SAVEDIR)

