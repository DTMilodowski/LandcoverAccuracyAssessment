"""
This library hosts a bunch of functions to accompany accuracy_assessment.py.  For example,
scripts to load in geotiffs etc.
"""
import numpy as np
from osgeo import gdal
import sys
import osr
def get_nodata_mask(MASK_FILE):
    #First step is to load the image, so that missing data can be added to the
    # mask
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    try:
        ds = gdal.Open(MASK_FILE)
    except RuntimeError, e:
        print 'unable to open ' + MASK_FILE
        print e
        sys.exit(1)

    mask= np.array(ds.GetRasterBand(1).ReadAsArray())    
    
    ds = None
    return mask

def load_GeoTIFF_band(File,band_number=1):
    
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    try:
        image_ds = gdal.Open(File)
    except RuntimeError, e:
        print 'unable to open ' + File
        print e
        sys.exit(1)
        
    source_band = image_ds.GetRasterBand(band_number)
    if source_band is None:
        print "BAND MISSING"
        #continue
    image = np.array(image_ds.GetRasterBand(band_number).ReadAsArray(),dtype=np.float64)
    return image

def load_GeoTIFF_band_and_georeferencing(File,band_number=1):
    
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    try:
        image_ds = gdal.Open(File)
    except RuntimeError, e:
        print 'unable to open ' + File
        print e
        sys.exit(1)
        
    source_band = image_ds.GetRasterBand(band_number)
    if source_band is None:
        print "BAND MISSING"
        #continue

    image = np.array(image_ds.GetRasterBand(band_number).ReadAsArray(),dtype=np.float64)
    geotrans = image_ds.GetGeoTransform()
    coord_sys = image_ds.GetProjectionRef()

    return image, geotrans, coord_sys


def save_GeoTIFF_band(array, geoTrans, coord_sys, OUT_FILE, nodata=-9999):
    array[np.isnan(array)]=nodata
    (NRows,NCols)=array.shape
    # Write output files
    print "\tWriting geotiff"

    outdriver = gdal.GetDriverByName('GTiff')
    outdriver.Register()
    # first set all the relevant geospatial information
    dataset_ds = outdriver.Create( OUT_FILE, NCols, NRows, 1, gdal.GDT_Int16)
    dataset_ds.SetGeoTransform( geoTrans )
    CoordSys= osr.SpatialReference()
    CoordSys.ImportFromWkt( coord_sys )
    dataset_ds.SetProjection( CoordSys.ExportToWkt() )
    # write array
    dataset_ds.GetRasterBand(1).WriteArray( array )
    dataset_ds.GetRasterBand(1).SetNoDataValue( nodata )
    
    dataset_ds = None


def load_RapidEye_forest_loss(ForestLossEarly_file,ForestLossLate_file):
    
    # Get ForestLoss arrays
    YEAR_array_Early = load_GeoTIFF_band(ForestLossEarly_file,1)
    MONTH_array_Early = load_GeoTIFF_band(ForestLossEarly_file,2)
    DAY_array_Early = load_GeoTIFF_band(ForestLossEarly_file,3)
    
    YEAR_array_Late = load_GeoTIFF_band(ForestLossLate_file,1)
    MONTH_array_Late = load_GeoTIFF_band(ForestLossLate_file,2)
    DAY_array_Late = load_GeoTIFF_band(ForestLossLate_file,3)

    return YEAR_array_Early, MONTH_array_Early, DAY_array_Early, YEAR_array_Late, MONTH_array_Late, DAY_array_Late


def convert_RapidEye_forest_loss_to_thematic_map(ForestLossEarly,ForestLossLate,Mask,StartYear=2010,EndYear=2015):
    (NRows,NCols) = ForestLossEarly.shape
    ChangeMap=np.zeros((NRows,NCols))
    ChangeMap[Mask==0]=np.nan

    temp = ForestLossEarly.copy()
    temp[ForestLossEarly<StartYear]=0
    temp[ForestLossEarly==1]=1
    ForestLossEarly=temp.copy()

    temp = ForestLossLate.copy()
    temp[ForestLossLate<StartYear]=0
    temp[ForestLossLate==1]=1
    ForestLossLate=temp.copy()
    
    ForestLossEarly[ForestLossEarly>EndYear]=1
    ForestLossLate[ForestLossLate>EndYear]=1

    # Now create "Date" array from Hansen array
    for i in range(0,NRows):
        progress = '%.2f' % (float(i)/float(NRows)*100)
        sys.stdout.write('\r'+progress+'%',)
        sys.stdout.flush()
        for j in range(0,NCols):
            if Mask[i,j]==1:
                # Check Rapideye
                # is pixel already deforested in both early and late?
                if ForestLossEarly[i,j] == 0 and ForestLossLate[i,j] == 0:
                    ChangeMap[i,j]=0
                # is pixel forested throughout?
                elif ForestLossEarly[i,j] == 1 and ForestLossLate[i,j] == 1:
                    ChangeMap[i,j]=1
                # Is change recorded in both RapidEye-based maps?
                elif ForestLossEarly[i,j] > 1 and ForestLossLate[i,j] > 1:
                    ChangeMap[i,j]=2
                # Otherwise, change period confidence level is undefined
                else:
                    ChangeMap[i,j]=np.nan
    return ChangeMap

def get_class_areas_from_RapidEye_pair(ForestLossEarly,ForestLossLate,Mask,StartYear=2010,EndYear=2015):
    class_count = np.zeros(3)
    temp = ForestLossEarly.copy()
    temp[ForestLossEarly<StartYear]=0
    temp[ForestLossEarly==1]=1
    temp[Mask==0]=np.nan
    ForestLossEarly=temp.copy()

    temp = ForestLossLate.copy()
    temp[ForestLossLate<=StartYear]=0
    temp[ForestLossLate==1]=1
    temp[Mask==0]=np.nan
    ForestLossLate=temp.copy()

    class_count[0]=float(np.sum(ForestLossEarly==0) + np.sum(ForestLossLate==0))/2
    class_count[1]=float(np.sum(ForestLossEarly==1) + np.sum(ForestLossLate==1))/2
    class_count[2]=float(np.sum(ForestLossEarly>1) + np.sum(ForestLossLate>1))/2

    return class_count

def get_class_areas_from_changemap(changemap,mask):
    class_count = np.zeros(3)
    temp = changemap.copy()
    temp[mask==0]=np.nan
    ForestLossEarly=temp.copy()


    class_count[0]=float(np.sum(ForestLossEarly==0) + np.sum(ForestLossLate==0))/2
    class_count[1]=float(np.sum(ForestLossEarly==1) + np.sum(ForestLossLate==1))/2
    class_count[2]=float(np.sum(ForestLossEarly==2) + np.sum(ForestLossLate==2))/2

    return class_count
