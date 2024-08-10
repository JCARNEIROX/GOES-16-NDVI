# Training: Python and GOES-R Imagery: Script 8 - Functions for download yyyymmddhhmn from AWS
# -----------------------------------------------------------------------------------------------------------
# Required modules
import os  # Miscellaneous operating system interfaces
import colorsys  # To make convertion of colormaps
import boto3  # Amazon Web Services (AWS) SDK for Python
from botocore import UNSIGNED  # boto3 config
from botocore.config import Config  # boto3 config
import math  # Mathematical functions
from osgeo import osr  # Python bindings for GDAL
from datetime import timedelta, datetime  # Basic Dates and time types
import cartopy, cartopy.crs as ccrs  # Plot maps
import cartopy.io.shapereader as shpreader  # Import shapefiles
import cartopy.feature as cfeature # features
import numpy as np  # Scientific computing with Python
from osgeo import gdal  # Python bindings for GDAL
import matplotlib.pyplot as plt  # Plotting library
import matplotlib
from netCDF4 import Dataset
from pyorbital import astronomy
from pyspectral.rayleigh import Rayleigh                     # Correção atmosférica no espectro visível 
from pyorbital.astronomy import get_alt_az
from pyorbital.orbital import get_observer_look


#======================================================================================================
# Input and output directories
dir_in = "G:\Drives compartilhados\GOES16_FDCF\Programacao\dir_in\\"
dir_main = "G:\Drives compartilhados\GOES16_FDCF\Programacao\\"
dir_out = dir_main + "output\\"
dir_libs = dir_main + "libs\\"
dir_shapefiles = dir_main + "shapefiles\\"
dir_colortables = dir_main + "colortables\\"
dir_logos = dir_main + "logos\\"


#Desired Extent
extent_br = [-90.0, -40.0, -20.0, 10.0]  # Min Lon, Max lon, Min lat, Max lat Brasil
extent_pantanal = [-59.3, -54.7, -22.3, -15.5] # Min Lon, Max lon, Min lat, Max lat Pantanal
#======================================================================================================

def loadCPT(path):
    global btemp, gtemp, xtemp, rtemp, bb, gg, rr, i
    try:
        f = open(path)
    except:
        print("File ", path, "not found")
        return None

    lines = f.readlines()

    f.close()

    x = np.array([])
    r = np.array([])
    g = np.array([])
    b = np.array([])

    colorModel = 'RGB'

    for l in lines:
        ls = l.split()
        if l[0] == '#':
            if ls[-1] == 'HSV':
                colorModel = 'HSV'
                continue
            else:
                continue
        if ls[0] == 'B' or ls[0] == 'F' or ls[0] == 'N':
            pass
        else:
            x = np.append(x, float(ls[0]))
            r = np.append(r, float(ls[1]))
            g = np.append(g, float(ls[2]))
            b = np.append(b, float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

        x = np.append(x, xtemp)
        r = np.append(r, rtemp)
        g = np.append(g, gtemp)
        b = np.append(b, btemp)

    if colorModel == 'HSV':
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
        r[i] = rr
        g[i] = gg
        b[i] = bb

    if colorModel == 'RGB':
        r = r / 255.0
        g = g / 255.0
        b = b / 255.0

    xNorm = (x - x[0]) / (x[-1] - x[0])

    red = []
    blue = []
    green = []

    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])

    colorDict = {'red': red, 'green': green, 'blue': blue}

    return colorDict

#função Problema Servidor
def download_cmi(yyyymmddhhmn, band, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'
    product_name = 'ABI-L2-CMIPF'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6C{int(band):02.0f}_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Band-{band}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]
            #Alteração problema servidor
            new_file_name = file_name.replace('OR','CG')
            # print(file_name)
            # print(new_file_name)

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')

            # Alteração problema servidor
            old_file = os.path.join(path_dest, f'{file_name}.nc')
            new_file = os.path.join(path_dest, f'{new_file_name}.nc')
            os.rename(old_file, new_file)

    return f'{file_name}.nc'
# -----------------------------------------------------------------------------------------------------------
def download_CMI(yyyydddhhmn, band, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%j')
    hour = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%H')
    min = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'
    product_name = 'ABI-L2-CMIPF'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6C{int(band):02.0f}_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyydddhhmn}, Band-{band}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
            
    # return caminho
    return f'{file_name}.nc'


# -----------------------------------------------------------------------------------------------------------
def download_prod(yyyymmddhhmn, product_name, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]
            # Alteração problema servidor
            new_file_name = file_name.replace('OR', 'CG')
            # print(file_name)
            # print(new_file_name)

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}{file_name}.nc exists')
            else:
                print(f'Downloading file {file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')

            # Alteração problema servidor
            old_file = os.path.join(path_dest, f'{file_name}.nc')
            new_file = os.path.join(path_dest, f'{new_file_name}.nc')
            os.rename(old_file, new_file)
    return f'{file_name}'

def download_DMW(yyyymmddhhmn,Ch,path_dest):
    product_name = 'ABI-L2-DMWF'
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6C{int(Ch):02.0f}_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}.nc'

def download_PROD(yyyydddhhmn, product_name, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%j')
    hour = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%H')
    min = datetime.strptime(yyyydddhhmn, '%Y%j%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyydddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}.nc'


def download_PROD_simposio(yyyydddhhmn, product_name, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = yyyydddhhmn[0:4]
    day_of_year = yyyydddhhmn[4:7]
    hour = yyyydddhhmn[7:9]
    min = yyyydddhhmn[9:]

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyydddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}.nc'


# -----------------------------------------------------------------------------------------------------------
def download_GLM(yyyymmddhhmnss, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%H')
    min = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%M')
    seg = datetime.strptime(yyyymmddhhmnss, '%Y%m%d%H%M%S').strftime('%S')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    product_name = "GLM-L2-LCFA"
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}_G16_s{year}{day_of_year}{hour}{min}{seg}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmnss}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}'

# -----------------------------------------------------------------------------------------------------------

# Functions to convert lat / lon extent to array indices
def geo2grid(lat, lon,nc):
    # Apply scale and offset
    xscale, xoffset = nc.variables['x'].scale_factor, nc.variables['x'].add_offset
    # print(xscale,xoffset)
    yscale, yoffset = nc.variables['y'].scale_factor, nc.variables['y'].add_offset
    # print(yscale,yoffset)

    x, y = latlon2xy(lat, lon)
    print(x,y)
    col = (x - xoffset) / xscale
    lin = (y - yoffset) / yscale
    return int(lin), int(col)

# -----------------------------------------------------------------------------------------------------------
def latlon2xy(lat, lon):
    # goes_imagery_projection:semi_major_axis
    req = 6378137  # meters
    #  goes_imagery_projection:inverse_flattening
    invf = 298.257222096
    # goes_imagery_projection:semi_minor_axis
    rpol = 6356752.31414  # meters
    e = 0.0818191910435
    # goes_imagery_projection:perspective_point_height + goes_imagery_projection:semi_major_axis
    H = 42164160  # meters
    # goes_imagery_projection: longitude_of_projection_origin
    lambda0 = -1.308996939

    # Convert to radians
    latRad = lat * (math.pi / 180)
    lonRad = lon * (math.pi / 180)

    # (1) geocentric latitude
    Phi_c = math.atan(((rpol * rpol) / (req * req)) * math.tan(latRad))
    # (2) geocentric distance to the point on the ellipsoid
    rc = rpol / (math.sqrt(1 - ((e * e) * (math.cos(Phi_c) * math.cos(Phi_c)))))
    # (3) sx
    sx = H - (rc * math.cos(Phi_c) * math.cos(lonRad - lambda0))
    # (4) sy
    sy = -rc * math.cos(Phi_c) * math.sin(lonRad - lambda0)
    # (5)
    sz = rc * math.sin(Phi_c)

    # x,y
    x = math.asin((-sy) / math.sqrt((sx * sx) + (sy * sy) + (sz * sz)))
    y = math.atan(sz / sx)

    return x, y

# -----------------------------------------------------------------------------------------------------------
# Function to convert lat / lon extent to GOES-16 extents
def convertExtent2GOESProjection(extent):
    # GOES-16 viewing point (satellite position) height above the earth
    GOES16_HEIGHT = 35786023.0
    # GOES-16 longitude position
    GOES16_LONGITUDE = -75.0

    a, b = latlon2xy(extent[1], extent[0])
    c, d = latlon2xy(extent[3], extent[2])
    return (a * GOES16_HEIGHT, c * GOES16_HEIGHT, b * GOES16_HEIGHT, d * GOES16_HEIGHT)

# -----------------------------------------------------------------------------------------------------------
def reprojectBruno(reproj_file, reproj_var, reproj_extent, reproj_resolution, path_dest):

    # test_file = reproj_file.split('\\')
    # test_file.reverse()
    # rtest_file = test_file[0].replace('.nc', f'_reproj_.nc')

    # if os.path.exists(f"{path_dest}{rtest_file}"):
    #     print('Arquivo já Reprojetado')
    #     caminho = f"{path_dest}{rtest_file}"
    #     print(caminho)
    #     return caminho
    
    
    def get_geot(ex, nlines, ncols):
        # Compute resolution based on data dimension
        resx = (ex[2] - ex[0]) / ncols
        resy = (ex[3] - ex[1]) / nlines
        return [ex[0], resx, 0, ex[3], 0, -resy]

    # GOES-16 Spatial Reference System
    source_prj = osr.SpatialReference()
    source_prj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 '
                               '+lat_0=0.0 +lon_0=-75 +sweep=x +no_defs')
    # Lat/lon WSG84 Spatial Reference System
    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Abrindo imagem com a biblioteca GDAL
    raw = gdal.Open(f'NETCDF:{reproj_file}:' + reproj_var, gdal.GA_ReadOnly)

    # Lendo os metadados do cabecalho
    if reproj_var == 'BCM' or 'Mask':  ### O arquivo Clear Sky não possui sacale e offset é um arquivo binário
        metadata = raw.GetMetadata()
        scale = 1
        offset = 0
        undef = float(metadata.get(reproj_var + '#_FillValue'))
        file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
        file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    else:  # Alteração João 26/08/2022
        metadata = raw.GetMetadata()
        scale = float(metadata.get(reproj_var + '#scale_factor'))
        offset = float(metadata.get(reproj_var + '#add_offset'))
        undef = float(metadata.get(reproj_var + '#_FillValue'))
        file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
        file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    # Setup projection and geo-transformation
    raw.SetProjection(source_prj.ExportToWkt())
    # raw.SetGeoTransform(raw.GetGeoTransform())
    GOES16_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]
    raw.SetGeoTransform(get_geot(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))

    # Compute grid dimension
    KM_PER_DEGREE = 111.32
    sizex = int(((reproj_extent[2] - reproj_extent[0]) * KM_PER_DEGREE) / reproj_resolution)
    sizey = int(((reproj_extent[3] - reproj_extent[1]) * KM_PER_DEGREE) / reproj_resolution)

    # Get memory driver
    driver = gdal.GetDriverByName('MEM')

    # Create grid
    grid = driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)

    # Setup projection and geo-transformation
    grid.SetProjection(target_prj.ExportToWkt())
    grid.SetGeoTransform(get_geot(reproj_extent, grid.RasterYSize, grid.RasterXSize))


    # Perform the projection/resampling
    gdal.ReprojectImage(raw, grid, source_prj.ExportToWkt(), target_prj.ExportToWkt(), gdal.GRA_NearestNeighbour,
                        options=['NUM_THREADS=ALL_CPUS'])
    
    
    # Close files
    raw = None
    del raw

    # Read grid data
    array = grid.ReadAsArray()

    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)

    # Aplicando scale, offset
    array = array * scale + offset


    grid.GetRasterBand(1).SetNoDataValue(-1)
    grid.GetRasterBand(1).WriteArray(array)

    # Define the parameters of the output file
    kwargs = {'format': 'netCDF',
              'dstSRS': target_prj,
              'outputBounds': (reproj_extent[0], reproj_extent[3], reproj_extent[2], reproj_extent[1]),
              'outputBoundsSRS': target_prj,
              'outputType': gdal.GDT_Float32,
              'srcNodata': undef,
              'dstNodata': 'nan',
              'resampleAlg': gdal.GRA_NearestNeighbour}

    reproj_file = reproj_file.split('\\')
    reproj_file.reverse()
    r_file = reproj_file[0].replace('.nc', f'_reproj_.nc')
    gdal.Warp(f'{path_dest}/{r_file}', grid, **kwargs)

    # os.remove(f'{input}{reproj_file[0]}')
    # return file_dtime, file_satellite, grid
    caminho = F'{path_dest}{r_file}'
    return caminho

# -----------------------------------------------------------------------------------------------------------
def Degrees(file_id):

    # proj_info = file_id.variables['goes_imager_projection']
    # lon_origin = proj_info.longitude_of_projection_origin
    # H = proj_info.perspective_point_height + proj_info.semi_major_axis
    # r_eq = proj_info.semi_major_axis
    # r_pol = proj_info.semi_minor_axis

    lon_origin = 0
    H = 35786023 + 6378137
    r_eq = 6378137
    r_pol = 6356752.31414

    # Data info
    lat_rad_1d = file_id.variables['x'][:]
    lon_rad_1d = file_id.variables['y'][:]

    # Create meshgrid filled with radian angles
    lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)

    # lat/lon calculus routine from satellite radian angle vectors
    lambda_0 = (lon_origin * np.pi) / 180.0

    a_var = np.power(np.sin(lat_rad), 2.0) + (np.power(np.cos(lat_rad), 2.0) * (
            np.power(np.cos(lon_rad), 2.0) + (((r_eq * r_eq) / (r_pol * r_pol)) * np.power(np.sin(lon_rad), 2.0))))
    b_var = -2.0 * H * np.cos(lat_rad) * np.cos(lon_rad)
    c_var = (H ** 2.0) - (r_eq ** 2.0)

    r_s = (-1.0 * b_var - np.sqrt((b_var ** 2) - (4.0 * a_var * c_var))) / (2.0 * a_var)
    

    s_x = r_s * np.cos(lat_rad) * np.cos(lon_rad)
    s_y = - r_s * np.sin(lat_rad)
    s_z = r_s * np.cos(lat_rad) * np.sin(lon_rad)

    Lat = (180.0 / np.pi) * (
        np.arctan(((r_eq * r_eq) / (r_pol * r_pol)) * ((s_z / np.sqrt(((H - s_x) * (H - s_x)) + (s_y * s_y))))))
    Lon = (lambda_0 - np.arctan(s_y / (H - s_x))) * (180.0 / np.pi)

    return Lat, Lon

# -----------------------------------------------------------------------------------------------------------
def save_txt(array, nome_arquivo, diretorio, type):
    # Type fdcf ou table
    # Checa se a matriz é vazia
    if type == 'fdcf' and len(array) == 0:
        return print(f'Matriz {nome_arquivo} Vazia')
    elif type == 'table' and np.sum(array) == 0:
        return print(f'Soma {nome_arquivo} Vazia')
    else:
        if type == 'fdcf':
            # criando nome do arquivo e diretório -- Mudar as barras para Linux
            with open(f'{diretorio}\\{nome_arquivo}.txt', 'w') as file:
                for valor in array:
                    valor = f'{valor[0]},{valor[1]}\n'
                    file.write(valor)
                file.close()
        else:
            with open(f'{diretorio}\\{nome_arquivo}.txt', 'w') as file:
                for valor in array:
                    valor = f'{int(valor)}\n'
                    file.write(valor)
                file.close()

        # Para ler os arquivos salvos com type fdcf:
        # array = np.loadtxt(diretory arquivo, 'float', delimiter=',')

        # Para ler os arquivos salvo com type table
        # data = []
        # data_table = open(f'{dir_out}txt\\fdcf_total_UF_20220719_br.txt', 'r')
        # for valor in data_table:
        #     valor = valor.split(',')
        #     sigla = valor[0]
        #     num = valor[1].replace('\n', '')
        #     data.append((sigla, num))

# -----------------------------------------------------------------------------------------------------------
def save_log_erro(array_errors, nome_arquivo, diretorio):
    if len(array_errors) == 0:
        pass
    else:
        # criando nome do arquivo e diretório -- Mudar as barras para Linux
        with open(f'{diretorio}\\{nome_arquivo}.txt', 'w') as file:
            for valor in array_errors:
                erro = f'{valor}\n'
                file.write(erro)

# -----------------------------------------------------------------------------------------------------------
def read_fdcf(directory, nome_arquivo):
    with open(f'{directory}\\{nome_arquivo}', 'r') as file:
        for linha in file.readlines():
            print(linha, end='')


# -----------------------------------------------------------------------------------------------------------
# Function to reproject the yyyymmddhhmn
def reproject(file_name, ncfile, array, extent, undef):
    # Read the original file projection and configure the output_Plot projection
    source_prj = osr.SpatialReference()
    source_prj.ImportFromProj4(ncfile.GetProjectionRef())

    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

    # Reproject the yyyymmddhhmn
    GeoT = ncfile.GetGeoTransform()
    driver = gdal.GetDriverByName('MEM')
    raw = driver.Create('raw', array.shape[0], array.shape[1], 1, gdal.GDT_Float32)
    raw.SetGeoTransform(GeoT)
    raw.GetRasterBand(1).WriteArray(array)

    # Define the parameters of the output_Plot file
    kwargs = {'format': 'netCDF',
              'srcSRS': source_prj,
              'dstSRS': target_prj,
              'outputBounds': (extent[0], extent[3], extent[2], extent[1]),
              'outputBoundsSRS': target_prj,
              'outputType': gdal.GDT_Float32,
              'srcNodata': undef,
              'dstNodata': 'nan',
              'resampleAlg': gdal.GRA_NearestNeighbour}
    # Write the reprojected file on disk
    gdal.Warp(file_name, raw, **kwargs)

# -----------------------------------------------------------------------------------------------------------
def reproject_CLSM(reproj_file, reproj_var, reproj_extent):
    global dir_in

    def get_geot(ex, nlines, ncols):
        # Compute resolution based on data dimension
        resx = (ex[2] - ex[0]) / ncols
        resy = (ex[3] - ex[1]) / nlines
        return [ex[0], resx, 0, ex[3], 0, -resy]

    def getScaleOffset(path, variable):
        nc = Dataset(path, mode='r')

        if (variable == "BCM") or (variable == "Phase") or (variable == "Smoke") or (variable == "Dust") or (
                variable == "Mask") or (variable == "Power"):
            scale = 1
            offset = 0
        else:
            scale = nc.variables[variable].scale_factor
            offset = nc.variables[variable].add_offset
        # scale = 0
        # offset = 0
        nc.close()
        return scale, offset


    # Open the file
    img = gdal.Open(f'NETCDF:{reproj_file}:' + reproj_var)

    # Data Quality Flag (DQF)
    dqf = gdal.Open(f'NETCDF:{reproj_file}:DQF')

    # Read the header metadata
    metadata = img.GetMetadata()
    undef = float(metadata.get(reproj_var + '#_FillValue'))
    dtime = metadata.get('NC_GLOBAL#time_coverage_start')

    # Load the data
    ds = (img.ReadAsArray(0, 0, img.RasterXSize, img.RasterYSize).astype(float))
    ds_dqf = dqf.ReadAsArray(0, 0, dqf.RasterXSize, dqf.RasterYSize).astype(float)

    # Apply NaN's where the quality flag is greater than 1
    ds[ds_dqf > 1] = np.nan

    # GOES-16 Spatial Reference System
    source_prj = osr.SpatialReference()
    source_prj.ImportFromProj4(
        '+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs')
    # Lat/lon WSG84 Spatial Reference System
    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Reprojetando
    GOES16_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]
    GeoT = get_geot(GOES16_EXTENT, img.RasterYSize, img.RasterXSize)
    driver = gdal.GetDriverByName('MEM')
    raw = driver.Create('raw', ds.shape[0], ds.shape[1], 1, gdal.GDT_Float32)
    raw.SetGeoTransform(GeoT)
    raw.GetRasterBand(1).WriteArray(ds)

    # Define the parameters of the output_Plot file
    kwargs = {'format': 'netCDF',
              'srcSRS': source_prj,
              'dstSRS': target_prj,
              'outputBounds': (reproj_extent[0], reproj_extent[3], reproj_extent[2], reproj_extent[1]),
              'outputBoundsSRS': target_prj,
              'outputType': gdal.GDT_Float32,
              'srcNodata': undef,
              'dstNodata': 'nan',
              'resampleAlg': gdal.GRA_NearestNeighbour}

    reproj_file = reproj_file.split('\\')
    reproj_file.reverse()
    r_file = reproj_file[0].replace('.nc', f'_reproj_.nc')
    gdal.Warp(f'{dir_in}/{r_file}', raw, **kwargs)

    # os.remove(f'{input}{reproj_file[0]}')
    # return file_dtime, file_satellite, grid
    caminho = F'{dir_in}{r_file}'
    return caminho

# -----------------------------------------------------------------------------------------------------------
def salva_diario_txt_fdcf(dir_txt, yyyymmddhhmn_ini, hours_loop):
    # Initial time and date
    yyyy = datetime.strptime(yyyymmddhhmn_ini, '%Y%m%d%H%M').strftime('%Y')
    mm = datetime.strptime(yyyymmddhhmn_ini, '%Y%m%d%H%M').strftime('%m')
    dd = datetime.strptime(yyyymmddhhmn_ini, '%Y%m%d%H%M').strftime('%d')
    hh = datetime.strptime(yyyymmddhhmn_ini, '%Y%m%d%H%M').strftime('%H')
    mn = datetime.strptime(yyyymmddhhmn_ini, '%Y%m%d%H%M').strftime('%M')

    date_ini = str(datetime(int(yyyy), int(mm), int(dd), int(hh), int(mn)))
    date_end = str(datetime(int(yyyy), int(mm), int(dd), int(hh), int(mn)) + timedelta(hours=hours_loop))

    array_diario = []
    while date_ini < date_end:
        file_id = f'fdcf_{date_ini}_br'
        array = np.loadtxt(f'{dir_txt}{file_id}.txt', 'float', delimiter=',')
        for p in array:
            array_diario.append(p)
        # Incrementa +10 min
        date_ini = str(datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S') + timedelta(minutes=10))

    file_id_diario = f'fdcf_{yyyy}{dd}_br'
    save_txt(array_diario, file_id_diario, dir_txt)

# -----------------------------------------------------------------------------------------------------------
def imprime_matriz(M):
    for linha in M:
        linha_string = [str(elem) for elem in linha]
        print("\t".join(linha_string))

# -----------------------------------------------------------------------------------------------------------
def save_array(array,name_file,diretorio):
    with open(f'{diretorio}{name_file}.txt','w+') as file:
        for linha in range(array.shape[0]):
            linha_string = [str(elem) for elem in array[linha]]
            file.write(str(linha_string))
        file.close()

# -----------------------------------------------------------------------------------------------------------
def reprojectDownload(reproj_file, reproj_var, reproj_extent, reproj_resolution, dir_out):
    def get_geot(ex, nlines, ncols):
        # Compute resolution based on data dimension
        resx = (ex[2] - ex[0]) / ncols
        resy = (ex[3] - ex[1]) / nlines
        return [ex[0], resx, 0, ex[3], 0, -resy]

    # GOES-16 Spatial Reference System
    source_prj = osr.SpatialReference()
    source_prj.ImportFromProj4(
        '+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs')
    # Lat/lon WSG84 Spatial Reference System
    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Abrindo imagem com a biblioteca GDAL
    raw = gdal.Open(f'NETCDF:{reproj_file}:' + reproj_var, gdal.GA_ReadOnly)

    # Lendo os metadados do cabecalho
    if reproj_var == 'BCM':  ### O arquivo Clear Sky não possui sacale e offset é um arquivo binário
        metadata = raw.GetMetadata()
        scale = 1
        offset = 0
        undef = float(metadata.get(reproj_var + '#_FillValue'))
        file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
        file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    else:  # Alteração João 26/08/2022
        metadata = raw.GetMetadata()
        scale = float(metadata.get(reproj_var + '#scale_factor'))
        offset = float(metadata.get(reproj_var + '#add_offset'))
        undef = float(metadata.get(reproj_var + '#_FillValue'))
        file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
        file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    # Setup projection and geo-transformation
    raw.SetProjection(source_prj.ExportToWkt())
    # raw.SetGeoTransform(raw.GetGeoTransform())
    GOES16_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]
    raw.SetGeoTransform(get_geot(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))

    # Compute grid dimension
    KM_PER_DEGREE = 111.32
    sizex = int(((reproj_extent[2] - reproj_extent[0]) * KM_PER_DEGREE) / reproj_resolution)
    sizey = int(((reproj_extent[3] - reproj_extent[1]) * KM_PER_DEGREE) / reproj_resolution)

    # Get memory driver
    driver = gdal.GetDriverByName('MEM')

    # Create grid
    grid = driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)

    # Setup projection and geo-transformation for img
    grid.SetProjection(target_prj.ExportToWkt())
    grid.SetGeoTransform(get_geot(reproj_extent, grid.RasterYSize, grid.RasterXSize))

    # Perform the projection/resampling for img
    gdal.ReprojectImage(raw, grid, source_prj.ExportToWkt(), target_prj.ExportToWkt(), gdal.GRA_NearestNeighbour,
                        options=['NUM_THREADS=ALL_CPUS'])

    # Close files
    raw = None
    del raw

    # Read grid data
    array = grid.ReadAsArray()

    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)

    # Aplicando scale, offset
    array = array * scale + offset

    # array[DQF > 1] = np.nan
    grid.GetRasterBand(1).SetNoDataValue(-1)
    grid.GetRasterBand(1).WriteArray(array)

    # Define the parameters of the output file
    kwargs = {'format': 'netCDF',
              'dstSRS': target_prj,
              'outputBounds': (reproj_extent[0], reproj_extent[3], reproj_extent[2], reproj_extent[1]),
              'outputBoundsSRS': target_prj,
              'outputType': gdal.GDT_Float32,
              'srcNodata': undef,
              'dstNodata': 'nan',
              'resampleAlg': gdal.GRA_NearestNeighbour}

    os.remove(reproj_file)
    reproj_file = reproj_file.split('\\')
    reproj_file.reverse()
    r_file = reproj_file[0].replace('.nc', f'_reproj_.nc')
    gdal.Warp(f'{dir_out}/{r_file}', grid, **kwargs)

    # os.remove(f'{input}{reproj_file[0]}')
    # return file_dtime, file_satellite, grid
    caminho = F'{dir_out}{r_file}'

    return caminho

# -----------------------------------------------------------------------------------------------------------
def plot_NDVI(array_ndvi,extent,name_img,dir_out):
    global dir_ndvi,dir_logos,dir_shapefiles
    # Choose the plot size (width x height, in inches)
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black',
                     facecolor='black')

    # Use the Geostationary projection in cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

    # Estados, coastline e shape Brasil
    estados = shpreader.Reader(f'{dir_shapefiles}gadm40_BRA_shp\gadm40_BRA_1.shp').geometries()
    ax.add_geometries(estados, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.7, zorder=4)
    ax.coastlines(resolution='10m', color='black', linewidth=0.5, zorder=4)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5, zorder=4)
    # # Add an ocean mask
    ocean = ax.add_feature(cfeature.OCEAN, facecolor='lightsteelblue', zorder=3)
    gl = ax.gridlines(color='white', linestyle='--', linewidth=0.25, alpha=0.7, xlocs=np.arange(-180, 180, 5),
                      ylocs=np.arange(-90, 90, 5), zorder=5)
    gl.top_labels = False
    gl.right_labels = False

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat

    # Custom Colormap
    colors = ["#0000FF", "#0030FF", "#0060FF", "#0090FF", "#00C0FF", "#003000", "#164A16", "#00C000", "#FFFF00",
              "#FF8000", "#FF0000"]
    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)  # Create a custom linear colormap

    img = ax.imshow(array_ndvi, origin='upper', vmin=-1, vmax=1, cmap=my_cmap, extent=img_extent)

    cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
    cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[-0.66, -0.33, 0, 0.33, 0.66])
    cb.ax.set_xticklabels(['-0.66', '-0.33', '0', '0.33', '0.66'])
    cb.ax.tick_params(axis='x', colors='black',
                      labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de cores

    # # Adicionando descricao da imagem
    # date = datetime.strptime(name_img, '%Y%m%d').strftime('%d-%b-%Y')
    description = f"Composição NDVI satélite GOES-16 {name_img}"
    institution = f'CEPAGRI - UNICAMP'
    # Criando novos eixos de acordo com a posicao da imagem
    cax1 = fig.add_axes(
        [ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
    cax1.patch.set_color('black')  # Alterando a cor do novo eixo
    cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
    cax1.text(0.85, 0.13, institution, color='yellow', size=10)  # Adicionando texto
    cax1.xaxis.set_visible(False)  # Removendo rotulos do eixo X
    cax1.yaxis.set_visible(False)  # Removendo rotulos do eixo Y

    # Adicionando os logos
    logo_noaa = plt.imread(dir_logos + 'NOAA_Logo.png')  # Lendo o arquivo do logo
    logo_goes = plt.imread(dir_logos + 'GOES_Logo.png')  # Lendo o arquivo do logo
    logo_cepagri = plt.imread(dir_logos + 'CEPAGRI-Logo.png')  # Lendo o arquivo do logo
    fig.figimage(logo_noaa, 32, 233, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_goes, 10, 150, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_cepagri, 10, 70, zorder=3, alpha=0.8, origin='upper')  # Plotando logo

    # Salvando a imagem de saída
    plt.savefig(f'{dir_out}NDVI_{name_img}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
    plt.close()

# -----------------------------------------------------------------------------------------------------------
def remap(path, variable, extent, resolution):
    
    def getGeoT(extent, nlines, ncols):
        # Compute resolution based on data dimension
        resx = (extent[2] - extent[0]) / ncols
        resy = (extent[3] - extent[1]) / nlines
        return [extent[0], resx, 0, extent[3] , 0, -resy]

    def getScaleOffset(path, variable):
        nc = Dataset(path, mode='r')
        
        if (variable == "BCM") or (variable == "Phase") or (variable == "Smoke") or (variable == "Dust") or (variable == "Mask") or (variable == "Power"): 
            scale  = 1
            offset = 0     
        else:
            scale = nc.variables[variable].scale_factor
            offset = nc.variables[variable].add_offset
        nc.close()
            
        return scale, offset
    
    # Default scale    
    scale = 1
    
    # Default offset
    offset = 0
    
    # GOES Extent (satellite projection) [llx, lly, urx, ury]
    GOES_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]
    
    # Setup NetCDF driver
    gdal.SetConfigOption('GDAL_NETCDF_BOTTOMUP', 'NO')
        
    if not (variable == "DQF"):              
        # Read scale/offset from file
        scale, offset = getScaleOffset(path, variable) 
      
    connectionInfo = f'NETCDF:\"' + path + '\"://' + variable
	
    # Lendo os metadados do cabecalho
    # raw = gdal.Open(connectionInfo, gdal.GA_ReadOnly) 
    raw = gdal.Open(f'NETCDF:{path}:' + variable, gdal.GA_ReadOnly)        
    
    # GOES-16 Spatial Reference System
    sourcePrj = osr.SpatialReference()
    sourcePrj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs')
    # Lat/lon WSG84 Spatial Reference System
    targetPrj = osr.SpatialReference()
    targetPrj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Setup projection and geo-transformation
    raw.SetProjection(sourcePrj.ExportToWkt())
    raw.SetGeoTransform(getGeoT(GOES_EXTENT, raw.RasterYSize, raw.RasterXSize))  

    # Compute grid dimension
    KM_PER_DEGREE = 111.32
    sizex = abs(int(((extent[2] - extent[0]) * KM_PER_DEGREE) / resolution))
    sizey = abs(int(((extent[3] - extent[1]) * KM_PER_DEGREE) / resolution))

    # print(sizex,sizey)
    
    # Get memory driver
    memDriver = gdal.GetDriverByName('MEM')
   
    # Create grid
    grid = memDriver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)
        
    # Setup projection and geo-transformation
    grid.SetProjection(targetPrj.ExportToWkt())
    grid.SetGeoTransform(getGeoT(extent, grid.RasterYSize, grid.RasterXSize))
    
    gdal.ReprojectImage(raw, grid, sourcePrj.ExportToWkt(), targetPrj.ExportToWkt(), gdal.GRA_NearestNeighbour, options=['NUM_THREADS=ALL_CPUS']) 

    # Read grid data
    array = grid.ReadAsArray()
    
    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)
    
    # Read as uint16
    array = array.astype(np.uint16)  
       
    # Apply scale and offset
    array = array * scale + offset

    # Get the raster 
    grid.GetRasterBand(1).WriteArray(array)

	# Close file
    raw = None
    del raw
	
    return grid

# -----------------------------------------------------------------------------------------------------------
def get_lat_lon_fromgrid(file_fdcf,extent):
    
    #Lendo os arquivos do cabeçalho do arquivo reprojetado
    data = Dataset(file_fdcf).variables['Band1'][:]
    lat,lon = data.shape

    #Calculando a escala das latitudes e longitudes
    lat_scale_factor = (extent[2]-extent[0])/lat
    lon_scale_factor = (extent[3]-extent[1])/lon
    # print(lat_scale_factor,lon_scale_factor)

    #Valores inciais de latitude e longitude
    lat_0 = extent[2]
    lon_0 = extent[3]

    #Informações fixas do satélite GOES
    lon_origin = -75
    H = 35786023 + 6378137
    r_eq =  6378137
    r_pol = 6356752.31414

    # Criando as matrizes com latitudes 
    lat_rad_1d = np.zeros(lat)
    lat_rad_1d[0]= lat_0
    for i in range(1,lat):
        lat_rad_1d[i] = lat_rad_1d[i-1]-lat_scale_factor


    # print(lat_rad_1d)

    # Matriz das longitudes
    lon_rad_1d = np.zeros(lon)
    lon_rad_1d[0]=lon_0
    for i in range(1,lon):
        lon_rad_1d[i] = lon_rad_1d[i-1]-lon_scale_factor
        
    # print(lon_rad_1d)

    # Create meshgrid filled with radian angles
    lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)

    # print(lat_rad)

    # lat/lon calculus routine from satellite radian angle vectors
    lambda_0 = (lon_origin * np.pi) / 180.0

    a_var = np.power(np.sin(lat_rad), 2.0) + (np.power(np.cos(lat_rad), 2.0) * (
            np.power(np.cos(lon_rad), 2.0) + (((r_eq * r_eq) / (r_pol * r_pol)) * np.power(np.sin(lon_rad), 2.0))))
    b_var = -2.0 * H * np.cos(lat_rad) * np.cos(lon_rad)
    c_var = (H ** 2.0) - (r_eq ** 2.0)

    r_s = (-1.0 * b_var - np.sqrt((b_var ** 2) - (4.0 * a_var * c_var))) / (2.0 * a_var)
    # print(f'r_s \n{r_s}')

    s_x = r_s * np.cos(lat_rad) * np.cos(lon_rad)
    s_y = - r_s * np.sin(lat_rad)
    s_z = r_s * np.cos(lat_rad) * np.sin(lon_rad)

    Lat = (180.0 / np.pi) * (
        np.arctan(((r_eq * r_eq) / (r_pol * r_pol)) * ((s_z / np.sqrt(((H - s_x) * (H - s_x)) + (s_y * s_y))))))
    Lon = (lambda_0 - np.arctan(s_y / (H - s_x))) * (180.0 / np.pi)

    # print(f's_x{s_x}')
    # print(f's_y{s_y}')
    # print(f's_z{s_z}')
    # print(f'r_eq{r_eq}')
    # print(f'r_pol{r_pol}')


    return Lat, Lon

###-----------------Funções Novo Processamento.py ------ By Guilhemre Moura--------------------------------------
###----------------- git: https://github.com/guimouraO1 ---------------------------------------------------------


def area_para_recorte(v_extent):
    # Area de interesse para recorte
    if v_extent == 'br':
        # Brasil # Esquerda , baixo,  direita, cima
        extent = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 2
    elif v_extent == 'sp':
        # São Paulo
        extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 1.0
    elif v_extent == 'pant':
        # São Paulo
        extent = [-59.3, -22.3, -54.7, -15.5] # Min lon, Min lat, Max lon, Max lat - Pantanal
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 2.0
    else:
        extent = [-115.98, -55.98, -25.01, 34.98]  # Min lon, Min lat, Max lon, Max lat
        resolution = 2.0
    return extent, resolution
###-----------------Funções Novo Processamento.py ------ By Guilhemre Moura--------------------------------------
###----------------- git: https://github.com/guimouraO1 ---------------------------------------------------------

def area_para_recorte(v_extent):
    # Area de interesse para recorte
    if v_extent == 'br':
        # Brasil # Esquerda , baixo,  direita, cima
        extent = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 4.0
    elif v_extent == 'sp':
        # São Paulo
        extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 1.0
    elif v_extent == 'pant':
        # São Paulo
        extent = [-59.3, -22.3, -54.7, -15.5]  # Min lon, Min lat, Max lon, Max lat
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 2
    else:
        extent = [-115.98, -55.98, -25.01, 34.98]  # Min lon, Min lat, Max lon, Max lat
        resolution = 2.0
    return extent, resolution

def adicionando_shapefile(v_extent, ax, colors='cyan'):
    if v_extent == 'br':
        # Adicionando o shapefile dos estados brasileiros
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil\\BR_UF_2020.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=0.7, zorder=8)
    elif v_extent == 'sp':
        # Adicionando o shapefile dos estados brasileiros e cidade de campinas
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil\\BR_UF_2020.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=0.7, zorder=8)
        shapefile = list(shpreader.Reader(dir_shapefiles + 'campinas\\campinas.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=1, zorder=8)

    elif v_extent == 'pant':
        # shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil\\BR_UF_2020.shp').geometries())
        # ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=0.7, zorder=8)
        biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
        pantanal = biomas[5]   
        ax.add_geometries(pantanal, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=0.7, zorder=8)

def adicionando_linhas(ax, colors='cyan'):
    # Adicionando  linhas dos litorais
    ax.coastlines(resolution='10m', color=colors, linewidth=0.5, zorder=8)
    # Adicionando  linhas das fronteiras
    ax.add_feature(cartopy.feature.BORDERS, edgecolor=colors, linewidth=0.5, zorder=8)
    # Adicionando  paralelos e meridianos
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=0.7, linestyle='--', linewidth=0.2, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), zorder=8)
    gl.top_labels = False
    gl.right_labels = False

def adicionando_descricao_imagem(description, institution, ax, fig, cruz=False):
    cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
    cax1.patch.set_color('black')  # Alterando a cor do novo eixo
    cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
    if cruz:
        cruzStr = '+'
        cax1.text(0.192, 0.13, cruzStr, color='red', size=12)  # Adicionando símbolo "+"
    cax1.text(0.85, 0.13, institution, color='yellow', size=10)  # Adicionando texto
    cax1.xaxis.set_visible(False)  # Removendo rótulos do eixo X
    cax1.yaxis.set_visible(False)  # Removendo rótulos do eixo Y

def adicionando_logos(fig):
    # Adicionando os logos
    logo_noaa = plt.imread(dir_logos + 'NOAA_Logo.png')  # Lendo o arquivo do logo
    logo_goes = plt.imread(dir_logos + 'GOES_Logo.png')  # Lendo o arquivo do logo
    logo_cepagri = plt.imread(dir_logos + 'CEPAGRI-Logo.png')  # Lendo o arquivo do logo
    fig.figimage(logo_noaa, 32, 240, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_goes, 10, 160, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_cepagri, 10, 80, zorder=3, alpha=0.8, origin='upper')  # Plotando logo


def apply_cira_stretch(band_data):
    
    log_root = np.log10(0.0223)
    denom = (1.0 - log_root) * 0.75
    band_data *= 0.01
    band_data = band_data.clip(np.finfo(float).eps)
    band_data = np.log10(band_data)
    band_data -= log_root
    band_data /= denom
    return 1 + band_data

def applying_rayleigh_correction(file_ch01, utc_time, lons, lats, sun_zenith, data_ch01, data_ch02, longitude):
    # Altitude do satélite
    sat_h = file_ch01.variables['goes_imager_projection'].perspective_point_height

    sunalt, suna = get_alt_az(utc_time, lons, lats)
    suna = np.rad2deg(suna)
    #sata, satel = get_observer_look(sat_lon, sat_lat, sat_alt, vis.attrs['start_time'], lons, lats, 0)
    sata, satel = get_observer_look(longitude, 0.0, sat_h, utc_time, lons, lats, 0)
    satz = 90 - satel

    # Correção de Rayleigh
    atmosphere = 'us-standard'
    aerosol_type = 'rayleigh_only'
    corrector = Rayleigh('GOES-16', 'abi', atmosphere=atmosphere, aerosol_type=aerosol_type)

    sata = sata % 360.
    suna = suna % 360.
    ssadiff = np.absolute(suna - sata)
    ssadiff = np.minimum(ssadiff, 360 - ssadiff)

    red = data_ch02 * 100

    refl_cor_band_c01 = corrector.get_reflectance(sun_zenith, satz, ssadiff, 'C01', redband=red)
    data_ch01 = data_ch01 - (refl_cor_band_c01 / 100)

    refl_cor_band_c02 = corrector.get_reflectance(sun_zenith, satz, ssadiff, 'C02', redband=red)
    data_ch02 = data_ch02 - (refl_cor_band_c02 / 100)

    return data_ch01, data_ch02 

def calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03):
        
    year = date.strftime('%Y')
    month = date.strftime('%m')
    day = date.strftime('%d')
    hour = date.strftime('%H')
    minutes = date.strftime('%M')
    
    # Criar as latitudes e longitudes com base na extensão
    lat = np.linspace(extent[3], extent[1], data_ch01.shape[0])
    lon = np.linspace(extent[0], extent[2], data_ch01.shape[1])
    xx,yy = np.meshgrid(lon,lat)
    lons = xx.reshape(data_ch01.shape[0], data_ch01.shape[1])
    lats = yy.reshape(data_ch01.shape[0], data_ch01.shape[1])

    # Obter o ano, mês, dia, hora e minuto para aplicar a correção zenital
    utc_time = datetime(int(year), int(month), int(day), int(hour), int(minutes))
    sun_zenith = np.zeros((data_ch01.shape[0], data_ch01.shape[1]))
    sun_zenith = astronomy.sun_zenith_angle(utc_time, lons, lats)

    # Aplicar a correção zenital do sol
    data_ch01 = (data_ch01)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch02 = (data_ch02)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch03 = (data_ch03)/(np.cos(np.deg2rad(sun_zenith)))

    return utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03


def process_truecolor(v_extent, yyyymmddhhmn):
    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)
    # Variable to remap
    variable = "CMI"
    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------

    # Faz o download dos arquivos das bandas
    ch01 = download_CMI(yyyymmddhhmn,1,f'{dir_in}band{str(1).zfill(2)}')
    file_ch01 = Dataset(f'{dir_in}band01\\{ch01}')
    ch02 = download_CMI(yyyymmddhhmn,2,f'{dir_in}band{str(2).zfill(2)}')
    file_ch02 = Dataset(f'{dir_in}band02\\{ch02}')
    ch03 = download_CMI(yyyymmddhhmn,3,f'{dir_in}band{str(3).zfill(2)}')
    file_ch03 = Dataset(f'{dir_in}band03\\{ch03}')
    ch13 = download_CMI(yyyymmddhhmn,13,f'{dir_in}band{str(13).zfill(2)}')
    file_ch13 = Dataset(f'{dir_in}band13\\{ch13}')

    # Lê a longitude central
    longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin
    # Lê a data do arquivo
    add_seconds = int(file_ch01.variables['time_bounds'][0])
    date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
    date_file = date.strftime('%Y%m%d_%H%M%S')
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')

    #------------------------------------------------------------------------------------------------------#
    #-------------------------------------------Reprojetando----------------------------------------------#
    #------------------------------------------------------------------------------------------------------#
    # reprojetando band 01
    grid = remap(f'{dir_in}band{str(1).zfill(2)}\\{ch01}', variable, extent, resolution)
    # Lê o retorno da função
    data_ch01 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # reprojetando band 02
    grid = remap(f'{dir_in}band{str(2).zfill(2)}\\{ch02}', variable, extent, resolution)
    # Lê o retorno da função
    data_ch02 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # reprojetando band 03
    grid = remap(f'{dir_in}band{str(3).zfill(2)}\\{ch03}', variable, extent, resolution)
    # Lê o retorno da função 
    data_ch03 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
     # reprojetando band 13
    grid = remap(f'{dir_in}band{str(13).zfill(2)}\\{ch13}', variable, extent, resolution)
    # Lê o retorno da função
    data_ch13 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    
    # Calculando correção zenith
    utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03 = calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03)

    # Aplicando a correção de Rayleigh
    data_ch01, data_ch02 = applying_rayleigh_correction(file_ch01, utc_time, lons, lats, sun_zenith, data_ch01, data_ch02, longitude)

    # Calculando as cores verdadeiras (True color)
    R = data_ch02
    G = (data_ch01 + data_ch02) / 2 * 0.93 + 0.07 * data_ch03 
    B = data_ch01
    # Aplicando o estiramento CIRA
    R = apply_cira_stretch(R)
    G = apply_cira_stretch(G)
    B = apply_cira_stretch(B)

    # Create the RGB
    RGB = np.stack([R, G, B], axis=2)		

    # If zenith angle is greater than 85°, the composite pixel is zero
    RGB[sun_zenith > 85] = 0
    # Create the mask for the regions with zero
    mask = (RGB == [0.0,0.0,0.0]).all(axis=2)
    # Apply the mask to overwrite the pixels
    RGB[mask] = [0,0,0]
    # Create the fading transparency between the regions with the sun zenith angle of 75° and 85°
    alphas = sun_zenith / 100
    min_sun_angle = 0.75
    max_sun_angle = 0.85
    # Normalize the transparency mask
    alphas = ((alphas - max_sun_angle) / (min_sun_angle - max_sun_angle))
    RGB = np.dstack((RGB, alphas))
    
    return RGB

   