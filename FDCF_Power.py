# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------
# Required modules
import matplotlib.pyplot as plt  # Plotting library
from osgeo import gdal  # Python bindings for GDAL

print('Script Iniciado...')
gdal.PushErrorHandler('CPLQuietErrorHandler')  # Ignore GDAL warnings

#======================================================================================================
# Input and output directories
dir_in = "G:\Meu Drive\Minhas_coisas\Faculdade\ProjetoSAE\GOES\clone_dir_servidor\goes\\"
dir_main = "G:\Meu Drive\Minhas_coisas\Faculdade\ProjetoSAE\GOES\clone_dir_servidor\Scripts\goes\\"
dir_out = dir_main + "output\\"
dir_libs = dir_main + "libs\\"
dir_shapefiles = dir_main + "shapefiles\\"
dir_colortables = dir_main + "colortables\\"
dir_logos = dir_main + "logos\\"
# # Desired extent
# extent = [-90.0, -40.0, -20.0, 10.0]  # Max Lat, Max lon, Min lat, Min Lon Brasil
extent = [-74.0, -10, -55.5, 2.5]  #Mato Grosso


# Datetime to process
# yyyymmddhhmn = '202209100300'
yyyymmddhhmn = '202308311600'      
product_name = "ABI-L2-FDCF"