# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------
# Required modules
import time as t
from netCDF4 import Dataset  # Read / Write NetCDF4 files
from datetime import timedelta, datetime  # Basic Dates and time types
import numpy as np  # Scientific computing with Python
from osgeo import gdal  # Python bindings for GDAL
import os  # Utilitario para trabalhar com chamadas de sistema
import re  # Utilitario para trabalhar com expressoes regulares
# -----------------------------------------------------------------------------------------------------------

start = t.time()
print('Script Iniciado...')
gdal.PushErrorHandler('CPLQuietErrorHandler')  # Ignore GDAL warnings

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
extent_br = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat - Brasil
extent_pantanal = [-59.3, -22.3, -54.7, -15.5] # Min lon, Min lat, Max lon, Max lat - Pantanal

## Dias Julianos 
ano = 2020
mes = 'Jul'
start_jul = 183
end_jul = 213


#Cria o arquivo de teste Para processamento do NDVI- Mensal:
NDVI_maxmonth = np.zeros((378,256))
# Setar os valores todos para NAN pois max(nan,value) = value 
NDVI_maxmonth[NDVI_maxmonth == 0]= np.nan


for d in range(start_jul,end_jul+1):
    dia = f"{str(d).zfill(3)}"
    if (d==187)or (d==189):
        continue
    print(f'dia {dia}')
    try:
        dir_ndvi = f'{dir_in}ndvi\\{ano}\\{mes}\\{dia}\\';os.makedirs(dir_ndvi, exist_ok=True)
        
        #Cria uma lista como os arquivos de NDVI sem nuvens do respectivo dia
        ndvi_list = [name for name in os.listdir(dir_ndvi) if os.path.isfile(os.path.join(dir_ndvi, name)) and
                     re.match('^ndvi_%s.+_pant_max.npy$'%ano, name)]
        ndvi_list.sort()
        

        #Date to save de max day file
        date_file = ndvi_list[0][ndvi_list[0].find("ndvi_") + 5:ndvi_list[0].find("ndvi_") + 13]

        for f in ndvi_list:
            print(f)
            ndvi_file = np.load(f'{dir_ndvi}{f}', allow_pickle=True)

            NDVI_maxmonth = np.fmax(NDVI_maxmonth, ndvi_file)

    except Exception as erro:
        print(f'{erro} dia{dia}')
        continue

#Salva o acumulado Mensal

dir_month = f'{dir_in}ndvi\\{ano}\\{mes}\\'
name_file = f'ndvi_{mes}_pant_max.npy'
NDVI_maxmonth.dump(f'{dir_month}{name_file}')

print(f'Tempo total {round((t.time() - start), 2)} segundos ')