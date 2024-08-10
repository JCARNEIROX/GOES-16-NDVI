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


# Variaveis quando for Jun de 2024
ano = '2024'
start_jul = 153
end_jul = 182
mes = 'Jun'

# # Variaveis quando for Jul de 2024
# ano = '2024'
# start_jul = 183
# end_jul = 213
# mes = 'Jul'

print('Processando dados diários da Região1')
reg = 'reg1'

# Processa NDVI de cada dia

for d in range(start_jul,end_jul+1):
    dia = f"{str(d).zfill(3)}"
    print(f'Dia {dia}')

    dir_ndvi = f'{dir_in}ndvi\\{ano}\\{mes}\\{dia}\\' ;os.makedirs(dir_ndvi, exist_ok=True)
    
    #Cria uma lista como os arquivos de NDVI sem nuvens do respectivo dia
    ndvi_list = [name for name in os.listdir(dir_ndvi) if os.path.isfile(os.path.join(dir_ndvi, name)) and
                    re.match('^ndvi_%s.+_%s.npy$'%ano,reg, name)]
    ndvi_list.sort()
    # print(ndvi_list)
    teste = np.load(f'{dir_ndvi}{ndvi_list[0]}', allow_pickle=True)
    NDVI_maxday = np.zeros(teste.shape)

    # Setar os valores todos para NAN pois max(nan,value) = value 
    NDVI_maxday[NDVI_maxday == 0]= np.nan
    del teste

    #Date to save de max day file
    date_file = ndvi_list[0][ndvi_list[0].find("ndvi_") + 5:ndvi_list[0].find("ndvi_") + 13]

    for f in ndvi_list:
        print(f)
        try:
            ndvi_file = np.load(f'{dir_ndvi}{f}', allow_pickle=True)

            NDVI_maxday = np.fmax(NDVI_maxday, ndvi_file)
        except Exception as erro:
            print(f'{erro} dia{dia}')
            continue

    NDVI_maxday.dump(f'{dir_ndvi}ndvi_{date_file}_{reg}_max.npy')

    

print('Processando dados diários da Região2')
reg = 'reg2'

# Processa NDVI de cada dia

for d in range(start_jul,end_jul+1):
    dia = f"{str(d).zfill(3)}"
    print(f'Dia {dia}')

    dir_ndvi = f'{dir_in}ndvi\\{ano}\\{mes}\\{dia}\\' ;os.makedirs(dir_ndvi, exist_ok=True)
    
    #Cria uma lista como os arquivos de NDVI sem nuvens do respectivo dia
    ndvi_list = [name for name in os.listdir(dir_ndvi) if os.path.isfile(os.path.join(dir_ndvi, name)) and
                    re.match('^ndvi_2020.+_%s.npy$'%reg, name)]
    ndvi_list.sort()
    # print(ndvi_list)
    teste = np.load(f'{dir_ndvi}{ndvi_list[0]}', allow_pickle=True)
    NDVI_maxday = np.zeros(teste.shape)

    # Setar os valores todos para NAN pois max(nan,value) = value 
    NDVI_maxday[NDVI_maxday == 0]= np.nan
    del teste

    #Date to save de max day file
    date_file = ndvi_list[0][ndvi_list[0].find("ndvi_") + 5:ndvi_list[0].find("ndvi_") + 13]

    for f in ndvi_list:
        print(f)
        try:
            ndvi_file = np.load(f'{dir_ndvi}{f}', allow_pickle=True)

            NDVI_maxday = np.fmax(NDVI_maxday, ndvi_file)
        except Exception as erro:
            print(f'{erro} dia{dia}')
            continue

    NDVI_maxday.dump(f'{dir_ndvi}ndvi_{date_file}_{reg}_max.npy')

    

print('Processando dados diários da Região3')
reg = 'reg3'

# Processa NDVI de cada dia

for d in range(start_jul,end_jul+1):
    dia = f"{str(d).zfill(3)}"
    print(f'Dia {dia}')

    dir_ndvi = f'{dir_in}ndvi\\{ano}\\{mes}\\{dia}\\' ;os.makedirs(dir_ndvi, exist_ok=True)
    
    #Cria uma lista como os arquivos de NDVI sem nuvens do respectivo dia
    ndvi_list = [name for name in os.listdir(dir_ndvi) if os.path.isfile(os.path.join(dir_ndvi, name)) and
                    re.match('^ndvi_2020.+_%s.npy$'%reg, name)]
    ndvi_list.sort()
    # print(ndvi_list)
    teste = np.load(f'{dir_ndvi}{ndvi_list[0]}', allow_pickle=True)
    NDVI_maxday = np.zeros(teste.shape)

    # Setar os valores todos para NAN pois max(nan,value) = value 
    NDVI_maxday[NDVI_maxday == 0]= np.nan
    del teste

    #Date to save de max day file
    date_file = ndvi_list[0][ndvi_list[0].find("ndvi_") + 5:ndvi_list[0].find("ndvi_") + 13]

    for f in ndvi_list:
        print(f)
        try:
            ndvi_file = np.load(f'{dir_ndvi}{f}', allow_pickle=True)

            NDVI_maxday = np.fmax(NDVI_maxday, ndvi_file)
        except Exception as erro:
            print(f'{erro} dia{dia}')
            continue

    NDVI_maxday.dump(f'{dir_ndvi}ndvi_{date_file}_{reg}_max.npy')


print('Processando dados diários da Região4')
reg = 'reg4'

# Processa NDVI de cada dia

for d in range(start_jul,end_jul+1):
    dia = f"{str(d).zfill(3)}"
    print(f'Dia {dia}')

    dir_ndvi = f'{dir_in}ndvi\\{ano}\\{mes}\\{dia}\\' ;os.makedirs(dir_ndvi, exist_ok=True)
    
    #Cria uma lista como os arquivos de NDVI sem nuvens do respectivo dia
    ndvi_list = [name for name in os.listdir(dir_ndvi) if os.path.isfile(os.path.join(dir_ndvi, name)) and
                    re.match('^ndvi_2020.+_%s.npy$'%reg, name)]
    ndvi_list.sort()
    # print(ndvi_list)
    teste = np.load(f'{dir_ndvi}{ndvi_list[0]}', allow_pickle=True)
    NDVI_maxday = np.zeros(teste.shape)

    # Setar os valores todos para NAN pois max(nan,value) = value 
    NDVI_maxday[NDVI_maxday == 0]= np.nan
    del teste

    #Date to save de max day file
    date_file = ndvi_list[0][ndvi_list[0].find("ndvi_") + 5:ndvi_list[0].find("ndvi_") + 13]

    for f in ndvi_list:
        print(f)
        try:
            ndvi_file = np.load(f'{dir_ndvi}{f}', allow_pickle=True)

            NDVI_maxday = np.fmax(NDVI_maxday, ndvi_file)
        except Exception as erro:
            print(f'{erro} dia{dia}')
            continue

    NDVI_maxday.dump(f'{dir_ndvi}ndvi_{date_file}_{reg}_max.npy')

    