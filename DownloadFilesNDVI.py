# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------
# Required modules
import time as t
from netCDF4 import Dataset  # Read / Write NetCDF4 files
import matplotlib.pyplot as plt  # Plotting library
import matplotlib
from datetime import timedelta, datetime  # Basic Dates and time types
import cartopy, cartopy.crs as ccrs  # Plot maps
import cartopy.io.shapereader as shpreader  # Import shapefiles
import cartopy.feature as cfeature # features
import numpy as np  # Scientific computing with Python
from osgeo import gdal  # Python bindings for GDAL
import os  # Utilitario para trabalhar com chamadas de sistema
import shutil
import re  # Utilitario para trabalhar com expressoes regulares
from utilities import download_CMI,download_PROD_simposio
from utilities import Degrees, reprojectBruno, save_txt  # Our function for reproject
from utilities import imprime_matriz

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

extent_br = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat - Brasil
extent_pantanal = [-59.3, -22.3, -54.7, -15.5] # Min lon, Min lat, Max lon, Max lat - Pantanal  

## Dias Julianos 
start_jul = 183
end_jul = 184
product_name = 'ABI-L2-ACMF'
mes = 'Jul'

# Range da quantidade de dias
for d in range(start_jul,end_jul+1):

    # Criação da nomenclatura dos dias para criação dos diretórios de armazenamento
    dia = f"{str(d).zfill(3)}"
    dir_ch2 = f'{dir_in}ch02\\{mes}\\{dia}\\';os.makedirs(dir_ch2, exist_ok=True)
    dir_ch3 = f'{dir_in}ch03\\{mes}\\{dia}\\';os.makedirs(dir_ch3, exist_ok=True)
    dir_clsm = f'{dir_in}clsm\\{mes}\\{dia}\\';os.makedirs(dir_clsm, exist_ok=True)
    
    # Hora ao qual se incia o acumulado diario neste caso 13-18 UTC
    hh= 13
    mn = 00
    for h in range(5):
        for m in range(6):
            yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
            print(yyyydddhhmn)
            #Download files CH2
            file_ch2 = download_CMI(yyyydddhhmn, 2, dir_ch2)
            reprojec_ch2 = ''
            if file_ch2 == -1:
                print(f'No file ch2 for the date {yyyydddhhmn}\n')
                mn += 10
                continue
            else:
                try:
                    print(f'reprojecting {file_ch2}')
                    reprojec_ch2 = reprojectBruno(f'{dir_ch2}{file_ch2}', 'CMI', extent_pantanal, 2, dir_ch2)
                    # os.remove(f'{dir_ch2}{file_ch2}')
                except Exception as erro:
                    os.remove(f'{dir_ch2}{file_ch2}')
                    yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
                    print(f'Arquivo ch02 {yyyydddhhmn} com erro {erro}')
                    mn += 10
                    continue

            # Download Files CH3        
            file_ch3 = download_CMI(yyyydddhhmn, 3, dir_ch3)
            reproject_ch3 = ''
            if file_ch3 == -1:
                os.remove(reprojec_ch2)
                print(f'No file ch3 for the date {yyyydddhhmn}\n')
                mn += 10
                continue
            else:
                try:
                    print(f'reprojecting {file_ch3}')
                    reproject_ch3 = reprojectBruno(f'{dir_ch3}{file_ch3}', 'CMI', extent_pantanal, 2, dir_ch3)
                    # os.remove(f'{dir_ch3}{file_ch3}')
                except Exception as erro:
                    os.remove(f'{dir_ch3}{file_ch3}')
                    os.remove(reprojec_ch2)
                    yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
                    print(f'Arquivo ch03 {yyyydddhhmn} com erro {erro}')
                    mn += 10
                    continue
                
            # Download Files Mask        
            file_mask = download_PROD_simposio(yyyydddhhmn, product_name, dir_clsm)
            reproject_file_mask = ''
            if file_mask == -1:
                os.remove(reprojec_ch2)
                os.remove(reproject_ch3)
                print(f'No file mask for the date {yyyydddhhmn}\n')
                mn += 10
                continue
            else:
                try:
                    print(f'reprojecting {file_mask}')
                    reproject_file_mask = reprojectBruno(f'{dir_clsm}{file_mask}', 'BCM', extent_pantanal, 2, dir_clsm)
                    # os.remove(f'{dir_clsm}{file_mask}')
                except Exception as erro:
                    os.remove(reprojec_ch2)
                    os.remove(reproject_ch3)
                    os.remove(f'{dir_clsm}{file_mask}')
                    yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
                    print(f'Arquivo clsm {yyyydddhhmn} com erro {erro}')
                    mn += 10
                    continue
            mn += 10
        mn -= 60
        hh += 1

        # Donwload file 18:00
        yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
        print(yyyydddhhmn)
        file_ch2 = download_CMI(yyyydddhhmn, 2, dir_ch2)
        reprojec_ch2 = ''
        if file_ch2 == -1:
            print(f'No file ch2 for the date {yyyydddhhmn}\n')
            continue
        else:
            try:
                print(f'reprojecting {file_ch2}')
                reprojec_ch2 = reprojectBruno(f'{dir_ch2}{file_ch2}', 'CMI', extent_pantanal, 2, dir_ch2)
                # os.remove(f'{dir_ch2}{file_ch2}')
            except Exception as erro:
                os.remove(f'{dir_ch2}{file_ch2}')
                yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
                print(f'Arquivo ch02 {yyyydddhhmn} com erro {erro}')
                continue

        file_ch3 = download_CMI(yyyydddhhmn, 3, dir_ch3)
        reproject_ch3 = ''
        if file_ch3 == -1:
            os.remove(reprojec_ch2)
            f'No file ch3 for the date {yyyydddhhmn}\n'
            continue
        else:
            try:
                print(f'reprojecting {file_ch3}')
                reproject_ch3 = reprojectBruno(f'{dir_ch3}{file_ch3}', 'CMI', extent_pantanal, 2, dir_ch3)
                # os.remove(f'{dir_ch3}{file_ch3}')
            except Exception as erro:
                os.remove(f'{dir_ch3}{file_ch3}')
                os.remove(reprojec_ch2)
                yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
                print(f'Arquivo ch03 {yyyydddhhmn} com erro {erro}')
                continue

        file_mask = download_PROD_simposio(yyyydddhhmn, product_name, dir_clsm)
        reproject_file_mask = ''
        if file_mask == -1:
            os.remove(reprojec_ch2)
            os.remove(reproject_ch3)
            f'No file mask for the date {yyyydddhhmn}\n'
            continue
        else:
            try:
                print(f'reprojecting {file_mask}')
                reproject_file_mask = reprojectBruno(f'{dir_clsm}{file_mask}', 'BCM', extent_pantanal, 2, dir_clsm)
                # os.remove(f'{dir_clsm}{file_mask}')
            except Exception as erro:
                os.remove(reprojec_ch2)
                os.remove(reproject_ch3)
                os.remove(f'{dir_clsm}{file_mask}')
                yyyydddhhmn = f'2020{dia}{hh}{str(mn).zfill(2)}'
                print(f'Arquivo clsm {yyyydddhhmn} com erro {erro}')
                continue


# #Este Trecho de Código desliga o computador automaticamente ao final do script
# print('Desligando em:')
# for i in range(10):
#     print(i)
# os.system("shutdown /s /t 1")
