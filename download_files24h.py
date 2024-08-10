# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')

# Required modules
import time as t
from osgeo import gdal  # Python bindings for GDAL
import os  # Utilitario para trabalhar com chamadas de sistema
from libs.utilities import download_CMI,download_PROD

print('Script Iniciado...')
start = t.time()
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

# # Desired extent
extent_br = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat - Brasil
extent_pantanal = [-59.3, -22.3, -54.7, -15.5] # Min lon, Min lat, Max lon, Max lat - Pantanal  

#Dia inicial e final em Juliano e o mês de interesse.
start_jul = 183
end_jul = 217
mes = 'Jul'
for d in range(start_jul,end_jul+1):

    # Criação da nomenclatura dos dias para criação dos diretórios de armazenamento
    dia = f"{str(d).zfill(3)}"
    dir_ch06 = f'{dir_in}ch06\\{mes}\\{dia}\\';os.makedirs(dir_ch06, exist_ok=True)
    # dir_ch07 = f'{dir_in}ch07\\{mes}\\{dia}\\';os.makedirs(dir_ch07, exist_ok=True)
    # dir_fdcf = f'{dir_in}fdcf\\{mes}\\{dia}\\';os.makedirs(dir_fdcf, exist_ok=True)

    # Hora ao qual se incia o acumulado diario neste caso 00:00-23:50 UTC
    hh= 00
    mn = 00
    for h in range(24):
        for m in range(6):
            yyyydddhhmn = f'2020{dia}{str(hh).zfill(2)}{str(mn).zfill(2)}'
            print(yyyydddhhmn)
            # Download files CH2
            # file_fdcf = download_PROD(yyyydddhhmn,'ABI-L2-FDCF',dir_ch06)
            # file_ch07 = download_CMI(yyyydddhhmn,7,dir_ch07)
            file_ch06 = download_CMI(yyyydddhhmn,6,dir_ch06)
            if file_ch06 == -1:
                print(f'No file ch06 for the date {yyyydddhhmn}\n')
                mn += 10
                continue
            mn += 10
        mn -= 60
        hh += 1

        

        