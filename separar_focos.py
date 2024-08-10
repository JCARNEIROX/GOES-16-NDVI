# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------
# Required modules
import time as t
from netCDF4 import Dataset  # Read / Write NetCDF4 files
import cartopy.io.shapereader as shpreader  # Import shapefiles
from shapely.geometry import Point
from osgeo import gdal  # Python bindings for GDAL
import os  # Utilitario para trabalhar com chamadas de sistema
import re  # Utilitario para trabalhar com expressoes regulares
from utilities import Degrees, save_txt  # Our function for reproject


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
ano = 2024
start_jul = 153
end_jul = 182
mes = 'Jun'

#Shapefille Pantanal
biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
pantanal = biomas[5]

for d in range(start_jul,end_jul+1):

    dia = str(d).zfill(3)
    print('Dia:',dia)

    dir_dia = f'{dir_in}fdcf\\{ano}\\{mes}\\{dia}\\'; os.makedirs(dir_dia, exist_ok=True)
    list_fdcf = [name for name in os.listdir(dir_dia) if os.path.isfile(os.path.join(dir_dia, name))
                      and re.match("^OR_ABI-L2-FDCF-M6_G16_s.+.nc$", name)]
    list_fdcf.sort()

    for f in list_fdcf:
        name_file = f[23:34]
        print(name_file)
        # # selecionando os pixels
        Fire_Mask = Dataset(f'{dir_dia}{f}')
        fire_mask_values = Fire_Mask.variables['Mask'][:, :]

        selected_fires = (fire_mask_values == 10) |(fire_mask_values == 13) | (fire_mask_values == 30)| (fire_mask_values == 33)
        # print(selected_fires.shape)
        # print(type(selected_fires))

        Lat, Lon = Degrees(Fire_Mask)

        # print(Lon.shape,Lat.shape)

        # separando Latitudes e Longitudes dos pontos
        p_lat = Lat[selected_fires]
        p_lon = Lon[selected_fires]

        matriz_pixels_fogo = []
        start = t.time()
        for i in range(len(p_lat)):
                if pantanal.covers(Point(p_lon[i], p_lat[i])):
                    p = (p_lat[i],p_lon[i])
                    matriz_pixels_fogo.append(p)
                    # if not p in matriz_pixels_fogo:
                        # matriz_diaria.append(p)

        print(f'Quantidade de focos: {len(matriz_pixels_fogo)}')
        # print(f'Tempo demorado {round(t.time() - start,2)}')

        # # Salvando a matriz dos pontos em txt
        print(f'Salvando Matriz txt file {name_file}')
        dir_txt = dir_dia+'txt\\'; os.makedirs(dir_txt, exist_ok=True)
        name_save = f"fdcf_{name_file}_pant"
        save_txt(matriz_pixels_fogo, name_save, dir_txt,'fdcf')