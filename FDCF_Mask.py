# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------
# Required modules
import time as t
from netCDF4 import Dataset  # Read / Write NetCDF4 files
import matplotlib.pyplot as plt  # Plotting library
from datetime import datetime  # Basic Dates and time types
import cartopy.crs as ccrs  # Plot maps
from osgeo import gdal  # Python bindings for GDAL
import numpy as np  # Scientific computing with Python
from libs.utilities import download_PROD
from libs.utilities import Degrees, reprojectBruno, save_txt  # Our function for reproject
import cartopy.io.shapereader as shpreader  # Import shapefiles
from shapely.geometry import Point
## Verificar esse trecho eu coloquei pq no tutorial eles disseram que era para fazer
# com que parassem os avisos da GDAL

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

extent_br = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat - Brasil
extent_pantanal = [-59.3, -22.3, -54.7, -15.5] # Min lon, Min lat, Max lon, Max lat - Pantanal



# Datetime to process
yyyymmddhhmn = '202007011300'      
product_name = "ABI-L2-FDCF"

#Download file fire_hotspot
file_fdcf = download_PROD(yyyymmddhhmn,product_name,dir_in+'fdcf\\')

#Shapefille Pantanal
biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
pantanal = biomas[5]
ax.add_geometries(pantanal, ccrs.PlateCarree(), edgecolor='orange', facecolor='none', linewidth=0.65, zorder=4)

# # selecionando os pixels
Fire_Mask = Dataset(f'{dir_in}fdcf\\{file_fdcf}')
fire_mask_values = Fire_Mask.variables['Mask'][:, :]

#Seleciona os pontos caracterizados pelo satelite, com maior probabilidade de serem realmente focos de incêncdios
selected_fires = (fire_mask_values == 10) |(fire_mask_values == 13) | (fire_mask_values == 30)| (fire_mask_values == 33)

Lat, Lon = Degrees(Fire_Mask)

# separando Latitudes e Longitudes dos pontos
p_lat = Lat[selected_fires]
p_lon = Lon[selected_fires]

matriz_pixels_fogo = []
for i in range(len(p_lat)):
        if pantanal.covers(Point(p_lon[i], p_lat[i])):
            p = (p_lat[i],p_lon[i])
            matriz_pixels_fogo.append(p)

#Reprojetando o arquivo 
file_reproject = remap(f'{dir_in}fdcf\\Jul\\207\\{file_fdcf}','Area',extent_pantanal,1)

data=file_reproject.ReadAsarray()

# #Salva o array em formato npy
data.dump(f'{dir_in}fdcf\\Jul\\207\\npy\\fdcf_20202070630_area.npy')


# Definindo tamanho da imagem de saida
d_p_i = 150
fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black',
                 facecolor='black')

ax = plt.axes(projection=ccrs.PlateCarree())
# Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
img_extent = [extent_pantanal[0], extent_pantanal[2], extent_pantanal[1], extent_pantanal[3]] # Min lon, Max lon, Min lat, Max lat

ax.set_extent(img_extent)

# # Plotandos os pontos de fogo
for i in range(len(matriz_pixels_fogo)):
    p_lon = matriz_pixels_fogo[i][1]
    p_lat = matriz_pixels_fogo[i][0]
    ax.plot(p_lon, p_lat, 'k+', ms=2.5, transform=ccrs.PlateCarree())

# Estados, coastline e shape Brasil
estados = shpreader.Reader(f'{dir_shapefiles}\\divisao_estados\\gadm40_BRA_1.shp').geometries()
ax.add_geometries(mt, ccrs.PlateCarree(), edgecolor='orange', facecolor='none', linewidth=0.65, zorder=4)
# ax.coastlines(resolution='10m', color='orange', linewidth=0.5)
# ax.add_feature(cartopy.feature.BORDERS, edgecolor='orange', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='black', alpha=0.7, linestyle='--', linewidth=0.2, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5))
gl.top_labels = False
gl.right_labels = False

#Plotando imagem
img = ax.imshow(data,origin='upper',cmap=cmap, extent=img_extent,vmin=vmin,vmax=vmax)

#Adicionando Colorbar
cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[0, 10, 13, 30, 33])
cb.ax.set_xticklabels(['0', '10', '13', '30', '33'])
cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de cores


# Adicionando descricao da imagem
description = f"GOES-16 Natural True Color,    Fire Hot Spot em {dd}-{datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%b')}-{yyyy} "
institution = f'CEPAGRI - UNICAMP'
cruz = '+'

# Criando novos eixos de acordo com a posicao da imagem
cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
cax1.patch.set_color('black')  # Alterando a cor do novo eixo
cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
cax1.text(0.193, 0.13, cruz, color='red', size=12)  # Adicionando simbolo "+"
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

#Salvando a imagem de saída
plt.savefig(f'{dir_out}fdcf\\FDCF_func_teste.png',bbox_inches='tight', pad_inches=0, dpi=d_p_i)
plt.show()

print(f'Tempo total {round((t.time() - start), 2)} segundos ')  
