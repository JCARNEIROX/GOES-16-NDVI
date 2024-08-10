# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------
# Required modules
import time as t
import matplotlib.pyplot as plt  # Plotting library
import cartopy.crs as ccrs  # Plot maps
import os  # Miscellaneous operating system interfaces
from osgeo import gdal  # Python bindings for GDAL
import numpy as np  # Scientific computing with Python
from libs.utilities import download_PROD
import cartopy.io.shapereader as shpreader  # Import shapefiles
from libs.remap import remap
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
yyyymmddhhmn = '202007011500'     
product_name = "ABI-L2-FDCF"
dia = 183


##------------------------------------
variable = 'Area' 
vmin = 0.004
vmax = 4
cmap = 'Reds'
# #Download file fire_hotspot
file_fdcf = download_PROD(yyyymmddhhmn,product_name,f'{dir_in}fdcf\\Jul\\207\\')

# # Reprjetando o arquivo da Àrea para recortar o pantanal
area_array_reproject = remap(f'{dir_in}fdcf\\Jul\\207\\{file_fdcf}','Area',extent_pantanal,1)

# #Salva o array em formato npy
area_array_reproject.dump(f'{dir_in}fdcf\\Jul\\207\\npy\\fdcf_20202070630_area.npy')


print('Plotando as imagens')
d_p_i = 150
fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black',
                facecolor='black')
ax = plt.axes(projection=ccrs.PlateCarree())


img_extent = [extent_pantanal[0], extent_pantanal[2], extent_pantanal[1], extent_pantanal[3]]  # Min Lat, Max Lat, Min Lon, Max Lon

#Adicionando shapefiles
#Shapefille Pantanal
biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
pantanal = biomas[5]
ax.add_geometries(pantanal, ccrs.PlateCarree(), edgecolor='orange', facecolor='none', linewidth=0.65, zorder=4)

xstep = (abs(extent_pantanal[0]) - abs(extent_pantanal[2]))/756
ystep = (abs(extent_pantanal[1]) - abs(extent_pantanal[3]))/512

gl = ax.gridlines(color='gray', linestyle='--', linewidth=0.25, alpha=0.4, xlocs=np.arange(extent_pantanal[0], extent_pantanal[2], xstep),
            ylocs=np.arange(extent_pantanal[1], extent_pantanal[3], ystep), zorder=5)

# Abri o arquivo npy
area_array_reproject = np.load(f'{dir_in}fdcf\\Jul\\207\\npy\\fdcf_20202070630_area.npy', allow_pickle=True)
print(area_array_reproject.shape)

#Plotando imagem
img = ax.imshow(area_array_reproject,origin='upper',cmap=cmap, extent=img_extent,vmin=vmin,vmax=vmax)

#Adicionando Colorbar
cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
cb = plt.colorbar(img, orientation="horizontal", cax=cax0,ticks=[0, 0.5, 1, 2, 3, 3.5])
cb.ax.set_xticklabels(['0', '0.5', '1', '2', '3', '3.5'])
cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de cores


#Salvando a imagem de saída
plt.savefig(f'{dir_out}fdcf\\Jul\\183\\fdcf_20201831500_pant_area1.png',bbox_inches='tight', pad_inches=0, dpi=d_p_i)
plt.show()



