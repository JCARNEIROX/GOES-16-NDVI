# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------
# Required modules
import matplotlib.pyplot as plt  # Plotting library
import cartopy.crs as ccrs  # Plot maps
import cartopy.io.shapereader as shpreader  # Import shapefiles
import numpy as np  # Scientific computing with Python
from osgeo import gdal  # Python bindings for GDAL
import os  # Utilitario para trabalhar com chamadas de sistema
import re  # Utilitario para trabalhar com expressoes regulares
from libs.utilities import save_txt  # Our function for reproject
import shapely

print('Script Iniciado...')
gdal.PushErrorHandler('CPLQuietErrorHandler')  # Ignore GDAL warnings

#======================================================================================================
# Input and output directories
dir_in = "G:\\Drives compartilhados\\GOES16_FDCF\\Programacao\\dir_in\\"
dir_main = "G:\\Drives compartilhados\GOES16_FDCF\\Programacao\\"
dir_out = dir_main + "output\\"
dir_libs = dir_main + "libs\\"
dir_shapefiles = dir_main + "shapefiles\\"
dir_colortables = dir_main + "colortables\\"
dir_logos = dir_main + "logos\\"

#Desired Extent
extent_br = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat - Brasil
extent_pantanal = [-59.3, -22.3, -54.7, -15.5] # Min lon, Min lat, Max lon, Max lat - Pantanal


## Dias Julianos 
# ano = 2020 
# start_jul = 183  
# end_jul = 213
# mes = 'Jul'

# #Fazendo o a matriz diaria de pontos 

# for d in range(start_jul,end_jul+1):
    
#     dia = str(d).zfill(3)
#     print('Dia:',dia)
#     dir_dia = f'{dir_in}fdcf\\{ano}\\{mes}\\{dia}\\'
#     dir_txt_dia = f'{dir_dia}txt\\'; os.makedirs(dir_txt_dia,exist_ok=True)
#     dir_out_dia = f'{dir_out}fdcf\\{mes}\\{dia}\\'; os.makedirs(dir_out_dia,exist_ok=True)
    
#     list_txt = [name for name in os.listdir(dir_txt_dia) if os.path.isfile(os.path.join(dir_txt_dia, name))
#                       and re.match("^fdcf_2020.+.txt$", name)]
#     list_txt.sort()

#     matriz_diaria = []
#     for file in list_txt:
#         array = np.loadtxt(f'{dir_txt_dia}{file}', 'float', delimiter=',')
#         for i in range(len(array)):
#             try:
#                 p_lon = array[i][1]
#                 p_lat = array[i][0]
#                 p = (p_lat, p_lon)
#                 # matriz_diaria.append(p)
#                 if not p in matriz_diaria:
#                     matriz_diaria.append(p)
#             except:
#                 p_lon = array[1]
#                 p_lat = array[0]
#                 p = (p_lat, p_lon)
#                 # matriz_diaria.append(p)
#                 if not p in matriz_diaria:
#                     matriz_diaria.append(p)
                
#     # # Salvando a matriz dos pontos em txt
#     print(f'Salvando Matriz diaria {d}')
#     name_save = f"fdcf_{dia}_pant"
#     save_txt(matriz_diaria, name_save, dir_dia+'txt','fdcf')


#======================================================================================================        
# #Fazendo o plot dos focos a cada 10 min

#Shapefille Pantanal
# estados = shpreader.Reader(f'{dir_shapefiles}\\divisao_estados\\gadm40_BRA_1.shp').geometries()
# biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
# pantanal = biomas[5]


# for d in range(start_jul,end_jul+1):
#     print('Dia:',d)
#     dir_dia = f'{dir_in}fdcf\\Jul\\{d}\\'; os.makedirs(dir_dia,exist_ok=True)
#     dir_txt_dia = f'{dir_dia}txt\\'; os.makedirs(dir_txt_dia,exist_ok=True)
#     dir_out_dia = f'{dir_out}fdcf\\Jul\\{d}\\'; os.makedirs(dir_out_dia,exist_ok=True)
    
#     list_txt = [name for name in os.listdir(dir_txt_dia) if os.path.isfile(os.path.join(dir_txt_dia, name))
#                       and re.match("^fdcf_2020.+.txt$", name)]
#     list_txt.sort()

    
#     for file in list_txt:
#         # Definindo tamanho da imagem de saida
#         d_p_i = 150
#         fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black',
#                         facecolor='black')
#         img_extent = [extent_pantanal[0], extent_pantanal[2], extent_pantanal[1], extent_pantanal[3]]
#         ax = plt.axes(projection=ccrs.PlateCarree())
#         ax.set_extent(img_extent)

#         #Adicionando estado
#         ax.add_geometries(pantanal, ccrs.PlateCarree(), edgecolor='orange', facecolor='none', linewidth=0.65, zorder=4)
#         gl = ax.gridlines(color='gray', linestyle='--', linewidth=0.25, alpha=0.4, xlocs=np.arange(-180, 180, 5),
#                     ylocs=np.arange(-90, 90, 5), zorder=5)
#         yyyymmddhhmn = file[5:16]
#         array = np.loadtxt(f'{dir_txt_dia}{file}', 'float', delimiter=',')
#         for i in range(len(array)):
#             try:
#                 p_lon = array[i][1]
#                 p_lat = array[i][0]
#                 p = (p_lat, p_lon)
#                 ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
#             except:
#                 p_lon = array[1]
#                 p_lat = array[0]
#                 p = (p_lat, p_lon)
#                 ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())

#         print('Salvando a imagem', file)
#         plt.savefig(f'{dir_out_dia}fdcf_{yyyymmddhhmn}_pant',bbox_inches='tight', pad_inches=0, dpi=d_p_i)
#         # plt.show()
#         plt.close()

#======================================================================================================
# #Fazendo o plot dos focos de um unico dia em uma unica imagem
#======================================================================================================

# # Shapefille Pantanal
# estados = shpreader.Reader(f'{dir_shapefiles}\\divisao_estados\\gadm40_BRA_1.shp').geometries()
# biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
# pantanal = biomas[5]


# for d in range(start_jul,end_jul+1):

#     print('Dia:',d)
#     dir_dia = f'{dir_in}fdcf\\Jul\\{d}\\'
#     dir_txt_dia = f'{dir_dia}txt\\'; os.makedirs(dir_txt_dia,exist_ok=True)
#     dir_out_dia = f'{dir_out}fdcf\\Jul\\{d}\\'; os.makedirs(dir_out_dia,exist_ok=True)
    
#     list_txt = [name for name in os.listdir(dir_txt_dia) if os.path.isfile(os.path.join(dir_txt_dia, name))
#                       and re.match("^fdcf_2020.+.txt$", name)]
#     list_txt.sort()

#     # Definindo tamanho da imagem de saida
#     d_p_i = 150
#     fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black',
#                     facecolor='black')
#     img_extent = [extent_pantanal[0], extent_pantanal[2], extent_pantanal[1], extent_pantanal[3]]
#     ax = plt.axes(projection=ccrs.PlateCarree())
#     ax.set_extent(img_extent)

#     #Adicionando estado
#     ax.add_geometries(pantanal, ccrs.PlateCarree(), edgecolor='orange', facecolor='none', linewidth=0.65, zorder=4)
#     gl = ax.gridlines(color='gray', linestyle='--', linewidth=0.25, alpha=0.4, xlocs=np.arange(-180, 180, 5),
#                 ylocs=np.arange(-90, 90, 5), zorder=5)

#     #Titulo da imagem
#     description = f'Focos de incêndio do dia {d}'
#     institution = "CEPAGRI - UNICAMP"
#     cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
#     cax1.patch.set_color('black')  # Alterando a cor do novo eixo
#     cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
#     cax1.text(0.78, 0.13, institution, color='yellow', size=10)  # Adicionando texto
#     cax1.xaxis.set_visible(False)  # Removendo rótulos do eixo X
#     cax1.yaxis.set_visible(False)  # Removendo rótulos do eixo Y

#     print('Plotando os pontos')
    
#     for file in list_txt:
#         yyyymmddhhmn = file[5:16]
#         array = np.loadtxt(f'{dir_txt_dia}{file}', 'float', delimiter=',')
#         for i in range(len(array)):
#             try:
#                 p_lon = array[i][1]
#                 p_lat = array[i][0]
#                 p = (p_lat, p_lon)
#                 ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                
#             except:
#                 p_lon = array[1]
#                 p_lat = array[0]
#                 p = (p_lat, p_lon)
#                 ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                
#     print(f'Salvando a imagem dia')
#     plt.savefig(f'{dir_out}fdcf\\fdcf_{d}_pant',bbox_inches='tight', pad_inches=0, dpi=d_p_i)
#     # plt.show()
#     plt.close()

#======================================================================================================
# #Fazendo o plot dos focos de incêndio de uma semana em uma unica imagem
#======================================================================================================

# #Shapefille Pantanal
# estados = shpreader.Reader(f'{dir_shapefiles}\\divisao_estados\\gadm40_BRA_1.shp').geometries()
# biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
# pantanal = biomas[5]   

# # # Definindo tamanho da imagem de saida
# d_p_i = 150
# fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black',
#                 facecolor='black')
# img_extent = [extent_pantanal[0], extent_pantanal[2], extent_pantanal[1], extent_pantanal[3]]
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_extent(img_extent)

# #Adicionando estado
# ax.add_geometries(pantanal, ccrs.PlateCarree(), edgecolor='orange', facecolor='none', linewidth=0.65, zorder=4)
# gl = ax.gridlines(color='gray', linestyle='--', linewidth=0.25, alpha=0.4, xlocs=np.arange(-180, 180, 5),
#             ylocs=np.arange(-90, 90, 5), zorder=5)

# #Titulo da imagem
# description = f'Focos de incêndio da semana {start_jul}-{end_jul}'
# institution = "CEPAGRI - UNICAMP"
# cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
# cax1.patch.set_color('black')  # Alterando a cor do novo eixo
# cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
# cax1.text(0.78, 0.13, institution, color='yellow', size=10)  # Adicionando texto
# cax1.xaxis.set_visible(False)  # Removendo rótulos do eixo X
# cax1.yaxis.set_visible(False)  # Removendo rótulos do eixo Y

# for d in range(start_jul,end_jul+1):
#     print('Dia:',d)
#     dir_dia = f'{dir_in}fdcf\\Jul\\{d}\\'
#     dir_txt_dia = f'{dir_dia}txt\\'; os.makedirs(dir_txt_dia,exist_ok=True)
#     dir_out_dia = f'{dir_out}fdcf\\Jul\\{d}\\'; os.makedirs(dir_out_dia,exist_ok=True)
    
#     list_txt = [name for name in os.listdir(dir_txt_dia) if os.path.isfile(os.path.join(dir_txt_dia, name))
#                       and re.match("^fdcf_2020.+.txt$", name)]
#     list_txt.sort()

#     print('Plotando os pontos')
    
#     for file in list_txt:
#         yyyymmddhhmn = file[5:16]
#         array = np.loadtxt(f'{dir_txt_dia}{file}', 'float', delimiter=',')
#         for i in range(len(array)):
#             try:
#                 p_lon = array[i][1]
#                 p_lat = array[i][0]
#                 p = (p_lat, p_lon)
#                 ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                
#             except:
#                 p_lon = array[1]
#                 p_lat = array[0]
#                 p = (p_lat, p_lon)
#                 ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                
# print(f'Salvando a imagem semanal')
# plt.savefig(f'{dir_out}fdcf\\fdcf_{start_jul}-{end_jul}_pant',bbox_inches='tight', pad_inches=0, dpi=d_p_i)
# # plt.show()
# plt.close()


#======================================================================================================
#Fazendo o plot dos focos de incêndio de um mês em uma única imagem
#======================================================================================================

# Shapefille Pantanal
estados = shpreader.Reader(f'{dir_shapefiles}\\divisao_estados\\gadm40_BRA_1.shp').geometries()
biomas = list(shpreader.Reader(f'{dir_shapefiles}\\Biomas_250mil\\lm_bioma_250.shp').geometries())
pantanal = biomas[5]   
    
# Definindo tamanho da imagem de saida
d_p_i = 150
fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black',
                facecolor='black')
img_extent = [extent_pantanal[0], extent_pantanal[2], extent_pantanal[1], extent_pantanal[3]]
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(img_extent)

#Adicionando estado
ax.add_geometries(pantanal, ccrs.PlateCarree(), edgecolor='orange', facecolor='none', linewidth=0.65, zorder=4)
# gl = ax.gridlines(color='gray', linestyle='--', linewidth=0.25, alpha=0.4, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), zorder=5)

# # #Adicionando areas de interesse Jul 2020
# regiao1 = [-58.4, -17.42, -58.06, -17] # Min lon, Min lat, Max lon, Max lat
# area1 = shapely.geometry.box(-58.4,-17.42,-58.06,-17)

# regiao2 = [-57.38, -17.25, -57.03, -16.71]
# area2 = shapely.geometry.box(-57.38,-17.25,-57.03, -16.71)

# regiao3 = [-57.71, -19.22, -57.27, -18.55]
# area3 = shapely.geometry.box(-57.71,-19.22,-57.27,-18.55)

# regiao4 = [-58.14, -19.93, -57.85, -19.41]
# area4 = shapely.geometry.box(-58.14,-19.93,-57.85,-19.41)


# ax.add_geometries(area1, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)
# ax.add_geometries(area2, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)
# ax.add_geometries(area3, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)
# ax.add_geometries(area4, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)

#Adicionando areas de interesse Jun 2024
# regiao1 = [-57.76, -19.78, -57.01, -28.37] # Min lon, Min lat, Max lon, Max lat
# area1 = shapely.geometry.box(-57.76,-19.78,-57.01,-28.37)

# regiao2 = [-57.65, -19.83, -56.99, -18.7]
# area2 = shapely.geometry.box(-57.65,-19.83,-56.99, -18.7)

# regiao3 = [-56.51, -19.8, -56.1, -19.47]
# area3 = shapely.geometry.box(-56.51,-19.8,-56.1,-19.47)

# ax.add_geometries(area1, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)
# ax.add_geometries(area2, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)
# ax.add_geometries(area3, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)

# #Adicionando areas de interesse Jul 2024
# regiao1 = [-58.47, -17.01, -58.33, -16.69] # Min lon, Min lat, Max lon, Max lat
# area1 = shapely.geometry.box(-58.47,-17.01,-58.33,-16.69)

# regiao2 = [-57.65, -19.83, -56.99, -18.7]
# area2 = shapely.geometry.box(-57.65,-19.83,-56.99, -18.7)

# regiao3 = [-56.51, -19.8, -56.1, -19.47]
# area3 = shapely.geometry.box(-56.51,-19.8,-56.1,-19.47)

# ax.add_geometries(area1, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)
# ax.add_geometries(area2, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)
# ax.add_geometries(area3, ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=0.65, zorder=4)


ano = 2024 
start_jun = 153
end_jun = 182
mes = 'Jun'

#Titulo da imagem
description = f'Focos de incêndio do mês de Junho {ano}'
institution = "CEPAGRI - UNICAMP"
cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
cax1.patch.set_color('black')  # Alterando a cor do novo eixo
cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
cax1.text(0.78, 0.13, institution, color='yellow', size=10)  # Adicionando texto
cax1.xaxis.set_visible(False)  # Removendo rótulos do eixo X
cax1.yaxis.set_visible(False)  # Removendo rótulos do eixo Y

cont1 =0
cont2 = 0
cont3 =0
cont4 =0
matriz_mensal = []
for d in range(start_jun,end_jun+1):
    dia = str(d).zfill(3)
    print('Dia:',dia)
    dir_dia = f'{dir_in}fdcf\\{ano}\\{mes}\\{dia}\\'
    dir_txt_dia = f'{dir_dia}txt\\'; os.makedirs(dir_txt_dia,exist_ok=True)
    dir_out_dia = f'{dir_out}fdcf\\{mes}\\{dia}\\'; os.makedirs(dir_out_dia,exist_ok=True)
    
    list_txt = [name for name in os.listdir(dir_txt_dia) if os.path.isfile(os.path.join(dir_txt_dia, name))
                      and re.match("^fdcf_%s.+.txt$"%ano, name)]
    list_txt.sort()

    print('Plotando os pontos')
    
    for file in list_txt:
        yyyymmddhhmn = file[5:16]
        array = np.loadtxt(f'{dir_txt_dia}{file}', 'float', delimiter=',')
        for i in range(len(array)):
            try:
                p_lon = array[i][1]
                p_lat = array[i][0]
                p = (p_lat, p_lon)
                
                # ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                # matriz_mensal.append(p)
                
                if not p in matriz_mensal:
                    ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                    # matriz_mensal.append(p)
                    # if area1.covers(Point(p_lon, p_lat)):
                    #     cont1+=1
                    # elif area2.covers(Point(p_lon, p_lat)):
                    #     cont2+=1
                    # elif area3.covers(Point(p_lon, p_lat)):
                    #     cont3+=1
                    # elif area4.covers(Point(p_lon, p_lat)):
                    #     cont4+=1
                
            except:
                p_lon = array[1]
                p_lat = array[0]
                p = (p_lat, p_lon)

                # ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                # matriz_mensal.append(p)
                
                if not p in matriz_mensal:
                    ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                    # matriz_mensal.append(p)
                    # if area1.covers(Point(p_lon, p_lat)):
                    #     cont1+=1
                    # elif area2.covers(Point(p_lon, p_lat)):
                    #     cont2+=1
                    # elif area3.covers(Point(p_lon, p_lat)):
                    #     cont3+=1
                    # elif area4.covers(Point(p_lon, p_lat)):
                    #     cont4+=1

ano = 2024 
start_jul = 183
end_jul = 213
mes = 'Jul'


for d in range(start_jul,end_jul+1):
    dia = str(d).zfill(3)
    print('Dia:',dia)
    dir_dia = f'{dir_in}fdcf\\{ano}\\{mes}\\{dia}\\'
    dir_txt_dia = f'{dir_dia}txt\\'; os.makedirs(dir_txt_dia,exist_ok=True)
    dir_out_dia = f'{dir_out}fdcf\\{mes}\\{dia}\\'; os.makedirs(dir_out_dia,exist_ok=True)
    
    list_txt = [name for name in os.listdir(dir_txt_dia) if os.path.isfile(os.path.join(dir_txt_dia, name))
                      and re.match("^fdcf_%s.+.txt$"%ano, name)]
    list_txt.sort()

    print('Plotando os pontos')
    
    for file in list_txt:
        yyyymmddhhmn = file[5:16]
        array = np.loadtxt(f'{dir_txt_dia}{file}', 'float', delimiter=',')
        for i in range(len(array)):
            try:
                p_lon = array[i][1]
                p_lat = array[i][0]
                p = (p_lat, p_lon)
                
                # ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                # matriz_mensal.append(p)
                
                if not p in matriz_mensal:
                    ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                    # matriz_mensal.append(p)
                    # if area1.covers(Point(p_lon, p_lat)):
                    #     cont1+=1
                    # elif area2.covers(Point(p_lon, p_lat)):
                    #     cont2+=1
                    # elif area3.covers(Point(p_lon, p_lat)):
                    #     cont3+=1
                    # elif area4.covers(Point(p_lon, p_lat)):
                    #     cont4+=1
                
            except:
                p_lon = array[1]
                p_lat = array[0]
                p = (p_lat, p_lon)

                # ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                # matriz_mensal.append(p)
                
                if not p in matriz_mensal:
                    ax.plot(p_lon, p_lat, 'r+', ms=1.5, transform=ccrs.PlateCarree())
                    # matriz_mensal.append(p)
                    # if area1.covers(Point(p_lon, p_lat)):
                    #     cont1+=1
                    # elif area2.covers(Point(p_lon, p_lat)):
                    #     cont2+=1
                    # elif area3.covers(Point(p_lon, p_lat)):
                    #     cont3+=1
                    # elif area4.covers(Point(p_lon, p_lat)):
                    #     cont4+=1



# dir_month = f"{dir_in}fdcf\\{ano}\\{mes}\\"
# name_save = f"fdcf_{mes}_pant"
# save_txt(matriz_mensal, name_save, dir_month,'fdcf')      
# print(f'Total de focos',len(matriz_mensal))

# print(f'Total de focos reg1 ',cont1)
# print(f'Total de focos reg2 ',cont2)
# print(f'Total de focos reg3 ',cont3)
# print(f'Total de focos reg4 ',cont4)

# print(f'Salvando a imagem dia')
# dir_artigo = r'G:\Drives compartilhados\\GOES16_FDCF\\Programacao\\output\\ndvi\\Artigo\\'
# plt.savefig(f'{dir_artigo}FDCF_{ano}{mes}_pant',bbox_inches='tight', pad_inches=0, dpi=d_p_i)


plt.show()

