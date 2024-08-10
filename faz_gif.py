# -----------------------------------------------------------------------------------------------------------
import sys
sys.path.insert(1, 'G:\Drives compartilhados\GOES16_FDCF\Programacao\Scripts')
# -----------------------------------------------------------------------------------------------------------

# Required modules
from osgeo import gdal  # Python bindings for GDAL
import os  # Utilitario para trabalhar com chamadas de sistema
import re  # Utilitario para trabalhar com expressoes regulares
import imageio.v2 as imageio

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

extent_pantanal = [-59.3, -15.5, -54.7, -22.3] # Max Lat, Max lon, Min lat, Min Lon
extent_br = [-90.0, -40.0, -20.0, 10.0]


## Dias Julianos 
mes = 'Jul'
start_jul = 214
end_jul = 217

#Shapefille Pantanal

for d in range(start_jul,end_jul+1):
    print('Dia:',d)
    dir_dia = f'{dir_in}fdcf\\{mes}\\{d}\\'
    dir_img_dia = f'{dir_out}fdcf\\{mes}\\{d}\\'; os.makedirs(dir_img_dia,exist_ok=True)
    
    list_img = [name for name in os.listdir(dir_img_dia) if os.path.isfile(os.path.join(dir_img_dia, name))
                      and re.match("^fdcf_2020.+._pant.png$", name)]
    list_img.sort()

    imagens =[]
    for name in list_img:
        imagens.append(imageio.imread(f'{dir_img_dia}{name}'))

    print(f'Salvando gif dia {d}')
    dir_gif = f'{dir_out}fdcf\\{mes}\\gif\\'; os.makedirs(dir_gif,exist_ok=True)
    name_gif = f'gif_fdcf_2020{d}.gif'
    imageio.mimsave(f'{dir_gif}{name_gif}',imagens,duration = 0.3)


