import math
import numpy as np
#import xraydb
import csv
import pandas as pd
import gzip
import shutil
import os
import sys

#np.set_printoptions(threshold=np.inf)
pd.set_option('display.max_rows',23500)



for phantom_num in range(52):
    ###################
    #parameter
    #phantom(zsize is differ from vcarious phantom!)

    csv_input4 = pd.read_csv(filepath_or_buffer='../phantom_information_isocenter.csv',header=None,encoding='utf-8-sig')#gender,age,height,weight,slice_size
    #csv_input4 = pd.read_csv(filepath_or_buffer='../phantom_information_isocenter_18pati_forpaper.csv',header=None,encoding='utf-8-sig')#gender,age,height,weight,slice_size
    gender=csv_input4.loc[phantom_num,0]
    age=csv_input4.values[phantom_num,1]
    height=csv_input4.values[phantom_num,2]
    weight=csv_input4.values[phantom_num,3]
    #zsize=csv_input4.values[phantom_num,4]
    reconstruction_slices=csv_input4.values[phantom_num,4]
    #sourcez=csv_input4.values[phantom_num,5]
    reconstruction_start_z=csv_input4.values[phantom_num,5]
    spectre=csv_input4.loc[phantom_num,7]
    #reconstruction_start_z=-15
    #reconstruction_slices=140
    recon_factor=1000
    #x-ray("energy" is for loop)
    #sourcez=25.0#10cm couch
    detector_size=400

    ######same parameter in this time.(2022/2/2)
    #phantom
    xsize=512
    ysize=512
    scale_x=0.1
    scale_y=0.1
    scale_z=0.2
    num_material=6
    #iso_inverse=zsize-sourcez/scale_z
    
    #x-ray
    projection_number=360
    #sourcex=25.6
    #sourcey=-74.4

    #GPU=0##In Dl-box(2 GPU machine is implemented), GPU is 0 or 1
    
    ################

    i = 0
    ww = []
    density = [] 

    """
    # on couch
    new_dir_path='%s_%dy/'%(gender,age)
    if not os.path.exists(new_dir_path):
        os.makedirs(new_dir_path)
    """
    
    new_dir_path='no_couch/%s_%dy/'%(gender,age)
    #new_dir_path='no_couch/18pati_forpaper/%s_%dy/'%(gender,age)
    if not os.path.exists(new_dir_path):
        os.makedirs(new_dir_path)


    #f = open(new_dir_path+'FDK_%dpro_%s%dy%dh%dw'%(projection_number,gender,age,height,weight)+'.txt','w')
    f = open(new_dir_path+'FDK_%dpro_%s%dy%dh%dw'%(projection_number,gender,age,height,weight)+'.txt','w')
    

    f.write('//////////////////////////////////////////////////////// \n')
    f.write('//\n')
    f.write('//    This is info file for CBCT reconstruction \n')
    f.write('//    Do NOT remove this header, blank and comment lines \n')
    f.write('//\n')
    f.write('//////////////////////////////////////////////////////// \n')
    f.write('\n')
    f.write('// X-ray Energy Spectrum\n')
    f.write('kV\n')
    f.write('\n')
    f.write('// total projections scanned\n')
    f.write('%d\n'%(projection_number))
    f.write('\n')
    f.write('// projection width \n')
    f.write('%d\n'%(detector_size))
    f.write('\n')
    f.write('// projection height \n')
    f.write('%d\n'%(detector_size))
    f.write('\n')
    f.write('// projection image file header bytes  \n')
    f.write('0\n')
    f.write('\n')
    f.write('// projection image data directory for input \n')
    #f.write('../../../../../../nfs_L78/Raytracing_output_ICRPdatabase_todaispectre/360projection/no_couch/%s_%dy/%s%dy%dh%dw/%dpro_raytracing_%s%dy%dh%dw\n' %(gender,age,gender,age,height,weight,projection_number,gender,age,height,weight))
    f.write('../../../../../../nfs_L78/Raytracing_output_ICRPdatabase_todaispectre/360projection/no_couch/%s_%dy/%s%dy%dh%dw/%dpro_%s%dy%dh%dw_raytracing+medianpredictedscatterepoch900+noiselog0.065\n' %(gender,age,gender,age,height,weight,projection_number,gender,age,height,weight))
    #f.write('../../../../../../nfs_L78/Raytracing_output_ICRPdatabase_todaispectre/360projection/no_couch/%s_%dy/%s%dy%dh%dw/log_noise0.044per_minx0.1_%dpro_raytracing_%s%dy%dh%dw\n' %(gender,age,gender,age,height,weight,projection_number,gender,age,height,weight))
    f.write('\n')
    f.write('// projection image file number digits (0-filled) \n')
    f.write('8\n')
    f.write('\n')
    f.write('// first projection image file number in HEX \n')
    f.write('hore\n')
    f.write('\n')
    f.write('// projection image file extension \n')
    f.write('raw\n')
    f.write('\n')
    f.write('// angle data file \n')
    f.write('%d_littele_slide.flexmap\n'%(projection_number))
    f.write('\n')
    f.write('// reconstruction size (square format)  \n')
    f.write('%d\n'%(xsize))
    f.write('\n')
    f.write('// reconstruction target class string (used for gating, etc.) \n')
    f.write('*\n')
    f.write('\n')
    f.write('// reconstruction start z [cm]\n')
    f.write('%d\n'%(reconstruction_start_z ))
    f.write('\n')
    f.write('// reconstruction z step[cm]\n')
    f.write('%f\n'%(scale_z))
    f.write('\n')
    f.write('/// reconstruction slices \n')
    f.write('%d\n'%(reconstruction_slices))
    f.write('\n')
    f.write('// reconstruction size factor [cm]\n')
    f.write('%f\n'%(scale_x))
    f.write('\n')
    f.write('// Factor for output CT value  100000 or  10000   or 1000   or  1\n')
    f.write('%d\n'%(recon_factor))
    f.write('\n')
    f.write('// reconstructed data file for output \n')
    f.write('360pro_FDK_scattermedianepoch900_log_noise0.065per_minx0.1_%s%dy%dh%dw_ICRP110_short.raw\n'%(gender,age,height,weight))
    #f.write('360pro_FDK_%s%dy%dh%dw_ICRP110_short.raw\n'%(gender,age,height,weight))
    f.write('\n')
    f.write('// Prior reconstruction image slices  \n')
    f.write('240\n')
    f.write('\n')
    f.write('//prior image shift (-x, -y, -z) (opposite of Couch Shift direction) \n')
    f.write('0.0 0.0 0.0 \n')
    
        
    f.close()
