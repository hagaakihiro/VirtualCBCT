# VirtualCBCT

This project provides a virtual cone-beam CT (CBCT) system with an elementary information for objects. This is composed of two basic codes,

material based forward projection algorithm (MBFPA) for CBCT geometry (GPU version),  
Feldkamp-Davis-Kress (FDK) algorithm (CPU version),

and high-resolution ICRP 110 phantoms (both male and female).
The scattered X-ray simulation with Unet model has beem also be shared;
Unet_scatter.py

These can be used freely, but plese cite the following papar for use.

Ref: Taisei Shimomura, et al., Virtual cone-beam computed tomography simulator with human phantom library and its application to the elemental material decomposition, Physica Medica, 2023

If you have a question, contact by e-mail. haga@tokushima-u.ac.jp

## 1: For use of material based forward projection algorithm (MBFPA) code
### 1-1: Preparation
Go to "Release_CBCT_Projection_CM/" It is required to prepare an object including the elementary information and an X-ray spectrum before runing the code. In this project, the elemental image of the cylindrical water phantoms with 20 cm and 25 cm diameters has been uploaded as the examples (waterphantom/waterphan_20cm and waterphantom/waterphan_25cm, both of which the datatype is the 16-bit unsigned, 512 × 512). Six elements, H, C, N, O, P, and Ca, are included, but H and O only are used for water phantom as well as N for Air.
A 100 kV spectrum is provided in "Spectrum/result_20220111_average.csv", which is the estimated spectrum in ELEKTA Synergy. One can apply any energy spectrum by modifying this file (first and second rows mean the first and last points in energy bin, and third row means the fraction for photons).

The information of the object and the spectrum must be indicated in the input txt file.  
Two templates of the input file were distributed in "txtfile/" folder  
(IR_waterphantom_36pro.txt for 36 projection angles with 10 degree interval and IR_waterphantom_360pro.txt for 360 projection angles with 1 degree interval).  

The geometry can be varied by modifying parameters listed in physParams.h, where the default values simulate projections in the geometry of Electa Synergy.

### 1-2: How to execute
Prepare elementary density image and incident X-ray spectrum as described above.
"make" produces the execute file "proj.exe". Then, execute  
./proj.exe [input.txt]  
ex. ./proj.exe txtfile/FDK_36pro_water20cm.txt  
This products the reprojection data (the output filename is defined in input file).
The datatype of the projection image is a 32-bit real, 400 × 400 (1-mm resolution) in the default geometory.


## 2: For use of FDK reconstruction code
### 2-1: Preparation
Go to "Release_recon_from_forward_pro_float_CM/". It is required for CBCT reconstruction to prepare the projection data produced by above reprojection code. The information of the projection must be indicated in the input txt file, such as "txtfile/FDK_36pro_water20cm.txt" and "txtfile/FDK_360pro_water25cm.txt". Note that the same geometry as that used in the MBFPA should be employed.

### 2-2: How to execute
Prepare projection image  
ex. Release_CBCT_Projection_CM/36pro_raytracing_waterphantom_20cm.raw  
make  
./cbct2 [input.txt]  
ex. ./cbct2 txtfile/FDK_36pro_water20cm.txt  
The output image is producted as "36pro_FDK_waterphantom_20cm_short.raw" (the filename should be indicated in the input txt).
The default datatype of the reconstructed image is a 16-bit short (factored by 1000), 270 × 270 (1-mm resolution).

## 3: High resolution human phantom
### 3-1: Distribution
Because the spatial resolution in the original ICRP 110 phantoms is relatively low (1.775 × 1.775 × 4.84 mm$^3$ for female and 2.137 × 2.137 × 8.0 mm$^3$ for male), both phantoms were enhanced with the higher resolution of 1.0 × 1.0 × 2.0 mm$^3$.
These are found as "ICRPphantoms_2mmslice/AF_wholebody_531x243_aaa_bbb.raw" and "ICRPphantoms_2mmslice/AM_wholebody_543x271_aaa_bbb.raw", where we divided as 9 parts of the whole body, because of the file size limitation in github. The datatype is a 8-bit. The anatomical ID number is assigned in each voxel as well as the original ICRP phantoms. One can marge the files to handle with them as the single human phantom.
### 3-2: Use for the reprojection
The elemental images are needed to create the CBCT reprojections using "Release_CBCT_Projection_CM/" code
(the example format in elemental image can be realized with "waterphantom").
To create such the elemental image data for human phantom, we prepared the program code "To_material_distribution_for512_random.cpp".
For use, first compile this program, and execute it. Then the six major elements (H, C, N, O, P, and Ca) images are created in AM_MD and AF_MD folders (the corresponding density images are also created in AM_ED and AF_ED folders). 

