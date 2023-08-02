# VirtualCBCT

This project provides a virtual cone-beam CT (CBCT) system with an elementary information for objects. This is composed of two basic codes,

material based forward projection algorithm (MBFPA) for CBCT geometry (GPU version),
Feldkamp-Davis-Kress (FDK) algorithm (CPU version),

and high-resolution ICRP 110 phantoms (both male and female).

These can be used freely, but plese cite the following papars for use.

Ref: Taisei Shimomura, et al., Virtual cone-beam computed tomography simulator with human phantom library and its application to the elemental material decomposition, Physica Medica, 2023

If you have a question, contact by e-mail. haga@tokushima-u.ac.jp

## 1: For use of material based forward projection algorithm (MBFPA) code
### 1-1: Preparation
Go to "Release_CBCT_Projection_CM/" It is required to prepare an object including the elementary information and an X-ray spectrum before runing the code. In this project, the elemental image of the cylindrical water phantoms with 20 cm and 25 cm diameters has been uploaded as the examples (waterphantom/waterphan_20cm and waterphantom/waterphan_25cm, both of which the datatype is the 16-bit unsigned, 512 × 512 x 200). Six elements, H, C, N, O, P, and Ca, are included, but H and O are used for water phantom as well as N for Air.
Only 100 kV spectrum is provided as "Spectrum/result_20220111_average.csv", which is the estimated spectrum in ELEKTA Synergy. One can apply any energy spectrum by modifying this file (first and second rows mean the first and last points in energy bin, and third row means the fraction for photons).

The information of the object and the spectrum must be indicated in the input txt file.
Two templates of the input file were distributed in "txtfile/" folder
(IR_waterphantom_36pro.txt for 36 projection angles with 10 degree interval and IR_waterphantom_360pro.txt for 360 projection angles with 1 degree interval).

The geometry can be varied by modifying parameters listed in physParams.h, where the default values simulate projections in the geometry of Electa Synergy.

### 1-2: How to execute
Prepare elementary density image and incident X-ray spectrum as described above.
"make" produces the execute file "proj.exe". Then, execute
./proj.exe [input.txt]
ex. ./proj.exe txtfile/FDK_36pro_water20cm.txt
This products the reprojection data (the filename is defined in input file).
The datatype of the projection image is a 32-bit real, 400 × 400 (1 mm resolution) in the default geometory.


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
The default datatype of the reconstructed image is a 16-bit short, 270 × 270 (1-mm resolution).
