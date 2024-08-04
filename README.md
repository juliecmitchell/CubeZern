# CubeZern

**************************************************************************************

CubeZern - (c) 2012

**************************************************************************************

This program or any other programs supplied with it are free to use/share for 
non-commerical purpose only.

All of this product's content was implemented by Atilla Sit.
MATLAB Release: 7.10 R2010a

If you use these programs, please cite the following reference:

**************************************************************************************

AUTHORS: Atilla Sit / Julie C Mitchell / George N Phillips, Jr / Stephen J Wright
TITLE: An Extension of 3D Zernike Moments for Shape Description and Retrieval of Maps 
Defined in Rectangular Solids
JOURNAL: Molecular Based Mathematical Biology (MBMB), 2012
https://doi.org/10.2478/mlbmb-2013-0004

**************************************************************************************

The CubeZern package contains the following files:

/CubeZern/README.txt
/CubeZern/zernike_main.m --> The central file of the product
/CubeZern/isosurfaces.m
/CubeZern/readCCP4.m
/CubeZern/reconstruction.m
/CubeZern/sample.png

/CubeZern/basis_functions/ballzern.m
/CubeZern/basis_functions/cubezern.m
/CubeZern/basis_functions/ballzern_upto_N=20.mat
/CubeZern/basis_functions/cubezern_upto_N=20.mat
/CubeZern/basis_functions/exponents_upto_N=20.mat
/CubeZern/basis_functions/indices_upto_N=20.mat

/CubeZern/symbolic_functions/symbolic_ballzern.m  --> Requires Symbolic Math Toolbox
/CubeZern/symbolic_functions/symbolic_cubezern.m  --> Requires Symbolic Math Toolbox
/CubeZern/symbolic_functions/symbolic_rectzern.m  --> Requires Symbolic Math Toolbox

/CubeZern/maps/emd_1148.map.gz
/CubeZern/maps/emd_2167.map.gz
/CubeZern/maps/emd_2227.map.gz
/CubeZern/maps/emd_5167.map.gz
/CubeZern/maps/emd_5381.map.gz

**************************************************************************************

To run the program, change the MATLAB current directory to /CubeZern/, run the 
central file "zernike_main.m", and input the requested information on the command window

SAMPLE USAGE:

On Matlab's command window, type

>> zernike_main

Enter the name of the MAP file including its location, e.g. ./maps/emd_5167.map.gz
Open file: ./maps/emd_5167.map.gz  <-- SAMPLE USER INPUT

Extracting density from the MAP file...
Done!
Elapsed time is 0.041096 seconds.

Enter the reconstruction order N where 0 <= N <= 20
N = 20  <-- SAMPLE USER INPUT

Forming the grid points...
Done!
Elapsed time is 0.578482 seconds.

Computing the geometric moments M...
Done!
Elapsed time is 0.077772 seconds.

Computing Zernike moments for BallZern of order up to 20...
Done!
Elapsed time is 0.323883 seconds.

Computing the reconstruction density using BallZern of order up to 20...
Done!
Elapsed time is 0.096080 seconds.

Computing Zernike moments for CubeZern of order up to 20...
Done!
Elapsed time is 0.355402 seconds.

Computing the reconstruction density using CubeZern of order up to 20...
Done!
Elapsed time is 0.155040 seconds.

RMSD values between the original map and the reconstruction obtained..
..by using BallZern = 0.0098
..by using CubeZern = 0.0077

CORR values between the original map and the reconstruction obtained..
..by using BallZern = 0.8637
..by using CubeZern = 0.9166

Enter the contour level for the map EMD-5167 as suggested by the EMDB
Contour level = 0.0315  <-- SAMPLE USER INPUT

Plotting isosurfaces for the map EMD-5167 and its reconstructions...
Done!

(To view a sample plot, open the file sample.png)

**************************************************************************************
