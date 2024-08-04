% MIT License
%
% Copyright (c) 2012 Julie C. Mitchell and the University of Wisconsin
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.


%
% All of codes were implemented by Atilla Sit.
% MATLAB Release: 7.10 R2010a
%
% If you use these programs, please cite the following reference:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS: Atilla Sit / Julie C Mitchell / George N Phillips, Jr /
% / Stephen J Wright
% TITLE: An Extension of 3D Zernike Moments for Shape Description and
% Retrieval of Maps Defined in Rectangular Solids
% JOURNAL: Molecular Based Mathematical Biology (MBMB), 2012
% http://versita.com/mbmb
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To download the current version of this code, please visit the
% following website:
% http://cubezern.mitchell-lab.org
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% zernike_main.m 
%
% This program reads a standard CCP4 .map file and extracts the map's
% density 'rho'; computes reconstructions of 'rho', of order up to a user
% inputted N, using the basis functions BallZern and CubeZern; computes the
% errors (RMSD and CORR) between the original map 'rho' and its
% reconstructions using BallZern and CubeZern and prints them on screen. It
% also produces comparative plots for the isosurfaces of the given map and
% its reconstructions computed using BallZern and CubeZern, and shows them
% in two different orthogonal views (Z and Y views) for the contour level
% given by the user
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with closing open figures, clearing memory and the command window
close all; clear; clc;

% Ask the user to input a file name for the .map file
fid = -1;
fprintf('\n');
msg = 'Enter the name of the MAP file including its location, e.g. ./maps/emd_5167.map.gz';
while fid < 0
   disp(msg);
   inputs.filename = input('Open file: ', 's');
   [fid,msg] = fopen(inputs.filename,'r');
end
clear fid msg;

% Unzip the .map.gz file, read the CCP4 file and extract the density 'rho',
% the filename 'emdid', the voxel dimensions 'Nx', 'Ny' and 'Nz'
fname = gunzip(inputs.filename);
[~, emdid, ~] = fileparts(char(fname));
clear fname;
emdid = ['EMD-',emdid(end-3:end)];
rho = readCCP4(inputs.filename);
[Nx, Ny, Nz] = size(rho);
rho = reshape(rho, 1, numel(rho));

% Ask the user to input an integer N for the reconstruction order between 0
% and Nmax (inclusive)
Nmax = 20;
fprintf('\nEnter the reconstruction order N where 0 <= N <= %d\n',Nmax);
inputs.N = input ('N = ');
while ( inputs.N < 0 || inputs.N > Nmax || isequal(inputs.N,floor(inputs.N))==0 ) 
  fprintf('Incorrect input, please try again\n');
  fprintf('\nEnter the reconstruction order N where 0 <= N <= %d\n',Nmax);
  inputs.N = input ('N = ');
end

% Compute reconstructions Ballrhohat and Cuberhohat of the given map rho,
% using the basis functions BallZern and CubeZern, respectively, of order
% up to N
minrho = min(rho);
maxrho = max(rho);
[BallMom, Ballrhohat, CubeMom, Cuberhohat] = reconstruction(rho, Nx, Ny, Nz, Nmax, inputs.N);

% Measure the errors (RMSD and CORR) between the original map and its
% reconstructions using BallZern and CubeZern and print them on screen
BallRMSD = sqrt( sum((rho - Ballrhohat).^2) / numel(rho) );
CubeRMSD = sqrt( sum((rho - Cuberhohat).^2) / numel(rho) );
fprintf('\nRMSD values between the original map and the reconstruction obtained..\n')
fprintf('..by using BallZern = %2.4f\n',BallRMSD)
fprintf('..by using CubeZern = %2.4f\n',CubeRMSD)
temp = corrcoef(rho,Ballrhohat);
BallCORR = temp(1,2);
temp = corrcoef(rho,Cuberhohat);
CubeCORR = temp(1,2);
clear temp;
fprintf('\nCORR values between the original map and the reconstruction obtained..\n')
fprintf('..by using BallZern = %2.4f\n',BallCORR)
fprintf('..by using CubeZern = %2.4f\n',CubeCORR)

% Ask the user to input the contour level for the isosurfaces suggested by
% the EMDB
fprintf('\nEnter the contour level for the map %s as suggested by the EMDB\n',emdid);
inputs.contourlevel = input ('Contour level = ');
while ( inputs.contourlevel < minrho || inputs.contourlevel > maxrho ) 
  fprintf('Incorrect input, please enter a contour level between %4.5f and %4.5f..\n',minrho,maxrho);
  fprintf('..or the one suggested by the EMDB\n');
  inputs.contourlevel = input ('Contour level = ');
end

% Produce comparative plots for the isosurfaces of the given map and its
% reconstructions computed using BallZern and CubeZern, and show them in
% two different orthogonal views (Z and Y views) for the contour level
% inputted by the user
[rho, Ballrhohat, Cuberhohat] = isosurfaces(rho, Ballrhohat, Cuberhohat, Nx,Ny,Nz, inputs, emdid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
