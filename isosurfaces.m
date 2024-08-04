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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% isosurfaces.m
%
% This function produces comparative plots for the isosurfaces of the
% given map and its reconstructions computed using BallZern and CubeZern,
% and shows them in two different orthogonal views (Z and Y views) for the
% contour level inputted by the user
%
% Inputs:
% - rho: The original density extracted from the MAP file
% - Ballrhohat: Reconstructed density computed using BallZern
% - Cuberhohat: Reconstructed density computed using CubeZern
% - Nx, Ny, Nz: Voxel dimensions
% - emdid: EMDB ID for the map
% - inputs: Structure of inputs
%      - filename: The filename (including its location) inputted by the
%      user
%      - N: The order of polynomials, up to which the reconstructions will
%      be computed [ 0 <= N <= Nmax in this case ]
%      - contourlevel: The contour level for which the isosurfaces will be
%      plotted
%
% Outputs:
% - rho: The original density extracted from the MAP file (Reshaped to a 3D
% array)
% - Ballrhohat: Reconstructed density computed using BallZern (Reshaped to
% a 3D array)
% - Cuberhohat: Reconstructed density computed using CubeZern (Reshaped to
% a 3D array)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho, Ballrhohat, Cuberhohat] ...
    = isosurfaces(rho, Ballrhohat, Cuberhohat, Nx, Ny, Nz, inputs, emdid)

fprintf('\nPlotting isosurfaces for the map %s and its reconstructions...\n',emdid);

a = Nx;  b = Ny;  c = Nz;

dx = (a-(-a))/Nx;  x = (-a:dx:a) + dx/2;  x(end) = [];  lbx = min(x);  ubx = max(x);
dy = (b-(-b))/Ny;  y = (-b:dy:b) + dy/2;  y(end) = [];  lby = min(y);  uby = max(y);
dz = (c-(-c))/Nz;  z = (-c:dz:c) + dz/2;  z(end) = [];  lbz = min(z);  ubz = max(z);
[x,y,z] = ndgrid(x,y,z);

rho = reshape(rho,Nx,Ny,Nz);
Ballrhohat = reshape(Ballrhohat,Nx,Ny,Nz);
Cuberhohat = reshape(Cuberhohat,Nx,Ny,Nz);

BallTitle = sprintf('BallZern  N=%d',inputs.N);
CubeTitle = sprintf('CubeZern  N=%d',inputs.N);

figure;

subplot(2,3,1); title({emdid;'Z-view'},'fontweight','b');
fv = isosurface(x,y,z,rho,inputs.contourlevel);
patch(fv,'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
view(0,90); daspect([1 1 1]);
axis([lbx ubx lby uby lbz ubz -Inf Inf])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
camlight; lighting phong;

Ballfv = isosurface(x,y,z,Ballrhohat,inputs.contourlevel);
subplot(2,3,2); title({BallTitle,'Z-view'},'fontweight','b');
patch(Ballfv,'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
view(0,90); daspect([1 1 1]);
axis([lbx ubx lby uby lbz ubz -Inf Inf])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
camlight; lighting phong;

Cubefv = isosurface(x,y,z,Cuberhohat,inputs.contourlevel);
subplot(2,3,3); title({CubeTitle,'Z-view'},'fontweight','b');
patch(Cubefv,'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
view(0,90); daspect([1 1 1]);
axis([lbx ubx lby uby lbz ubz -Inf Inf])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
camlight; lighting phong;

subplot(2,3,4); title({emdid;'Y-view'},'fontweight','b');
patch(fv,'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
view(180,0); daspect([1 1 1]); camroll(-90);
axis([lbx ubx lby uby lbz ubz -Inf Inf])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
camlight; lighting phong;

subplot(2,3,5);  title({BallTitle,'Y-view'},'fontweight','b');
patch(Ballfv,'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
view(180,0); daspect([1 1 1]); camroll(-90);
axis([lbx ubx lby uby lbz ubz -Inf Inf])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
camlight; lighting phong;

subplot(2,3,6); title({CubeTitle,'Y-view'},'fontweight','b');
patch(Cubefv,'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
view(180,0); daspect([1 1 1]); camroll(-90);
axis([lbx ubx lby uby lbz ubz -Inf Inf])
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
camlight; lighting phong;

fprintf('Done!\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
