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
% reconstruction.m
%
% This function computes reconstructions of a given density map, of order
% up to N, using the basis functions BallZern and CubeZern
%
% Inputs:
% - rho: The original density extracted from the MAP file
% - Nx, Ny, Nz: Voxel dimensions
% - Nmax: The maximal degree of polynomials which have been computed and
% stored in the folder ./basis_functions/ [ Nmax = 20 in this work ]
% - N: The order of polynomials, up to which the reconstructions will be
% computed [ 0 <= N <= Nmax in this case ]
%
% Outputs:
% - BallMom: Zernike moments for BallZern (coefficients in the expansion)
% - Ballrhohat: Reconstructed density computed using BallZern
% - CubeMom: Zernike moments for CubeZern (coefficients in the expansion)
% - Cuberhohat: Reconstructed density computed using CubeZern
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BallMom, Ballrhohat, CubeMom, Cuberhohat] = reconstruction(rho, Nx, Ny, Nz, Nmax, N)

% Form grids in the cube inscribed inside the unit ball
fprintf('\nForming the grid points...\n');
tic;

L = (N+1)*(N+2)/2;
firstTwoExponents = zeros(L,2);
l = 0;
for n = 0:N
    for r = 0:n
        for s = 0:n
            if (r+s == n)
                l = l + 1;
            	firstTwoExponents(l,:) = [r s];
            end
        end
    end
end

filename = sprintf('basis_functions/exponents_upto_N=%d.mat',Nmax);
load(filename);

K = (N+1)*(N+2)*(N+3)/6;
exponents = exponents(1:K,:);

exponents = [exponents, zeros(K,1)];
for I = 1:K
    if exponents(I,4) == 0
        temp1 = [exponents(I,1) exponents(I,2)];
        for l = 1:L
            temp2 = [firstTwoExponents(l,1)  firstTwoExponents(l,2)];
            if 1 == isequal(temp1, temp2)
                break;
            end
        end
        for J = I:K
            temp2 = [exponents(J,1) exponents(J,2)];
            if 1 == isequal(temp1, temp2)
                exponents(J,4) = l;
            end
        end
    end
end

aa = 1/sqrt(3);   bb = 1/sqrt(3);   cc = 1/sqrt(3);
dx = (aa-(-aa))/Nx;   dy = (bb-(-bb))/Ny;   dz = (cc-(-cc))/Nz;
xc = -aa:dx:aa;   yc = -bb:dy:bb;   zc = -cc:dz:cc;

xmid = 0.5*( xc(1:Nx) + xc(2:Nx+1) );
ymid = 0.5*( yc(1:Ny) + yc(2:Ny+1) );
zmid = 0.5*( zc(1:Nz) + zc(2:Nz+1) );

x = zeros(1,Nx);  y = zeros(1,Ny);
xMidyMid = zeros(Nx,Ny,L);
for l = 1:L
    r = firstTwoExponents(l,1);  s = firstTwoExponents(l,2);
    x(l,:) = (xc(2:Nx+1).^(r+1) - xc(1:Nx).^(r+1)) / (r+1);
    y(l,:) = (yc(2:Ny+1).^(s+1) - yc(1:Ny).^(s+1)) / (s+1);
    xMidyMid(:,:,l) = (xmid.^r)' * (ymid.^s);
end
clear xc yc xmid ymid firstTwoExponents;

z = zeros(K,Nz);
zMid = zeros(K,Nz);
for I = 1:K    
    t = exponents(I,3);
    z(I,:) = (zc(2:Nz+1).^(t+1) - zc(1:Nz).^(t+1)) / (t+1);
    zMid(I,:) = zmid.^t;
end
clear zc zmid;

fprintf('Done!\n');
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modify the density so that it has zero mean 
origmean = mean(rho);
rho = rho - origmean;
rho = reshape(rho,Nx,Ny,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the geometric moments
fprintf('\nComputing the geometric moments M...\n');
tic;

xyrho = zeros(L,Nz);
for l = 1:L
    xy = x(l,:)' * y(l,:);
    for k = 1:Nz
        xyrho(l,k) = sum(sum( xy .* rho(:,:,k) ));
    end
end
clear xy rho;

M = zeros(K,1);
for I = 1:K
	M(I) = sum( z(I,:) .* xyrho(exponents(I,4),:) );
end
clear xyrho z;

fprintf('Done!\n');
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Zernike moments for BallZern (coefficients in the expansion)
fprintf('\nComputing Zernike moments for BallZern of order up to %d...\n',N);
tic;
filename = sprintf('basis_functions/ballzern_upto_N=%d.mat',Nmax);
load(filename);
BallZern = BallZern(1:K,1:K);
BallMom = (3/(4*pi)) * BallZern' * M;
fprintf('Done!\n');
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing the reconstruction density using BallZern
fprintf('\nComputing the reconstruction density using BallZern of order up to %d...\n',N);
tic;

BallZern = BallZern * BallMom;

BzMid = zeros(l,Nz);
for I = 1:K
    BzMid(exponents(I,4),:) = BzMid(exponents(I,4),:) + BallZern(I) * zMid(I,:);
end
clear BallZern;

Ballrhohat = zeros(Nx,Ny,Nz);
for k = 1:Nz
    temp = zeros(Nx,Ny);
    for l = 1:L
        temp = temp + BzMid(l,k) * xMidyMid(:,:,l);
    end
    Ballrhohat(:,:,k) = temp;
end
clear BzMid temp;

Ballrhohat = reshape(real(Ballrhohat),1,numel(Ballrhohat));
Ballrhohat = Ballrhohat - mean(Ballrhohat) + origmean;
fprintf('Done!\n');
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Zernike moments for CubeZern (coefficients in the expansion)
fprintf('\nComputing Zernike moments for CubeZern of order up to %d...\n',N);
tic;
filename = sprintf('basis_functions/cubezern_upto_N=%d.mat',Nmax);
load(filename);
CubeZern = CubeZern(1:K,1:K);
CubeMom = (1/(8*aa*bb*cc)) * CubeZern' * M;
fprintf('Done!\n');
toc;
clear M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing the reconstruction density using CubeZern
fprintf('\nComputing the reconstruction density using CubeZern of order up to %d...\n',N);
tic;

CubeZern = CubeZern * CubeMom;

CzMid = zeros(l,Nz);
for I = 1:K
    CzMid(exponents(I,4),:) = CzMid(exponents(I,4),:) + CubeZern(I) * zMid(I,:);
end
clear CubeZern zMid exponents;

Cuberhohat = zeros(Nx,Ny,Nz);
for k = 1:Nz
    temp = zeros(Nx,Ny);
    for l = 1:L
        temp = temp + CzMid(l,k) * xMidyMid(:,:,l);
    end
    Cuberhohat(:,:,k) = temp;
end
clear CzMid temp xMidyMid;

Cuberhohat = reshape(real(Cuberhohat),1,numel(Cuberhohat));
Cuberhohat = Cuberhohat - mean(Cuberhohat) + origmean;
fprintf('Done!\n');
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
