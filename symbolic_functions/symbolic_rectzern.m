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
% symbolic_rectzern.m
%
% This program generates symbolic expessions of Zernike polynomials
% orthonormal within the rectangular prism inscribed inside the unit ball
% for Table II in the article (requires MATLAB's Symbolic Math Toolbox)

% Start with clearing memory and the command window
clear; clc;

Nmax = 20;
filename = sprintf('../basis_functions/ballzern_upto_N=%d.mat',Nmax);
load(filename);

N = 2;
K = (N+1)*(N+2)*(N+3)/6;

RectZern = sym( BallZern(1:K,1:K) );
exponents = exponents(1:K,:);
indices = indices(1:K,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute symbolic expressions for the Zernike functions orthonormal within
% the rectangular prism inscribed inside the unit ball with dimensions
% a:b:c

% Dimensions of the rectangular prism
syms('a','positive');
syms('b','positive');
syms('c','positive');
c = sqrt(1-(a^2+b^2));

% A precomputation to speed up the Gram-Schmidt process (see equation (43)
% in the article)
A = sym(zeros(K,K));
for tau = 1:K
    rtau = exponents(tau,1);
    stau = exponents(tau,2);
    ttau = exponents(tau,3);
    for zeta = 1:K
        rtauplusrzeta = rtau + exponents(zeta,1);
        if mod(rtauplusrzeta,2) == 0
            stauplusszeta = stau + exponents(zeta,2);
            if mod(stauplusszeta,2) == 0
                ttauplustzeta = ttau + exponents(zeta,3);
                if mod(ttauplustzeta,2) == 0
                	A(tau,zeta)  = ( (a^rtauplusrzeta) * (b^stauplusszeta) * (c^ttauplustzeta) ) ...
                        / ( (rtauplusrzeta+1) * (stauplusszeta+1) * (ttauplustzeta+1) );
                end
            end
        end
    end
end

% Recursive Gram-Schmidt process
for J = 2:K
    for I = 1:J-1
        coef = simplify( RectZern(:,J)' * simplify( A * RectZern(:,I) ) );
        RectZern(:,J) = simplify( RectZern(:,J) - coef * RectZern(:,I) );
    end
    nrmsq = simplify( RectZern(:,J)' * simplify( A * RectZern(:,J) ) );
    nrm = sqrt(nrmsq);
    RectZern(:,J) = simplify( RectZern(:,J) / nrm );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write the symbolic Zernike functions CubeZern on screen in pretty form
% Note that some of the expressions can manually be simplied further

syms('x','real');
syms('y','real');
syms('z','real');

for J = 1:K
    temp = 0;
    for I = 1:K
        r = exponents(I,1);   s = exponents(I,2);   t = exponents(I,3);
        temp = simplify( temp + RectZern(I,J) * x^r * y^s * z^t );
    end
    R = simplify(temp);
    n = indices(J,1);   l = indices(J,2);   m = indices(J,3);
    fprintf('\n\nR_{%d} = R_{%d%d}^{%d} = \n',J,n,l,m)
    pretty(R)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
