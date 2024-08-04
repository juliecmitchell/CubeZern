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
% cubezern.m
%
% This program computes the new Zernike polynomials numerically that are
% orthonormal within the cube inscribed inside the unit ball up to the
% expansion order Nmax; and stores them in the file
% 'cubezern_upto_N=Nmax.mat' (See Section 3 in the article for related
% formulas)

% Dimensions of the rectangular prism inscribed in the unit ball
a = 1/sqrt(3);
b = 1/sqrt(3);
c = 1/sqrt(3);  % We use only the cube in this work

Nmax = 20;
Kmax = (Nmax+1)*(Nmax+2)*(Nmax+3)/6;

filename = sprintf('ballzern_upto_N=%d.mat',Nmax);
load(filename);
filename = sprintf('exponents_upto_N=%d.mat',Nmax);
load(filename);
filename = sprintf('indices_upto_N=%d.mat',Nmax);
load(filename);

% A precomputation to speed up the Gram-Schmidt process (see equation (43)
% in the article)

A = zeros(Kmax,Kmax);
for tau = 1:Kmax
    rtau = exponents(tau,1);
    stau = exponents(tau,2);
    ttau = exponents(tau,3);
    for zeta = 1:Kmax
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

% Use the nonrecursive matrix approach to compute CubeZern of order up to
% N = Nmax

% First, orthonormalize the polynomials with respect to the standard inner
% product <P,Q> = P'*Q to ensure that the monomial coefficients in BallZern
% are controlled (i.e. they do not grow rapidly)
C = BallZern' * BallZern;
Q = chol(C);
M = inv(Q');
NewBallZern = zeros(Kmax,Kmax);
for mu = 1:Kmax
    for nu = 1:mu
        NewBallZern(:,mu) = NewBallZern(:,mu) + M(mu,nu) * BallZern(:,nu);
    end
end

% Now use the methodology explained in Section 3.2 of the article to
% orthonormalize the polynomials with respect to the inner product
% <P,Q> = (1/(8*a*b*c)) * int_[-a,a] int_[-b,b] int_[-c,c] {P(x,y,z)*conj(Q(x,y,z))} dz dy dx
C = NewBallZern' * A * NewBallZern;
Q = chol(C);
M = inv(Q');
CubeZern = zeros(Kmax,Kmax);
for mu = 1:Kmax
    for nu = 1:mu
        CubeZern(:,mu) = CubeZern(:,mu) + M(mu,nu) * NewBallZern(:,nu);
    end
end

% Save the Zernike polynomials in cube of order up to Nmax into
% cubezern_upto_N=Nmax.mat
filename = sprintf('cubezern_upto_N=%d.mat',Nmax);
save(filename, 'CubeZern');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Check one by one
% syms('x','real');
% syms('y','real');
% syms('z','real');
% syms temp;
% K = 20;
% for J = 1:K
%     temp = 0;
%     for I = 1:K
%         temp = temp + CubeZern(I,J)*x^(exponents(I,1))*y^(exponents(I,2))*z^(exponents(I,3));
%     end
%     simplify(temp)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
