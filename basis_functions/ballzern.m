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
% ballzern.m
%
% This program computes the classical Zernike polynomials numerically that
% are orthonormal within the unit ball, up to the expansion order Nmax and
% stores them in the file 'ballzern_upto_N=Nmax.mat' (See Section 2 in the
% article for related formulas)

Nmax = 20;

% Initial variables
N = 0;
exponents = [0 0 0];  % list of exponents r,s,t
indices = [0 0 0];  % list of indices n,l,m
BallZern = 1;  % The very first Zernike polynomial: Z_{1} = Z_{000} = 1
Kprev = 1;  % Current number of Zernike polynomials

% Note that {The set of Zernike polynomials of order up to N} = {The set of
% Zernike polynomials of order up to N-1} U {The set of Zernike polynomials
% of order = N}, so we can only compute the polynomials of order = N, and
% add those to the set of polynomials computed for N-1.

for N = 1:Nmax
    
    % Add the new exponents r,s,t in {x^r}*{y^s}*{z^t} to the list
    % 'exponents', where r,s,t satisfy r+s+t = N
    for r = 0:N
        for s = 0:N
            for t = 0:N
                if (r+s+t == N)
                    exponents = [exponents; r s t];
                end
            end
        end
    end
    
    % Add the new indices n,l,m for the Zernike functions Z_{nlm} to the
    % list of 'indices', where n = N, n-l is even, l = 0,...,n, and 
    % m = -l,...,l
    for l = 0:N
        if mod(N-l,2) == 0
            indices = [indices; N l 0];
            for m = 1:l
                indices = [indices; N l  m];
                indices = [indices; N l -m];
            end
        end
    end
    
    % For a given N, total number of monomials {x^r}*{y^s}*{z^t} such that
    % 0 <= r+s+t <= N is equal to K = (N+1)*(N+2)*(N+3)/6; and the total
    % number of Zernike functions of order up to N is also K, so BallZern
    % will be a K-by-K matrix for each N.
    
    K = (N+1)*(N+2)*(N+3)/6;
    
    % Update BallZern by embedding the previously computed set of Zernike
    % polynomials of order up to N-1 to new set
    Temp = zeros(K,K);
    Temp(1:Kprev,1:Kprev) = BallZern;
    BallZern = Temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Next loop computes the set of polynomials of order 'equal to N' using
    % the equations (10) and (11) in the article, and updates BallZern
    
    for J = Kprev+1:K
        
        n = indices(J,1);
        l = indices(J,2);
        m = indices(J,3);
        
        if m >= 0
            
            lminusm = l - m;
            k = (n-l)/2;
            twok = 2*k;
            coeff1 = (-1)^(k+m) * sqrt((2*l+1)*(2*l+4*k+3)/3) * sqrt(factorial(l+m)) * sqrt(factorial(lminusm)) * (factorial(m) / 2^(m+twok)) ;
            kplusl = k + l;
            twokplusl = 2*kplusl;
            ubformu = floor( lminusm / 2 );
            t1 = l - m;
            
            for nu = 0:k
                twonu = 2*nu;
                lplusnu = l + nu;
                kpluslplusnu = kplusl + nu;
                coeff2 = ( coeff1 * (-1)^nu * prod( 2*lplusnu+2 : 2*kpluslplusnu+1 ) / factorial(k-nu) ) / prod( lplusnu+1 : kpluslplusnu ) ;
                t2 = t1 + twonu;
                
                for alpha = 0:nu
                    coeff3 = coeff2 / factorial(alpha);
                    numinusalpha = nu - alpha;
                    r1 = 2*alpha;
                    t3 = t2 - r1;
                    
                    for beta = 0:numinusalpha
                        coeff4 = ( coeff3 / factorial(beta) ) / factorial(numinusalpha-beta) ;
                        twobeta = 2*beta;
                        s2 = m + twobeta;
                        t4 = t3 - twobeta;
                        
                        for u = 0:m
                            coeff5 = ( coeff4 * (-1i)^u / factorial(u) ) / factorial(m-u) ;
                            r2 = r1 + u;
                            s3 = s2 - u;
                            
                            for mu = 0:ubformu
                                twomu = 2*mu;
                                coeff6 = ( ( coeff5 * (-1)^mu /  2^(twomu) ) / factorial(m+mu) ) / factorial(lminusm-twomu) ;
                                s4 = s3 + twomu;
                                t  = t4 - twomu;
                                
                                for tau = 0:mu
                                    r = r2 + 2*tau;
                                    s = s4 - 2*tau;
                                    
                                    for I = 1:K
                                        
                                        if 1 == isequal( [r s t] , exponents(I,1:3) )
                                            coeff = ( coeff6 / factorial(tau) ) / factorial(mu-tau) ;
                                            BallZern(I,J) = BallZern(I,J) + coeff;
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                end
            end
            
        end
        
    end
    
    % Calculate the Zernike polynomials with negative m from the ones with
    % positive m using the symmetry relation given in the article
    for J = 1:K
        m = indices(J,3);
        if m < 0                               
            BallZern(:,J) =  (-1)^(abs(m)) * conj(BallZern(:,J-1)) ;
        end
    end
    
    Kprev = K;
    
end

% Save the Zernike polynomials in ball of order up to Nmax into
% ballzern_upto_N=Nmax.mat
filename = sprintf('ballzern_upto_N=%d.mat',Nmax);
save(filename, 'BallZern');
filename = sprintf('exponents_upto_N=%d.mat',Nmax);
save(filename, 'exponents');
filename = sprintf('indices_upto_N=%d.mat',Nmax);
save(filename, 'indices');

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
%         temp = temp + BallZern(I,J)*x^(exponents(I,1))*y^(exponents(I,2))*z^(exponents(I,3));
%     end
%     simplify(temp)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
