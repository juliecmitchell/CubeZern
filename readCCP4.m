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
% readCCP4.m 
%
% This function reads a standard CCP4 .map file and extracts the density
% 'rho' as a 3D array
%
% Inputs:
% - gzname: The name of the source file (including its .gz extension and
% file location)
%
% Outputs:
% - rho: The original density extracted from the MAP file
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rho = readCCP4(gzname)

fprintf('\nExtracting density from the MAP file...\n');
tic; % Starting time

fname = gunzip(gzname);
fname = char(fname);

% Open CCP4 file
fid = fopen(fname,'r');

% Grid size of the map
Nx = fread(fid,1,'long');
Ny = fread(fid,1,'long');
Nz = fread(fid,1,'long');

% Data type (mode) which can be 0, 1 or 2 for CCP4 maps stored in EMDB
data_type = fread(fid,1,'long');

% Go to beginning of the data
fseek(fid, 1024, 'bof');

% If stored as 8-bit signed byte (from -128 to 127)
if data_type == 0
    rho = fread(fid, Nx*Ny*Nz, 'int8');
end

% If stored as 16-bit integer (from -32768 to 32767)
if data_type == 1
    rho = fread(fid, Nx*Ny*Nz, 'int16');
end

% If stored as 32-bit floating point number
if data_type == 2
    rho = fread(fid, Nx*Ny*Nz, 'float32');
end

fclose(fid);
delete(fname);

rho = reshape(rho, Nx, Ny, Nz);

fprintf('Done!\n');
toc; % Finishing time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
