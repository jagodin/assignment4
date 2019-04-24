function vol(n1,n2,val)
% vol.m:
% Add stamp for voltage source to the global circuit representation
%            
%               Vsrc=val (Volts)
%       n1 o-----(+  -)-----o n2
%               ------->
%                 Isrc
%         
% ELEC4506, Lab-2
% Author:
% Date:
%--------------------------------------------------------------------------
% define global variables
global G b C;

d = size(G,1); % current size of the MNA
xr = d+1; % new row
b(xr) = 0; % add new row
% Matlab automatically increases the size of a matrix
% if you use an index that is bigger than the current size.
G(xr,xr) = 0; % add new row/column
C(xr,xr) = 0; % add new row/column

if (n1 ~= 0)
    G(n1,xr) = 1;
    G(xr,n1) = 1;
end

if (n2 ~= 0)
    G(n2,xr) = -1;
    G(xr,n2) = -1;
end
b(xr) = val;
%End
