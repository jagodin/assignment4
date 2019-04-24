function ccvs(nd1,nd2,ni1,ni2,val)
% vcvs.m
% Adds stamp for a current-controlled voltage-source to the global circuit representation
% 
%
%   ni1 -------o+          |----------o nd1
%              |           |
%              |          /+\     Vni1 - Vni2 = 0
%              |         /   \    Vnd1 - Vnd2 = Iccvs*G
%        Iccvs |         \   /    Vnd1 - Vnd2 - Iccvs*G = 0
%              V          \ / 
%              |           |
%   ni2 -------o-          |----------o nd2
%
% The nodes across the dependent source are nd1 and nd2 (positive voltage at nd1)
% The independent nodes are ni1 and ni2 (positive voltage at ni1).
%
%          
% ELEC4700, PA9
% Author: Jacob Godin
% Date: 2019/03/19
%--------------------------------------------------------------------------
% define global variables
global G C b;
d = size(G,1); % current size of the MNA
xr1 = d+1;      % two new rows/columns
xr2 = xr1+1;
b(xr2) = 0;     % add two new rows
% Matlab automatically increases the size of a matrix
% if you use an index that is bigger than the current size.


% A(row,col)

G(xr2,xr2) = 0; % add two new rows/columns
C(xr2,xr2) = 0;

G(xr1,xr1) = -val;

if (nd1 ~= 0)
    G(xr1,nd1) = 1;
    G(nd1,xr2) = -1;
end
if (nd2 ~= 0)
    G(xr1,nd2) = -1;
    G(nd2,xr2) = 1;
end

if (ni1 ~= 0)
    G(xr2,ni1) = 1;
    G(ni1,xr1) = 1;
end
if (ni2 ~= 0)
    G(ni2,xr1) = -1;
    G(xr2,ni2) = -1;
end
%END
