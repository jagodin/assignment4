function vcvs(nd1,nd2,ni1,ni2,val)
% vcvs.m
% Adds stamp for a voltage-controlled voltage-source to the global circuit representation
% 
%
%   ni1 -------o+          |----------o nd1
%                          |
%                         /+\
%                      | /   \    Vnd1 - Vnd2 = val*(Vni1 - Vni2)
%                Ivcvs | \   /
%                      V  \-/ 
%                          |
%   ni2 -------o-          |----------o nd2
%
% The nodes across the dependent source are nd1 and nd2 (positive voltage at nd1)
% The independent nodes are ni1 and ni2 (positive voltage at ni1).
%   Vnd1 - Vnd2 = val*(Vni1 - Vni2)
%
%          
% ELEC4506, Lab-2
% Author:
% Date:
%--------------------------------------------------------------------------
% define global variables
global G C b;
d = size(G,1); % current size of the MNA
xr = d+1;      % new row/column
b(xr) = 0;     % add new row
% Matlab automatically increases the size of a matrix
% if you use an index that is bigger than the current size.
G(xr,xr) = 0; % add new row/column
C(xr,xr) = 0;

if (nd1 ~= 0)
    G(nd1,xr) = 1;
    G(xr,nd1) = 1;
end
if (nd2 ~= 0)
    G(nd2,xr) = -1;
    G(xr,nd2) = -1;
end

if (ni1 ~= 0)
    G(xr,ni1) = G(xr,ni1)-val;
end
if (ni2 ~= 0)
    G(xr,ni2) = G(xr,ni2)+val;
end
%END
