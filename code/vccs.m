function vccs(nd1,nd2,ni1,ni2,val)
% vccs.m
% Adds stamp for a voltage-controlled current-source to the global circuit representation
% 
%
%   ni1 -------o+          |----------o nd1
%                          |
%                         / \
%                        / | \    Ind1 to Ind2 = val*(Vni1 - Vni2)
%                Ivcvs   \ | /
%                         \V/ 
%                          |
%   ni2 -------o-          |----------o nd2
%
% Stamp of a voltage controlled current source.
% The dependent nodes are nd1 and nd2 (positive current from nd1 to nd2)
% The independent nodes are ni1 and ni2 (positive voltage at ni1).
% Ind1 to Ind2 = val*(Vni1 - Vni2)
%
%          
% ELEC4506, Lab-2
% Author: Jacob Godin
% Date: 2018/10/01
%--------------------------------------------------------------------------
% define global variables
global G;

if (nd1 ~= 0)
    G(nd1,ni1) = G(nd1,ni1) + val;
    G(nd2,ni2) = G(nd2,ni2) + val;
end
if (nd2 ~= 0)
    G(nd1,ni2) = G(nd1,ni2) - val;
    G(nd2,ni1) = G(nd2,ni1) - val;
end

%END
