function ind(n1,n2,val)
% ind.m:
% Adds stamp for inductor to the global C-Matrix in circuit representation!
% 
% ind(n1,n2,val):
%                     L=val (H)
%               n1 o---////---o n2   
%
% ELEC4506, Lab-2
% Author: Jacob Godin
% Date: 2018/10/01
%--------------------------------------------------------------------------
% define global variables
global G b C;

d = size(G,1); % current size of the MNA
xr = d+1; % new row
b(xr) = 0; % add new row

G(xr,xr) = 0; % add new row/column
C(xr,xr) = 0; % add new row/column

C(xr, xr) = -val;

if (n1 ~= 0) && (n2 ~= 0)
    G(n1, xr) = 1;
    G(xr, n1) = 1;
    G(n2, xr) = -1;
    G(xr, n2) = -1;
end

if (n2 == 0) || (n1 == 0)
    G(n1, xr) = 1;
    G(xr, n1) = 1;
end

end