function cur(n1,n2,val)
% cur.m:
% Adds stamp for current to the global b-Matrix in circuit representation!
% 
% cur(n1,n2,val):
%                     J=val (Ampere)
%               n1 o---OO---o n2   
%                    ----->
%
% ELEC4506, Lab-2
% Author: Jacob Godin
% Date: 2018/10/01
%--------------------------------------------------------------------------
% define global variables
global b;

if (n1 ~= 0)
    b(n1,1) = -val;
end

if (n2 ~= 0)
    b(n2,1) = val;
end

end