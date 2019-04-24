function cap(n1,n2,val)
% cap.m:
% Adds stamp for capacitor to the global C-Matrix in circuit representation!
% 
% cap(n1,n2,val):
%                     C=val (Farad)
%               n1 o---||---o n2   
%
% ELEC4506, Lab-2
% Author:
% Date:
%--------------------------------------------------------------------------
% define global variables
global C;

if (n1 ~= 0)
    C(n1,n1) = C(n1,n1) + val;
end

if (n2 ~= 0)
    C(n2,n2) = C(n2,n2) + val;
end

if (n1 ~= 0) && (n2 ~= 0)
    C(n1,n2) = C(n1,n2) - val;
    C(n2,n1) = C(n2,n1) - val;
end
%END
