%% ** Effective Moduli for Symmetric Laminates **
% Pablo Tejeda - Stephen Burt - Ryan Murphy
% 
function [Ex,Ey, vxy, vyx, Gxy] = EffectiveModuliSymmetricLaminates(A, t, n)
invA = inv(A);
h = t*n;
Ex  = 1/h/invA(1,1);
Ey  = 1/h/invA(1,1);
vxy = -invA(1,2)/invA(1,1);
vyx = -invA(1,2)/invA(2,2);
Gxy = 1/invA(3,3)/h;
end






