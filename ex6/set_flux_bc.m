function [flux_x_ex, flux_x_im] = set_flux_bc(flux_x_ex, flux_x_im,Hu,u_old,Gu,PHI,g_u)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

flux_x_ex(1) = Hu(2)*u_old(3) - Hu(1)*g_u;
flux_x_im(1) = Hu(2).*PHI(2).*Gu(2) - Hu(1)*PHI(1)*g_u;

end

