function [flux_x_ex, flux_x_im] = set_flux_bc(flux_x_ex, flux_x_im,Hu,u_old,Gu,PHI,g_q)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

flux_x_ex(1) = Hu(2)*u_old(3) - g_q;
flux_x_im(1) = Hu(2).*PHI(2).*Gu(2) - g_q;

end

