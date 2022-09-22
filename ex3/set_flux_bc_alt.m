function [flux_x_ex, flux_x_im] = set_flux_bc_alt(flux_x_ex, flux_x_im,g_q1,g_q2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%set boundary fluxes with flow rates (e.g. flux = g_q)
flux_x_ex(1) = g_q1;
flux_x_im(1) = g_q1;

% flux_x_ex(end) = g_q2;
% flux_x_im(end) = g_q2;

end