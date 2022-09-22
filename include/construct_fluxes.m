function [flux_x_ex,flux_x_im] = construct_fluxes(Hu,u_old,Gu,PHI)
%CONSTRUCT_FLUXES Construct fluxes for the mass conservation equation
%   flux_x_ex is the explicit flux in x direction
%   flux_x_im is part of the implicit flux

flux_x_ex = Hu(2:end)            .*u_old(3:end-1) - Hu(1:end-1)              .*u_old(2:end-2);
flux_x_im = Hu(2:end).*PHI(2:end).*Gu(2:end)      - Hu(1:end-1).*PHI(1:end-1).*Gu(1:end-1);

end

