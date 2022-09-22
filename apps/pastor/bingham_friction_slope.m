function [S_MD] = bingham_friction_slope(Hu,u_old,g,rho,tau_y,K,wd_tol)
%BINGHAM_FRICTION_SLOPE Friction term for a Bingham fluid stress
%   Approximates the stress of a Bingham fluid via a stress term. tau_y
%   is the yield stress and K is the kinematic viscosity/consistency index
%   (consistency index is used for Herschel Bulkley models). For chanel
%   flow/SVE eqns we replace Hu with R, the hydraulic radius

F = 1./Hu;

%shear rate (approx. assumes a parabolic profile in vertical)
shear_rate = 3*u_old.*F;

%add yield stress and shear stress
tau = tau_y + K*shear_rate;

%friction slope (to be added to momentum eqn u_t + uu_x =...)
S_MD = tau ./ (rho*g*Hu);

end

