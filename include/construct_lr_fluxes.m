function [flux_ex,flux_im] = construct_lr_fluxes(Hu_ex,Hu_im,u_old,Gu,PHI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

flux_ex = Hu_ex.*u_old(2:end-1);
flux_im = Hu_im.*PHI.*Gu;

end