
%Based on example 3.2.1 from SWASHES 2013 paper
%This is a subcritical flow with non-trivial bathymetry
%The steady-state water column height is specified, and the bathymetry is solved for
%This bathymetry is then used as input for the numerical solver.

%VARIABLES
%H:   total elevatation of water, H = eta + B
%eta: water elevation above the plane y = 0
%B:   depth of terrain below the plane y = 0
%u:   depth-averaged velocity in x-direction

%EQUATIONS
%Mass: H_t + (H*u)_x = 0
%Momentum: H*(u_t + uu_x) = -g*H*eta_x - gamma*u

addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include')

clear all
close all

%parameters
t_f = 5.0*10^3;
Lx = 1000; %m
g = 9.81; %m/s^2
C_f = (0.033)^2; %m^(-1/3), Manning friction coefficient

%boundary and initial data
H0 = 0.0;
u0 = 0.0;
g_q = 2.0;
bdry = 2;

%Tolerances
n_tol  = 1.0e-6;
wd_tol = 1.0e-5;
max_iter = 10;

%timestepping parameters
u_max = 1; %sqrt(g*eta0) + g_q/g_h;
CFL = 0.06;
theta = 0.5;

%number of cells in each direction
N = 2^8;

%step sizes
dx = Lx/N;
dt = CFL * dx / u_max;
t_incr = 50;
stride = ceil(t_incr/dt);

nu = 0.5*dx^2; %numerical viscosity

%cell centers
xc = ((1:N)' - 0.5) * dx;

%edges 
xe = (0:N)' * dx;
      
%analytical ss solns and bathymetry
[H_act,dH_dx, u_act] = analytical_soln(Lx,xe,g,g_q);
B = solve_bathymetry(H_act,u_act,dH_dx,g,C_f,dx);
eta_act = H_act - B;
g_H = H_act(end);
      
%avg to get a cell centered versions
B   = (B(1:end-1) + B(2:end))/2;
eta_act = (eta_act(1:end-1) + eta_act(2:end))/2;
H_act   = (H_act(1:end-1) + H_act(2:end))/2;
g_eta = eta_act(end);
      
%initialize total water column and surface elevation
H = H0*ones(N,1);
eta = H - B;
      
%axis
eta_axis = [-2*dx Lx+2*dx min(-B) 1.5*max(eta)];
u_axis   = [-2*dx Lx+2*dx -2.5*max(u_act) 1.5*max(u_act)];
q_axis   = [-2*dx Lx+2*dx -2.5*g_q 1.5*g_q];

%initialize bathymetry
B = [B(1); B; B(end)];

%initialize variables
u = u0 * ones(N+1,1);

%create ghost cell arrays
eta_new = [eta(1); eta; g_eta];
H_new = eta_new + B;
u_new = [u(1); u; u(end)];

%intiialize variables for previous time steps
eta_old = eta_new;
H_old = H_new;
u_old = u_new;

u_old_old = u_old;

Hu = calculate_edge_elevations(u_old,H_old);
gamma = zeros(N+1,1);

%total flux
Q_h = Lx * mean(eta);
Q_q = 0;

t_idx = 0;
t = 0;
while t <= t_f
    
    %set ghost cells
    eta_old = set_eta_gc(eta_old,g_eta);
    u_old = set_u_gc(u_old);
        
    %plot
    plot_state(xc,xe,t,t_idx,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Construct u at non-grid points for semi-Lagrangian approach
    u_star = stelling_flux(H_old,Hu,u_old,dx,dt,wd_tol);
    
    %gradient of eta at t_n
    d_eta_n_dx = (eta_old(2:end) - eta_old(1:end-1))/dx;
    
    %Laplacian of u,v at t_n
    Lap_u = Lap_1d(u_star,dx);

    %Casulli's "G" operator
    Gu = construct_Gu(u_star,Lap_u,d_eta_n_dx,g,nu,dt,theta);   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_old = eta_old + B;
    Hu = calculate_edge_elevations(u_star,H_old);
      
    %wet and dry vectors
    w = find(Hu > wd_tol);
    d = find(Hu <= wd_tol);
    
    %friction
    gamma(w) = manning_friction(Hu(w),u_old(w+1),C_f,g);
    gamma(d) = zeros(size(d));
    
    %prefactor used for T and velocity update
    PHI = zeros(N+1,1);
    PHI(w) = Hu(w) ./ (Hu(w) + dt*gamma(w));
    
    %construct V old (linear)
    V_old = H_old * dx;
    
    %construct fluxes
    [flux_ex,flux_im] = construct_lr_fluxes(Hu,Hu,u_old,Gu,PHI);
    [flux_ex,flux_im] = set_flux_bc_alt(flux_ex, flux_im,g_q);
    [flux_diff_ex] = flux_difference(flux_ex);
    [flux_diff_im] = flux_difference(flux_im);
        
    %update volume
    V_new = V_old(2:end-1) - dt*flux_diff_ex;
    
    %reconstruct total elevation H
    H_new(2:end-1) = V_new / dx;
    
    %calculate surface elevation eta
    eta_new(2:N+1) = H_new(2:end-1) - B(2:end-1);
    
    %nonlinear
    T = construct_T(Hu,eta_old,theta,dt,dx,g,PHI);
    b = construct_b(V_old,flux_diff_ex,flux_diff_im,dt,theta);
    [T,b] = add_elevation_bc(Hu,eta_old,PHI,T,b,g,dt,theta,bdry);
    [eta_update,r,iter] = nonlinear_solve_db(xc,dx,eta_old,T,b,n_tol,wd_tol,B,max_iter);
    eta_new(2:end-1) = eta_update;
    
    eta_new = set_eta_gc(eta_new,g_eta);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_eta_n1_dx = (eta_new(2:end) - eta_new(1:end-1))/dx;
    
    %update u
    RHS = Gu - dt*g*theta*d_eta_n1_dx;      
    u_new(2:end-1) = PHI.*RHS;
    
    %set u BC
    u_new = set_u_bc(u_new,dx,g_q,Hu(1));
    u_new(1) = u_old(1);
    u_new(end) = u_old(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update time
    t_idx = t_idx + 1;
    t     = t + dt;
        
    %update previous values to new values
    u_old_old = u_old;
    eta_old = eta_new;
    u_old = u_new;
    
end


%error in SS soln
L_2_e  = [norm(eta_act - eta_new(2:end-1)) norm(u_act - u_new(2:end-1))] * sqrt(dx)
L_oo_e = [max(abs(eta_act - eta_new(2:end-1))) max(abs(u_act - u_new(2:end-1)))]


subplot(1,3,1)
hold on
plot(xc, eta_act)
plot(xc,g_eta*ones(N,1),'--')
plot(xc, -B(2:end-1))
subplot(1,3,2)
hold on
plot(xe, u_act)
subplot(1,3,3)
hold on
plot(xe, g_q*ones(N+1,1))

%read it Clawpack data
fileID = fopen('pyclaw_soln_alt.txt','r');
data = fscanf(fileID,'%f',inf);
row = ceil(length(data) / 3);
data = reshape(data,3,row);
data = data';
row = ceil(row/2);
data = data(1:row,:);
xd = (1:row)*Lx/row;

%plot for export
figure(2)
lw = 3;
fs = 24;
plot(xc,H_new(2:end-1),'LineWidth',lw+1)
hold on
plot(xc,-B(2:end-1),'--','LineWidth',lw)
plot(xc,H_act,'--','LineWidth',lw)
plot(xd,data(:,1),'LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('H','FontSize',fs)
legend('numerical','bathymetry','actual','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)

figure(3)
lw = 3;
fs = 24;
plot(xe,Hu.*u_new(2:end-1),'LineWidth',lw+1)
hold on
plot(xe,g_q*ones(N+1,1),'--','LineWidth',lw)
plot(xd,data(:,2),'LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('Hu','FontSize',fs)
legend('numerical','actual','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)
