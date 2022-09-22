%Based on Pastor et al (2008)
%Dam break problem with an extress stress term to mimic non-Newtonian
%behavior, a Bingham fluid
%This case has non-trivial bathymetry, either linear or quadratic.

%VARIABLES
%H:   total elevatation of water, H = eta + B
%eta: water elevation above the plane y = 0
%B:   depth of terrain below the plane y = 0
%u:   depth-averaged velocity in x-direction

%EQUATIONS
%Mass: H_t + (H*u)_x = 0
%Momentum: H*(u_t + uu_x) = -g*H*eta_x - gamma*u - g*H*S_MD,
%where S_MD = tau_MD / (rho*g*H)

addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include')
addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include/adaptive_timestepping')

clear all
close all

%parameters
t_f = 40;%;
Lx = 200; %m
g = 9.81; %m/s^2

rho = 1.0e3; %kg/m^3
tau_y = 1.0e4; %Pa, yield stress
M = 3.0e3; %Pa * s^flow_index, consistency index
C_f = (0.0)^2; %m^(-1/3), Manning friction coefficient


%boundary and initial data
x0 = 0.15*Lx; %100.0;
H0 = 10.0;
u0 = 0.0;
g_q = 0.0;

%Tolerances
n_tol  = 1.0e-10;
wd_tol = 1.0e-2;
max_iter = 10;

%timestepping parameters
u_max = sqrt(g*H0); %sqrt(g*eta0) + g_q/g_h;
CFL = 1.0e-1;
theta = 0.5;

%number of cells in each direction
N = 2^8;

%step sizes
dx = Lx/N;
dt = CFL * dx / u_max;
dt_max = 10*dt;
dt_min = dt/10;
t_incr = t_f/100;
stride = ceil(t_incr/dt);

nu = 0.5*dx^2; %numerical viscosity

%cell centers
xc = ((1:N)' - 0.5) * dx;

%edges 
xe = (0:N)' * dx;

%initialize total water column and surface elevation
H = H0*ones(N,1) .* (xc < x0) + .0*H0*ones(N,1) .* (xc >= x0);
[B,dB_dx] = bathymetry(xc);
eta = H - B;
      
%axis
eta_axis = [-2*dx Lx+2*dx -1.1*max(B) 1.5*max(eta)];
u_axis   = [-2*dx Lx+2*dx -5 50];
q_axis   = [-2*dx Lx+2*dx -10 100];

%initialize bathymetry
B = [B(1); B; B(end)];

%initialize variables
u = u0 * ones(N+1,1);

%create ghost cell arrays
eta_new = [eta(1); eta; eta(end)];
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

max_u = 0;
shear_rate = zeros(N+1,1);

t_idx = 0;
t = 0;
while t <= t_f
    
    %set ghost cells
    eta_old = set_eta_gc(eta_old);
    u_old = set_u_gc(u_old);
    
    %plot
    plot_state(xc,xe,t,t_idx,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Construct u at non-grid points for semi-Lagrangian approach
    u_star = stelling_flux(H_old,Hu,u_old,dx,dt,wd_tol);
    
    %gradient of eta at t_n
    d_eta_n_dx = (eta_old(2:end) - eta_old(1:end-1))/dx;
    
    %Laplacian of u,v at t_n
    Lap_u = Lap_1d(u_old,dx);

    %Casulli's "G" operator
    Gu = construct_Gu(u_star,Lap_u,d_eta_n_dx,g,nu,dt,theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_old = eta_old + B;
    Hu = calculate_edge_elevations(u_old,H_old);
      
    %wet and dry vectors
    w = find(Hu > wd_tol);
    
    %friction
    gamma = zeros(N+1,1);
    gamma(w) = manning_friction(Hu(w),u_old(w+1),C_f,g);
    
    %implicit non-Newtonian slope
    H_inv = zeros(N+1,1);
    H_inv(w) = 1./Hu(w);
    Gu = Gu - dt*tau_y/rho.*H_inv;
    Bingham_coeff = 3*M/rho.*H_inv;
    
    %prefactor used for T and velocity update
    PHI = zeros(N+1,1);
    PHI(w) = Hu(w) ./ (Hu(w) + dt*gamma(w) + dt*Bingham_coeff(w));

    %construct V old (linear)
    V_old = H_old * dx;
    V_old(2:end-1) = total_vol(xc,eta_old(2:end-1),dx,wd_tol);
    
    %construct fluxes
    [flux_ex,flux_im] = construct_lr_fluxes(Hu,Hu,u_old,Gu,PHI);
    [flux_ex,flux_im] = set_flux_bc_alt(flux_ex,flux_im,0);
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
    [eta_update,r,iter] = nonlinear_solve(xc,dx,eta_old,T,b,n_tol,wd_tol,max_iter);
    eta_new(2:end-1) = eta_update;
    
    H_new = eta_new + B;
    eta_new = set_eta_gc(eta_new);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_eta_n1_dx = (eta_new(2:end) - eta_new(1:end-1))/dx;
    
    %update u
    RHS = Gu - dt*g*theta*d_eta_n1_dx;      
    u_new(2:end-1) = PHI.*RHS;
    
    %set u BC
    u_new = set_u_bc(u_new,dx,g_q,Hu(1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update time
    t_idx = t_idx + 1;
    t     = t + dt;
        
    %update previous values to new values
    u_old_old = u_old;
    eta_old = eta_new;
    u_old = u_new;
       
end

%find shock position
xs = dx*find(H_new(2:end-1)>wd_tol*10,1,'last');

%approx of soln near the final front posn
xs_idx = ceil(xs/dx);
[H_act] = analytical_soln(xc,xs,H_new(xs_idx),B,dB_dx,tau_y,rho,g);
eta_act = H_act - B(2:end-1);

figure(1)
subplot(1,3,1)
hold on
plot(xc, eta_act)
plot(xc, -B(2:end-1))
subplot(1,3,2)
hold on
subplot(1,3,3)
hold on
plot(xe, g_q*ones(N+1,1))

%plot for export
figure(2)
lw = 3;
fs = 24;
plot(xc,H_old(2:end-1),'LineWidth',lw+1)
hold on
plot(xc,H_act,'--','LineWidth',lw+1)
plot(xc,-B(2:end-1),'--','LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('H','FontSize',fs)
legend('numerical','actual','bathymetry','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)

figure(3)
lw = 3;
fs = 24;
plot(xe,Hu.*u_new(2:end-1),'LineWidth',lw+1)
hold on
plot(xe,g_q*ones(N+1,1),'--','LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('Hu','FontSize',fs)
legend('numerical','actual','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)
axis(q_axis)
