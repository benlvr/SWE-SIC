
%Based on example 1 from Savant et al (2011?)
%This is an immersed bump with no diffusion, friction

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
t_f = 200;
Lx = 25; %m
g = 9.81; %m/s^2

%boundary and initial data
eta0 = 0.33;
u0 = 0.0;
g_h = eta0;
g_q = 0.18;
bdry = 2; %boundary where g_h is specified

%Tolerances
n_tol  = 1.0e-6;
wd_tol = 1.0e-10;
max_iter = 10;

%timestepping parameters
u_max = sqrt(g*eta0) + g_q/g_h;
CFL = 0.25; %stelling up to CFL=0.4 (reasonable offshoots)
theta = 0.5;

%number of cells in each direction
N = 2^8;

%step sizes
dx = Lx/N;
dt = CFL * dx / u_max;
t_incr = 10;
stride = ceil(t_incr/dt);

nu = 0; %dx^2; %numerical viscosity

%axis
eta_axis = [-2*dx Lx+2*dx -0.5*g_h 1.5*g_h];
% u_axis   = [-2*dx Lx+2*dx -0.5*g_q/g_h 1.5*g_q/g_h];
u_axis   = [-2*dx Lx+2*dx -0.5*g_q/g_h 4.5*g_q/g_h];
q_axis   = [-2*dx Lx+2*dx -0.5*g_q 1.5*g_q];

%cell centers
xc = ((1:N)' - 0.5) * dx;

%edges 
xe = (0:N)' * dx;

%initialize bathymetry
[B,dB_dx] = bathymetry(xc);
B = [B(1); B; B(end)];

%initialize variables
eta = eta0 * ones(N,1);
u = u0     * ones(N+1,1);

%create ghost cell arrays
eta_new = [eta0; eta; g_h];
H_new = eta_new + B;
u_new = [g_q/eta0; u; 0.0];

%intiialize variables for previous time steps
eta_old = eta_new;
H_old = H_new;
u_old = u_new;

u_old_old = u_old;

Hu = calculate_edge_elevations(u_old,H_old);

%total flux
Q_h = trapz(H_old(2:end-1))*dx;
Q_q = trapz(Hu.*u)*dx;

max_CFL = 0;

t_idx = 0;
t = 0;
while t <= t_f
    
    %set ghost cells
    eta_old = set_eta_gc(eta_old,B,g_h);
    u_old = set_u_gc(u_old);
        
    %plot
    plot_state(xc,xe,t,t_idx,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)
   
    max_CFL = max(max(abs(u_old))*dt/dx,max_CFL);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Construct u at non-grid points for semi-Lagrangian approach
    u_star = construct_u_star(u_old,u_old_old,dx,dt);
    
    %gradient of eta at t_n
    d_eta_n_dx = (eta_old(2:end) - eta_old(1:end-1))/dx;    
    
    %Laplacian of u,v at t_n
    Lap_u = Lap_1d(u_star,dx);

    %Casulli's "G" operator
    Gu = construct_Gu(u_star,Lap_u,d_eta_n_dx,g,nu,dt,theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_old = eta_old + B;
    Hu = calculate_edge_elevations(u_star,H_old);
    
    %wet vector
    w = find(Hu > wd_tol);
    
    %friction
    gamma = zeros(N+1,1);
    
    %prefactor used for T and velocity update
    PHI = zeros(N+1,1);
    PHI(w) = Hu(w) ./ (Hu(w) + dt*gamma(w));
    
    %construct V old (linear)
    V_old = H_old * dx;
    V_old(2:end-1) = total_vol(xc,eta_old(2:end-1),dx,wd_tol);

    %construct fluxes
    g_q2 = g*theta*dt*Hu(end)*PHI(end)*eta_old(end);
    [flux_ex,flux_im] = construct_lr_fluxes(Hu,Hu,u_old,Gu,PHI);
    [flux_ex,flux_im] = set_flux_bc_alt(flux_ex, flux_im,g_q,g_q2);
    [flux_diff_ex] = flux_difference(flux_ex);
    [flux_diff_im] = flux_difference(flux_im); %+ theta*dt*g_q2;
        
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
    [eta_update,r,iter] = nonlinear_solve(xc,dx,eta_old,T,b,n_tol,wd_tol,max_iter);
    eta_new(2:end-1) = eta_update;
    
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
    
    %total flux
    [Q_h, Q_q] = update_total_flux(Hu, u_new, B, Q_h, Q_q, g, dt);
end

%flux
Q_h_num = dx*trapz(H_new(2:end-1));
Q_q_num = dx*trapz(Hu.*u_new(2:end-1));
mass_change = [abs(Q_h_num - Q_h) abs(Q_q_num - Q_q)];

%read it Clawpack data
fileID = fopen('pyclaw_soln.txt','r');
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
plot(xd,data(:,1),'LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('H','FontSize',fs)
legend('numerical','bathymetry','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)
axis(eta_axis);


figure(3)
lw = 3;
fs = 24;
plot(xe,Hu.*u_new(2:end-1),'LineWidth',lw+1)
hold on
plot(xe,g_q*ones(N+1),'--','LineWidth',lw)
plot(xd,data(:,2),'LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('Hu','FontSize',fs)
legend('numerical','actual','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)
