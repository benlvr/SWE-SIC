
%Based on example 3.1.1 from SWASHES by Delestre et al (2016)
%This is an immersed bump with no diffusion, friction

%VARIABLES
%H:   total elevatation of water, H = eta + B
%eta: water elevation above the plane y = 0
%B:   depth of terrain below the plane y = 0
%u:   depth-averaged velocity in x-direction

%EQUATIONS
%Mass: H_t + (H*u)_x = 0
%Momentum: H*(u_t + uu_x) = -g*H*eta_x - gamma*u

clear all
close all

addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include')

%parameters
t_f = 15;
Lx = 25; %m
g = 9.81; %m/s^2
bdry = 2;

%Newton tolerance
n_tol  = 1.0e-8;
wd_tol = 1.0e-8;
max_iter = 10;

%timestepping parameters
u_max = 1;
CFL = 5;
theta = 0.5;

%number of cells in each direction
N = 2^7;

%step sizes
dx = Lx/N;
dt = CFL * dx / u_max;
t_incr = 0.1;
stride = ceil(t_incr/dt);

%numerical viscosity
nu = 0.0;

eta_axis = [0 Lx -1 3];
u_axis = [0 Lx -1 1];

%cell centers
xc = ((1:N)' - 0.5) * dx;

%edges 
xe = (0:N)' * dx;

%initialize bathymetry
B = bathymetry(xc);
B = [B(1); B; B(end)];

%initialize variables
eta0 = 0.5;
eta = eta0 * ones(N,1);
u = zeros(N+1,1);

%boundary conditions
g_eta1 = eta0;
g_eta2 = eta0;
g_u = 0;

%create ghost cell arrays
eta_new = [g_eta1; eta; g_eta2];
H_new = eta_new + B;
u_new = [0; u; 0];

%intiialize variables for previous time steps
eta_old = eta_new;
H_old = H_new;
u_old = u_new;

gamma = zeros(N+1,1);
Hu = calculate_edge_elevations(u_old,H_old);
u_old_old = u_new;

t_idx = 0;
t = 0;
while t <= t_f
    
    %set ghost cells
    eta_old = set_eta_gc(eta_old,g_eta1,g_eta2);
    u_old = set_u_gc(u_old);
    
    %plot
    if mod(t_idx,stride) == 0
        drawnow     
        plot(xc,eta_new(2:end-1))
        axis(eta_axis)    
        xlabel('x')
        ylabel('\eta')
        title(t)
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Construct u at non-grid points for semi-Lagrangian approach
%     u_star = construct_u_star(u_old,u_old_old,dx,dt);
    u_star = stelling_flux(H_old,Hu,u_old,dx,dt,wd_tol);
    
    %gradient of eta at t_n
    d_eta_n_dx = (eta_old(2:end) - eta_old(1:end-1))/dx;
    
    
    %Laplacian of u,v at t_n
    Lap_u = Lap_1d(u_star,dx);

    %Casulli's "G" operator
    Gu = u_star(2:end-1) + dt*(nu*Lap_u - g*(1-theta)*d_eta_n_dx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_old = eta_old + B;
    Hu = calculate_edge_elevations(u_old,H_old);
    
    %wet and dry vectors
    w = find(Hu > wd_tol);  
    d = find(Hu <= wd_tol);  
    
    %prefactor used for T and velocity update
    PHI = zeros(N+1,1);
    PHI(w) = Hu(w) ./ (Hu(w) + dt*gamma(w));
    
    %construct V old (linear)
    V_old = H_old * dx;
    V_old(2:end-1) = total_vol(xc,eta_old(2:end-1),dx,wd_tol);
    
%     %construct fluxes
    [flux_x_ex,flux_x_im] = construct_fluxes(Hu,u_old,Gu,PHI);
    [flux_x_ex,flux_x_im] = set_flux_bc(flux_x_ex,flux_x_im,Hu,u_old,Gu,PHI,g_eta1*g_u);
        
    %update volume
    V_new = V_old(2:end-1) - dt*flux_x_ex;
    
    %reconstruct total elevation H
    H_new(2:end-1) = V_new / dx;
    
    %calculate surface elevation eta
    eta_new(2:N+1) = H_new(2:N+1) - B(2:N+1);
    
    T = construct_T(Hu,eta_old,theta,dt,dx,g,PHI);
    b = construct_b(V_old,flux_x_ex,flux_x_im,dt,theta);
%     [T,b] = add_elevation_bc(Hu,eta_old,PHI,T,b,g,dt,theta,1);
    [T,b] = add_elevation_bc(Hu,eta_old,PHI,T,b,g,dt,theta,2);
    [eta_update,r,iter] = nonlinear_solve(xc,dx,eta_old,T,b,n_tol,wd_tol,max_iter);
%     eta_new(2:end-1) = eta_update;
    
    eta_new = set_eta_gc(eta_new,g_eta1,g_eta2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_eta_n1_dx = (eta_new(2:end) - eta_new(1:end-1))/dx;
    
    %update u
    RHS = Gu - dt*g*theta*d_eta_n1_dx;
    u_new(2:end-1) = PHI.*RHS;
    
    %set u BC (on boundaries 2 and 4, boundaries 1 and 3 already set by ghost cells)
    u_new = set_u_bc(u_new,g_u);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update time
    t_idx = t_idx + 1;
    t     = t + dt;
        
    %update previous values to new values
    eta_old = eta_new;
    u_old_old = u_old;
    u_old = u_new;
    
end

g_q = 0;
eta_act = eta;
u_act = zeros(N+1,1);
e_eta = abs(eta_act - eta_old(2:end-1));
max(e_eta)

H_act = B(2:end-1) + eta_act;

hold on
plot(xc,-B(2:end-1),'--')
plot(xc,eta_act)

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
plot(xc,H_act,'--','LineWidth',lw)
plot(xd,data(:,1),'LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('H','FontSize',fs)
% legend('numerical','actual','bathymetry','FontSize',fs)
legend('numerical','bathymetry','actual','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)
axis(eta_axis);


figure(3)
lw = 3;
fs = 24;
plot(xe,u_new(2:end-1),'LineWidth',lw+1)
hold on
plot(xe,u_act,'--','LineWidth',lw)
plot(xd,data(:,2),'LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('u','FontSize',fs)
% legend('numerical','actual','bathymetry','FontSize',fs)
legend('numerical','actual','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)
axis(u_axis);
