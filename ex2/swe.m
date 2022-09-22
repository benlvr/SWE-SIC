%Based on example 4.1.1 from SWASHES by Delestre et al (2016)
%This is a flat bottom with no diffusion, friction


clear all
close all

addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include')
addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include/limiters')

%parameters
t_f = 6;
Lx = 10; %m
g = 9.81; %m/s^2

%boundary and initial data
% u0 = 0.0;
x0 = Lx/2;
h_l = 0.005;
h_r = 0.00;
u_l = 0.0;
u_r = 0.0;
g_1 = h_l;
g_2 = h_r;
bdry = 2; %boundary where we apply an elevation bc

%Newton tolerance
n_tol  = 1.0e-10;
wd_tol = 1.0e-10;
max_iter = 10;

%timestepping parameters
u_max = sqrt(g*max(h_l,h_r));
CFL = .025;
theta = 0.5;

%number of cells in each direction
N = 2^8;

%step sizes
dx = Lx/N;
dt = CFL * dx / u_max;
t_incr = t_f/100;
stride = ceil(t_incr/dt);

%numerical viscosity
nu = 0.0*dx^2;

%axes
eta_axis = [0 Lx -min(h_r,h_l) 2*max(h_r,h_l)];
u_axis   = [0 Lx -0.25*u_max 2.5*u_max];
q_axis   = [0 Lx -0.25*u_max u_max];

%cell centers
xc = ((1:N)' - 0.5) * dx;

%edges 
xe = (0:N)' * dx;

%initialize bathymetry
B = bathymetry(xc);
B = [B(1); B; B(end)];

%initialize variables
eta = h_l * (xc <= x0) + h_r*(xc >= x0);
H = eta;
u = u_l * (xe <= x0) + u_r*(xe >= x0);

%create ghost cell arrays
eta_new = [h_l-B(1); eta; h_r-B(end)];
H_new = eta_new + B;
u_new = [u_l; u; u_r];

%intiialize variables for previous time steps
eta_old = eta_new;
H_old = H_new;
u_old = u_new;

u_old_old = u_old;

gamma = zeros(N+1,1);
Hu = calculate_edge_elevations(u_old,H_old);
w = find(Hu > wd_tol); 

%total flux
Q_h = sum(H)*dx;
Q_q = trapz(Hu.*u)*dx - g*trapz(Hu.*(eta_new(2:end) - eta_new(1:end-1)));

t_idx = 0;
t = 0;
while t <= t_f
    
    %set ghost cells
    eta_old = set_eta_gc(eta_old,dx,t,g_1,g_2);
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
            
    %wet vectors
    w = find(Hu > wd_tol); 
          
    %friction
    gamma = zeros(N+1,1);
    
    %prefactor used for T and velocity update
    PHI = zeros(N+1,1);
    PHI(w) = Hu(w) ./ (Hu(w) + dt*gamma(w));
    
    %construct V old (linear)
    V_old = H_old * dx;
    V_old(2:end-1) = total_vol(xc,eta_old(2:end-1),dx,wd_tol);
    
    %alternative
    [flux_ex,flux_im] = construct_lr_fluxes(Hu,Hu,u_old,Gu,PHI);
    [flux_ex,flux_im] = set_flux_bc_alt(flux_ex,flux_im,h_l*u_l);
    [flux_diff_ex] = flux_difference(flux_ex);
    [flux_diff_im] = flux_difference(flux_im);
        
    %update volume
    V_new = V_old(2:end-1) - dt*flux_diff_ex;
    
    %reconstruct total elevation H
    H_new(2:end-1) = V_new / dx;
    
    %calculate surface elevation eta
    eta_new(2:end-1) = H_new(2:end-1) - B(2:end-1);
    
    T = construct_T(Hu,eta_old,theta,dt,dx,g,PHI);
    b = construct_b(V_old,flux_diff_ex,flux_diff_im,dt,theta);
    [T,b] = add_elevation_bc(Hu,eta_old,PHI,T,b,g,dt,theta,bdry);
    [eta_update,r,iter] = nonlinear_solve(xc,dx,eta_old,T,b,n_tol,wd_tol,max_iter);    
    eta_new(2:end-1) = eta_update;
    
    eta_new = set_eta_gc(eta_new,dx,t+dt,g_1,g_2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_eta_n1_dx = (eta_new(2:end) - eta_new(1:end-1))/dx;
    
    %update u
    RHS = Gu - dt*g*theta*d_eta_n1_dx;
    u_new(2:end-1) = PHI.*RHS;
    
    %set u BC
    u_new = set_u_bc(u_new,dx,u_l);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update time
    t_idx = t_idx + 1;
    t     = t + dt;
        
    %update previous values to new values
    u_old_old = u_old;
    eta_old = eta_new;
    u_old = u_new;
    
end

%analytic soln
[eta_act,u_act] = analytic_soln(h_l,h_r,x0,xc,xe,t,g,n_tol);

%error in SS soln
L_2_e  = [norm(eta_act - eta_new(2:end-1)) norm(u_act - u_new(2:end-1))] * sqrt(dx)
L_oo_e = [max(abs(eta_act - eta_new(2:end-1))) max(abs(u_act - u_new(2:end-1)))]

%plot final state
plot_state(xc,xe,t,0,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)

%flux
Q_h_num = dx*sum(H_new(2:end-1));
Q_q_num = dx*trapz(Hu.*u_new(2:end-1)) - g*dx*sum(Hu.*d_eta_n_dx);
mass_change = [abs(Q_h_num - Q_h) abs(Q_q_num - Q_q)]

[~,u_act_c] = analytic_soln(h_l,h_r,x0,xc,xc,t,g,n_tol);
[Fr_act] = calculate_froude(eta_act,u_act_c,g);
max_Fr_act = max(Fr_act);

%read it Clawpack data
if h_r > 0
    fileID = fopen('pyclaw_soln.txt','r');
else
    fileID = fopen('pyclaw_dry_soln.txt','r');
end
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
plot(xc,eta_new(2:end-1),'LineWidth',lw+1)
hold on
plot(xc,-B(2:end-1),'--','LineWidth',lw)
plot(xc,eta_act,'--','LineWidth',lw)
plot(xd,data(:,1),'LineWidth',lw)
xlabel('x direction','FontSize',fs)
ylabel('H','FontSize',fs)
% legend('numerical','actual','bathymetry','FontSize',fs)
legend('numerical','bathymetry','actual','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)

figure(3)
lw = 3;
fs = 24;
plot(xe,Hu.*u_new(2:end-1),'LineWidth',lw+1)
hold on
plot(xc,eta_act.*u_act_c,'--','LineWidth',lw)
plot(xd,data(:,2),'LineWidth',lw-1)
xlabel('x direction','FontSize',fs)
ylabel('Hu','FontSize',fs)
legend('numerical','actual','Clawpack','FontSize',fs)
title('t=',num2str(t),'FontSize',fs)
set(gca,'fontsize',fs)
