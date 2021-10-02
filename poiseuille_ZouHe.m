% function [error_r,error_L2,time] = poiseuille_ZouHe( NumCell )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is developed to test the Lattice
% Boltzmann method by the benchmarks of Poiseuille steady plan channel flow
% with analytical solution of u_x(y) = 4*uMax*(y-y^2). It's worth to
% mention that the code is flexible to handle problems from dimensional
% system to dimensionless system and dicrete system, and all the variables
% and constants are provided with both physical expression and rescaled
% expression. The boundary condition used in this code is Zou/He BC
% suggested by Qisu Zou and Xiaoyi He in 1996, with bounceback rule applied
% for non-equilibrium part of density distribution. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:  Qisu Zou and Xiaoyi He, On pressure and velocity flow
% boundary conditions and bounceback for the Lattice Boltzamann BGK model.
% http://www.lbmethod.org/numerics:codes?s[]=cylinder     for cylinder.m

% Copyright (C) 2010 at CMCS EPFL 
% Dr. Matteo Astorinoat E-mail: matteo.astorino@epfl.ch
% Chen Peng E-mail: peng.chen@epfl.ch
% Advisor : Prof. Alfio Quarteroni E-mail: alfio.quarteroni@epfl.ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
 close all
 clc

time = cputime;
%%
% =====================IN DIMENSIONAL SYSTEM==============================
% all the variables and constants are in the physical system with units
% namely we have 

% 1. rho    macroscopic density 
% 2. m      molecular mass
% 3. k_B    Boltzmann constant
% 4. R      specific gas constant
% 5. T      absolute temperture, also written as theta
% 6. x      spatial variable, vector 
% 7. t      temporal variable 
% 8. u      macroscopic velocity, vector 
% 9. c      molecular velocity, vector
% 10.l_0    characteristic length
% 11.t_0    characteristic time
% 12.T_0    characteristic temperture, also theta_0
% 13.u_0    characteristic velocity u_0 = l_0/t_0
% 14.rho_0  characteristic density
% 15.c_s    sound velocity c_s=sqrt(R*T_0)
% 16.nu     shear viscosity
% 17.Re     Reynolds number Re = l_0*l_0/(t_0*nu)
% 18.f      density distribution function depending on rho,m,k_B,R,T,c,u
% 19.tau    relaxation time, tau = nu/c_s^2 + 1./2
% characteristic variables are chosen according to the physical system
%  flow domain looks as follows 
%  |-----------------|
%  |                 |
%  |-----------------|
% the physical system is CGS(cm-s) system
rho_p   = 1;       % density, unit g/cm^3
rho_p_out = 1;     % outflow density
mu_p    = 0.1;     % shear viscocity, unit
nu_p    = mu_p./rho_p;   % shear viscocity, unit poise or g/(cm*s)
lx_p    = 2;       % length in x-direction, unit cm
ly_p    = 1;       % length in y-direction, unit cm
uMax_p  = 1;       % maximal velocity, unit cm/s
maxT_p  = 10;      % final time unit s
rho_p0  = 1;       % characteristic density, unit g/cm^3
l_p0    = 1;       % characteristic length, unit cm
u_p0    = uMax_p;  % characteristic velocity, unit cm^2/s
t_p0    = l_p0./u_p0;% characteristic time, unit s
Re      = u_p0*l_p0./nu_p; % Reynolds number 



%%
% ====================IN DIMENSIONLESS SYSTEM=============================
% all the variables and constants are rescaled by characteristic variables
% so that all the characteristic variables = 1 in the dimensionless system
% especially the following variables

% 1. t_d0    = 1;         % characteristic time,          fixed to be 1
% 2. l_d0    = 1;         % characteristic length,        fixed to be 1
% 3. u_d0    = 1;         % characteristic velocity,      fixed to be 1    
% 4. rho_d    = 1;         % density, can be changed for compressible flows
% 5. T_d      = 1;         % temperture, can be changed for thermal dynamics
% 6. lx_d     = 1;         % length in x-direction, can be changed
% 7. ly_d     = 1;         % length in y-direction, can be changed
% 8. lz_d     = 1;         % length in z-direction, can be changed, for 3D
% 9. ux_d     = 1;         % velocity in x-direction, can be changed
% 10.uy_d     = 1;         % velocity in y-direction, can be changed
% 11.uz_d     = 1;         % velocity in z-direction, can be changed, for 3D
% 12.nu     = 1./Re;     % shear viscosity, Reynolds number fixed

rho_d   = rho_p/rho_p0;       % density
rho_d_out = rho_p_out/rho_p0; % outflow density 
nu_d    = 1./Re;   % shear viscocity
lx_d    = lx_p/l_p0;       % length in x-direction
ly_d    = ly_p/l_p0;       % length in y-direction
uMax_d  = uMax_p/u_p0;     % maximal velocity
maxT_d  = maxT_p/t_p0;     % final time
l_d0    = 1;         % characteristic length,        fixed to be 1
u_d0    = 1;         % characteristic velocity,      fixed to be 1 
t_d0    = 1;         % characteristic time,          fixed to be 1
dpdx_d  = -8/Re;

%%
% ==================IN DISCRETE SYSTEM====================================
N       = 32;%NumCell;           % number of cells
dx = l_d0./N;           % grid resolution
% uMax_lb= 0.01;        % fix uMax_lb
% dt = uMax_lb*dx/uMax_d; % time linearly depend on spatial resolution
dt  = dx^2;             % time resolution: conserve the shear viscocity
nu_lb = nu_d*dt/dx^2;   % lattice viscocity           
tau = 3*nu_lb+0.5;      % relaxation time
uMax_lb= dt/dx*uMax_d;  % lattice maximal velocity < c_s sound velocity
% rho_lb = rho_d;         % lattice density
rho_lb_out = rho_d_out; % lattice outflow density


% dimensionless system for the grid
gridx_d = 0:dx:lx_d;                    % grid in x-direction
gridy_d = 0:dx:ly_d;                    % grid in y-direction
[y_d,x_d] = meshgrid(gridy_d,gridx_d);  % get coordinate of matrix indices
% % lb system for the grid
% gridx_lb = gridx_d/dx;                    % grid in x-direction
% gridy_lb = gridy_d/dx;                    % grid in y-direction
% x_lb = x_d/dx;
% y_lb = y_d/dx;
dpdx_lb =dpdx_d*(dt/dx)^2;
rho_lb = rho_lb_out+3*(lx_d-x_d)*(-dpdx_lb);



[Nx,Ny] = size(x_d);     % get size of matrix
in = 1;  out = Nx;       % in and out boundary x-direction indices
col = 2:Ny-1;            % in and out boundary y-direction indices
up = Ny; down = 1;       % up and down wall boundary y-direction indices
row = 2:Nx-1;            % up and down wall boundary x-direction indices
maxT_lb   = maxT_d/dt;   % terminal time steps


% D2Q9 LATTICE CONSTANTS
%  7  3  6
%   \ | /
% 4 - 1 - 2 
%   / | \
%  8  5  9
w   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];% quadrature weights
cx  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];% rescaled by c_l = dx/dt
cy  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];% rescaled by c_l = dx/dt
opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7]; % opposite indices


% INITIAL CONDITION: Poiseuille profile at equilibrium
ux_a = 4*uMax_lb*(y_d-y_d.*y_d); % analytical velocity in x-direction
uy_a = zeros(Nx,Ny);      % analytical velocity in y-direction
ux = ux_a;                % initial velocity in x-direction
uy = uy_a;                % initial velocity in y-direction
% ux = zeros(Nx,Ny);      % alternative for initial velocity = 0
errorPost = 1;            % posterior error



% INITIAL DISTRIBUTION AS EQUILIBRIUM DISTRIBUTION
for i=1:9
    cu = 3*(cx(i)*ux+cy(i)*uy);
    fIn(i,:,:) = rho_lb*w(i) .* ...
                   ( 1  + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2) );
end

% TEMPORAL MACROSCOPIC VARIABLES 
rho_t = sum(fIn);                                   
ux_t  = reshape ((cx * reshape(fIn,9,Nx*Ny)),1,Nx,Ny) ./rho_t;
uy_t  = reshape ((cy * reshape(fIn,9,Nx*Ny)),1,Nx,Ny) ./rho_t;

%%
% MAIN LOOP (TIME CYCLES)
tPlot  = 5000;             % plotting time steps 
tic
cycle = 2;                % time cycle
while cycle < maxT_lb && errorPost > 1e-6 
    time = cputime;
    % MACROSCOPIC VARIABLES
    rho = rho_t;
    ux  = ux_t;
    uy  = uy_t;
    
    % MACROSCOPIC BOUNDARY CONDITIONS
      % Inlet: Poiseuille profile
    ux(:,in,:)  = 4*uMax_lb*(gridy_d - gridy_d.*gridy_d);
    uy(:,in,:) = 0;
    rho(:,in,col) = 1 ./ (1-ux(:,in,col)) .* ( ...
        sum(fIn([1,3,5],in,col)) + 2*sum(fIn([4,7,8],in,col)) );
    
%       % Outlet: Poiseuille profile
%     ux(:,out,:)  = 4*uMax_lb*(gridy_d - gridy_d.*gridy_d);
%     uy(:,out,:) = 0;
%     rho(:,out,col) = 1 ./ (1+ux(:,out,col)) .* ( ...
%         sum(fIn([1,3,5],out,col)) + 2*sum(fIn([2,6,9],out,col)) );    
    
      % Outlet: Constant pressure
    rho(:,out,:) = 1;
    ux(:,out,col) = -1 + 1 ./ (rho(:,out,col)) .* ( ...
        sum(fIn([1,3,5],out,col)) + 2*sum(fIn([2,6,9],out,col)) );
    uy(:,out,col)  = 0;
      
      % UpWall: Dirichlet boundary condition
    ux(:,:,up)  = 0;
    uy(:,:,up)  = 0;
    rho(:,row,up) = 1./(1 + uy(:,row,up)).*( ...
        sum(fIn([1,2,4],row,up)) + 2*sum(fIn([3,6,7],row,up)) );
%     rho(:,in,up) = 1/(1+1+sqrt(2))*(rho(:,in,up-1)+...
%         rho(:,in+1,up)+sqrt(2)*rho(:,in+1,up-1)); % 1st order
    rho(:,in,up) = rho(:,in,up-1);
    
      % DownWall: Dirichlet boundary condition
    ux(:,:,down) = 0;
    uy(:,:,down) = 0;    
    rho(:,row,down) = 1./(1 - uy(:,row,down)).*( ...
        sum(fIn([1,2,4],row,down)) + 2*sum(fIn([5,8,9],row,down)) );
%     rho(:,in,down) = 1/(1+1+sqrt(2))*(rho(:,in,down+1)+...
%                 rho(:,in+1,down)+sqrt(2)*rho(:,in+1,down+1)); % 1st order
    rho(:,in,down) = rho(:,in,down+1);
    %%
    % MICROSCOPIC BOUNDARY CONDITIONS
      % INLET (Zou/He BC)
    fIn(2,in,col) = fIn(4,in,col) + 2/3*rho(:,in,col).*ux(:,in,col); 
    fIn(6,in,col) = fIn(8,in,col) + 1/2*(fIn(5,in,col)-fIn(3,in,col)) ... 
                                    + 1/2*rho(:,in,col).*uy(:,in,col) ...
                                    + 1/6*rho(:,in,col).*ux(:,in,col); 
    fIn(9,in,col) = fIn(7,in,col) + 1/2*(fIn(3,in,col)-fIn(5,in,col)) ... 
                                    - 1/2*rho(:,in,col).*uy(:,in,col) ...
                                    + 1/6*rho(:,in,col).*ux(:,in,col); 
      % OUTLET (Zou/He BC)
    fIn(4,out,col) = fIn(2,out,col) - 2/3*rho(:,out,col).*ux(:,out,col); 
    fIn(8,out,col) = fIn(6,out,col) + 1/2*(fIn(3,out,col)-fIn(5,out,col)) ... 
                                      - 1/2*rho(:,out,col).*uy(:,out,col) ...
                                      - 1/6*rho(:,out,col).*ux(:,out,col); 
    fIn(7,out,col) = fIn(9,out,col) + 1/2*(fIn(5,out,col)-fIn(3,out,col)) ... 
                                      + 1/2*rho(:,out,col).*uy(:,out,col) ...
                                      - 1/6*rho(:,out,col).*ux(:,out,col); 
                                  
      % UP WALL (Zou/He BC)                             
    fIn(5,row,up) = fIn(3,row,up) - 2/3*rho(:,row,up).*uy(:,row,up); 
    fIn(8,row,up) = fIn(6,row,up) + 1/2*(fIn(2,row,up)-fIn(4,row,up)) ... 
                                    - 1/2*rho(:,row,up).*ux(:,row,up) ...
                                    - 1/6*rho(:,row,up).*uy(:,row,up); 
    fIn(9,row,up) = fIn(7,row,up) - 1/2*(fIn(2,row,up)-fIn(4,row,up)) ... 
                                    + 1/2*rho(:,row,up).*ux(:,row,up) ...
                                    - 1/6*rho(:,row,up).*uy(:,row,up);
      % DOWN WALL (Zou/He BC)     
    fIn(3,row,down) = fIn(5,row,down) + 2/3*rho(:,row,down).*uy(:,row,down); 
    fIn(6,row,down) = fIn(8,row,down) - 1/2*(fIn(2,row,down)-fIn(4,row,down)) ... 
                                    + 1/2*rho(:,row,down).*ux(:,row,down) ...
                                    + 1/6*rho(:,row,down).*uy(:,row,down);  
    fIn(7,row,down) = fIn(9,row,down) + 1/2*(fIn(2,row,down)-fIn(4,row,down)) ... 
                                    - 1/2*rho(:,row,down).*ux(:,row,down) ...
                                    + 1/6*rho(:,row,down).*uy(:,row,down);
      % INLET CORNERS (Zou/He BC): consistent with Poiseuille profile
    % since the density is unknown, we extrapolate density in 2nd order
    % rho(:,in,up) = 2*rho(:,in,up-1) - rho(:,in, up-2); % 2nd order
%     rho(:,in,up) = rho(:,in,up-1); % 1st order
%     rho(:,in,up) = 1/2*(rho(:,in,up-1)+rho(:,in+1,up)); % 1st order
    fIn(:,in,up) = fIn(opp,in,up);
    fIn([7,9],in,up) = 1./2*(rho(:,in,up) - ...
                                sum(fIn([1,2,3,4,5,6,8],in,up)));
    %rho(:,in,down) = 2*rho(:,in,down+1) - rho(:,in,down+2); %2nd order
%     rho(:,in,down) = rho(:,in,down+1); % 1st order
%     rho(:,in,down) = 1/2*(rho(:,in,down+1)+rho(:,in+1,down)); % 1st order
    fIn(:,in,down) = fIn(opp,in,down);
    fIn([6,8],in,down) = 1./2*(rho(:,in,down) - ...
                                sum(fIn([1,2,3,4,5,7,9],in,down)));
      
      % OUTLET CORNERS (Zou/He BC): consistent with Constant pressure
    fIn(:,out,up) = fIn(opp,out,up);
    fIn([6,8],out,up) = 1./2*(rho(:,out,up) - ...
                                sum(fIn([1,2,3,4,5,7,9],out,up)));
    fIn(:,out,down) = fIn(opp,out,down);
    fIn([7,9],out,down) = 1./2*(rho(:,out,down) - ...
                                sum(fIn([1,2,3,4,5,6,8],out,down)));

%%

    % COLLISION STEP
    for i=1:9
       cu = 3*(cx(i)*ux+cy(i)*uy);
       fEq(i,:,:)  = rho .* w(i) .* ...
                       ( 1 + cu + 1/2*(cu.*cu)  - 3/2*(ux.^2+uy.^2) );
       fOut(i,:,:) = fIn(i,:,:) - 1./tau .* (fIn(i,:,:)-fEq(i,:,:));
    end

%%
   
    % STREAMING STEP
      % INNER PART
    for i = 1:9
        fIn(i,row+cx(i),col+cy(i)) = fOut(i,row,col); 
    end

      % INLET and OUTLET
    for i = [1,2,3,5,6,9]
        fIn(i,in+cx(i),col+cy(i)) = fOut(i,in,col);
        fIn(opp(i),out+cx(opp(i)),col+cy(opp(i))) = fOut(opp(i),out,col);
    end
      % UP WALL and DOWN WALL
    for i = [1,2,4,5,8,9]
        fIn(i,row+cx(i),up+cy(i)) = fOut(i,row,up);
        fIn(opp(i),row+cx(opp(i)),down+cy(opp(i))) = fOut(opp(i),row,down);
    end
    
      % INLET-UP and OUTLET-DOWN CORNER
    for i = [1,2,5,9]
        fIn(i,in+cx(i),up+cy(i)) = fOut(i,in,up);
        fIn(opp(i),out+cx(opp(i)),down+cy(opp(i))) = fOut(opp(i),out,down);
    end
      % INLET-DOWN and OUTLET-UP CORNER
    for i = [1,2,3,6]
        fIn(i,in+cx(i),down+cy(i)) = fOut(i,in,down);
        fIn(opp(i),out+cx(opp(i)),up+cy(opp(i))) = fOut(opp(i),out,up);
    end
    %%
        
    % TEMPORAL MACROSCOPIC VARIABLES 
    rho_t = sum(fIn);                                   
    ux_t  = reshape ((cx * reshape(fIn,9,Nx*Ny)),1,Nx,Ny) ./rho_t;
    uy_t  = reshape ((cy * reshape(fIn,9,Nx*Ny)),1,Nx,Ny) ./rho_t;
    
    %posterior error between velocity in the current step and in the next.
    errorPost = sqrt(sum(sum((ux_t - ux).^2 + (uy_t - uy).^2))./...
                                        sum(sum(ux.^2 + uy.^2)) );
    
    rho = rho_t;
    ux  = ux_t;
    uy  = uy_t;                                    
                                    
    % VISUALIZATION
    if (mod(cycle,tPlot)==1)
%         % velocity
%         ux_physical=u_p0*(dx/dt)*ux;
%         uy_physical=u_p0*(dx/dt)*uy;
%         
%         u = reshape( sqrt(ux_physical.^2+uy_physical.^2),Nx,Ny);
%         imagesc(u');colorbar
        % density
        rho = reshape(rho, Nx, Ny);
        imagesc(rho');colorbar
        axis equal off; drawnow
        toc
        pause
    end

    time_physical = cycle*dt*t_p0;
    cycle = cycle + 1;
end

    % VISUALIZATION
        % velocity
        ux_physical=u_p0*(dx/dt)*ux;
        uy_physical=u_p0*(dx/dt)*uy;
        u = reshape( sqrt(ux_physical.^2+uy_physical.^2),Nx,Ny);
        figure
        imagesc(u');colorbar
        axis equal off; drawnow
        
        % density
        rho = reshape(rho, Nx, Ny);
        figure
        imagesc(rho');colorbar
        axis equal off; drawnow

ux = reshape(ux_t,Nx,Ny);
uy = reshape(uy_t,Nx,Ny);

% the relative error for LB velocity is defined as 
error_r = sqrt(sum(sum( (ux_a - ux).^2 + (uy_a - uy).^2 ))./...
                                        sum(sum(ux.^2 + uy.^2)) );
disp('The relative error between analytical velocity and LB velocity is : ')   
disp(error_r)

% the L2 error for LB velocity is defined as 
error_L2 = sqrt(sum(sum((ux_a - ux).^2 + (uy_a - uy).^2))/(dt/dx)^2/(Nx*Ny));
disp('The L2 error between analytical velocity and LB velocity is :')
disp(error_L2)

disp('cputime is : ')
time = cputime-time          % cputime spent on the running the programme