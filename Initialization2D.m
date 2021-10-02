function [Lattice,fIn,fEq,fOut] = Initialization2D(Lattice,BCs,ICs,LS2D,Char)
%% Initialization for Lattice Boltzmann Method

% initial conditions of velocity 
Lattice.UxInit     = reshape(ICs.Ux(Lattice.x,Lattice.y,Lattice.UMax),1,Lattice.Nx,Lattice.Ny);
Lattice.UyInit     = reshape(ICs.Uy(Lattice.x,Lattice.y,Lattice.UMax),1,Lattice.Nx,Lattice.Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions of pressure, given as analytical pressure
% or can be obtained by sovling Poisson problem or by MRT model
Lattice.dPdx   = ICs.dPdx*(Lattice.UMax)^2;
Lattice.Pres   = -(Lattice.Lx - Lattice.x)*Lattice.dPdx;

% initial conditions of density 
Lattice.RhoInit = reshape(3*Lattice.Pres + Lattice.Rho,1,Lattice.Nx,Lattice.Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIAL DISTRIBUTION AS EQUILIBRIUM DISTRIBUTION
Lattice.Ux   = Lattice.UxInit;
Lattice.Uy   = Lattice.UyInit;
Lattice.Rho  = Lattice.RhoInit;

for i=1:length(LS2D.cx)
    cu = 3*(LS2D.cx(i)*Lattice.Ux+LS2D.cy(i)*Lattice.Uy);
    fEq(i,:,:) = Lattice.Rho*LS2D.w(i) .* ...
        ( 1  + cu + 1/2*(cu.*cu) - 3/2*(Lattice.Ux.^2+Lattice.Uy.^2) );
end
fIn  = fEq;
fOut = fEq;


%% 5. Main function
cycle = 1;
PlotStep = 10;
StopTime = 1;
ErrorRho = 1;
while cycle < StopTime && ErrorRho > 1e-8
    Lattice.RhoInit = Lattice.Rho;
   %% 6. Boundary Conditions
    % Location = {'Boundary', 'Left', 'Right', 'Up', 'Down'}
    Lattice = BoundaryConditions2D(Lattice,fIn,BCs,'Boundary');
    
    %% 7. Boundary Dynamics 
    % BoundaryDynamics2D= { 'ZouHe2D', 'Regularized2D','FiniteDifference2D'}
    % Location = {'Boundary', 'Left', 'Right', 'Up', 'Down'}
    fIn = BoundaryDynamics2D('Regularized2D',Lattice,BCs,fIn,fEq,LS2D,'Boundary');
    
    %% 8. Collision Step
    [fOut,fEq] = Collision2D(fIn,fEq,fOut,Lattice,LS2D);
    
    %% 9. Streaming Step
    [fIn,Lattice] = Streaming2D(fIn,fOut,Lattice,LS2D);
    % initial conditions of velocity 
    Lattice.Ux     = Lattice.UxInit;
    Lattice.Uy     = Lattice.UyInit;
    %% 10. Visualization
    % available variables, {pressure,velocity}
    if (mod(cycle,PlotStep)==1)
        Visualization2D(Lattice,Char);   
            ErrorRho = sqrt(sum(sum((Lattice.Rho - Lattice.RhoInit).^2))/(Lattice.Nx*Lattice.Ny))
    end    
    %%

    cycle = cycle + 1;
end 
