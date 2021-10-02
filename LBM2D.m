clc
clear all
close all
%% 1. Problem Definition
% Physical system is defined and rescaled to dimensionless system, together
% with boundary conditions and initial conditions as well as the 
% characteristic constants for the problem 

% [Dimensionless, BCs,ICs, Char] = CavityFlow2DProblemDefinition;
% [Dimensionless, BCs,ICs,Char] = PoiseuilleFlow2DProblemDefinition;

[DimensionlessP,DimensionlessQ,DimensionlessC, BCsP, BCsQ,BCsC, ICsP,ICsQ,ICsC,CharP,CharQ,CharC] = Aneurysm2DProblemDefinition;

%% 2. Lattice Structure
LS2D = LatticeStructure2D;

%% 3. Discretization
% construct the spatial lattice for every block with given resolution
% the default temporal resolution is square of spatial one in single block
% while for multiblock we can choose the Order as {1,0} other than square 2
NumCell = 32; % number of cells for edge with characterics length
Order   = 2;  % temporal resolution wrt spatial dt = dx^{Order},{2,1,0}

% block and lattice can be built and extended to multiblock or multigrid model

% Block   = [0,Dimensionless.Lx;0,Dimensionless.Ly];
% [Lattice,BCs,Char] = Discretization2D(NumCell,Order,Block,Dimensionless,BCs,Char);
BlockP   = [0,DimensionlessP.Lx;0,DimensionlessP.Ly];
[LatticeP,BCsP,CharP] = Discretization2D(NumCell,Order,BlockP,DimensionlessP,BCsP,CharP,LS2D);

BlockQ   = [DimensionlessP.Lx,DimensionlessP.Lx+DimensionlessP.Ly;-3*DimensionlessP.Ly,DimensionlessP.Lx-3*DimensionlessP.Ly];
[LatticeQ,BCsQ,CharQ] = Discretization2D(NumCell,Order,BlockQ,DimensionlessQ,BCsQ,CharQ,LS2D);

BlockC   = [DimensionlessP.Lx,DimensionlessP.Lx+2*DimensionlessP.Ly;0,2*DimensionlessP.Ly];
[LatticeC,BCsC,CharC] = Discretization2D(NumCell,Order,BlockC,DimensionlessC,BCsC,CharC,LS2D);


%% 4. Initialization
% [Lattice,fIn,fEq,fOut] = Initialization2D(Lattice,BCs,ICs,LS2D,Char);

[LatticeP,fInP,fEqP,fOutP] = Initialization2D(LatticeP,BCsP,ICsP,LS2D,CharP);
[LatticeQ,fInQ,fEqQ,fOutQ] = Initialization2D(LatticeQ,BCsQ,ICsQ,LS2D,CharQ);
[LatticeC,fInC,fEqC,fOutC] = Initialization2D(LatticeC,BCsC,ICsC,LS2D,CharC);

fInC(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny) = fInP(:,LatticeP.Nx,:);
fInC(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny) = fInQ(:,:,1);

LatticeC.Ux(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny) = LatticeP.Ux(:,LatticeP.Nx,:);
LatticeC.Uy(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny) = LatticeP.Uy(:,LatticeP.Nx,:);
LatticeC.Rho(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny)= LatticeP.Rho(:,LatticeP.Nx,:);  

LatticeC.Ux(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny) = LatticeQ.Ux(:,:,1);
LatticeC.Uy(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny) = LatticeQ.Uy(:,:,1);
LatticeC.Rho(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny)= LatticeQ.Rho(:,:,1);     

%% 5. Main function
PlotStep  = 10;
cycle = 1;
while cycle < LatticeP.TimeMax
    %% 6. Boundary Conditions
    % Location = {'Boundary', 'Left', 'Right', 'Up', 'Down'}
%     Lattice  = BoundaryConditions2D(Lattice,fIn,BCs,'Boundary');

    LatticeC = BoundaryConditions2D(LatticeC,fInC,BCsC,'Curve'); 
    
    LatticeC.Ux(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny) = LatticeP.Ux(:,LatticeP.Nx,:);
    LatticeC.Uy(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny) = LatticeP.Uy(:,LatticeP.Nx,:);
    LatticeC.Rho(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny)= LatticeP.Rho(:,LatticeP.Nx,:);  

    LatticeC.Ux(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny) = LatticeQ.Ux(:,:,1);
    LatticeC.Uy(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny) = LatticeQ.Uy(:,:,1);
    LatticeC.Rho(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny)= LatticeQ.Rho(:,:,1);  

    LatticeP = BoundaryConditions2D(LatticeP,fInP,BCsP,'Boundary');
    
    LatticeP.Ux(:,LatticeP.Nx,:) = LatticeC.Ux(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny);
    LatticeP.Uy(:,LatticeP.Nx,:) = LatticeC.Uy(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny);
    LatticeP.Rho(:,LatticeP.Nx,:)= LatticeC.Rho(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny);  
    
    
    LatticeQ = BoundaryConditions2D(LatticeQ,fInQ,BCsQ,'Boundary');
    
    LatticeQ.Ux(:,:,1) = LatticeC.Ux(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny);
    LatticeQ.Uy(:,:,1) = LatticeC.Uy(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny);
    LatticeQ.Rho(:,:,1)= LatticeC.Rho(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny);      
    

    
    %% 7. Boundary Dynamics 
    % BoundaryDynamics2D= { 'ZouHe2D', 'Regularized2D','FiniteDifference2D'}
    % Location = {'Boundary', 'Left', 'Right', 'Up', 'Down'}
%     fIn = BoundaryDynamics2D('Regularized2D',Lattice,BCs,fIn,fEq,LS2D,'Boundary');

%     fInC = BoundaryDynamics2D('ZouHe2D',LatticeC,BCsC,fInC,fEqC,LS2D,'Curve');

    fInC(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny) = fInP(:,LatticeP.Nx,:);  
    fInC(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny) = fInQ(:,:,1);  

    fInP = BoundaryDynamics2D('Regularized2D',LatticeP,BCsP,fInP,fEqP,LS2D,'Boundary');
    fInP(:,LatticeP.Nx,:) = fInC(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny);   
    
    fInQ = BoundaryDynamics2D('ZouHe2D',LatticeQ,BCsQ,fInQ,fEqQ,LS2D,'Boundary'); 
    fInQ(:,:,1) = fInC(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny);    
    

%     fIn = BoundaryDynamics2D('HalfWay2D',Lattice,BCs,fIn,fEq,LS2D,'Left');
%     fIn = BoundaryDynamics2D('HalfWay2D',Lattice,BCs,fIn,fEq,LS2D,'Right');  
%     fIn = BoundaryDynamics2D('ZouHe2D',Lattice,BCs,fIn,fEq,LS2D,'Up');
%     fIn = BoundaryDynamics2D('HalfWay2D',Lattice,BCs,fIn,fEq,LS2D,'Down');      
    
%     fIn = BoundaryDynamics2D('FiniteDifference2D',Lattice,BCs,fIn,fEq,LS2D,'Left');
%     fIn = BoundaryDynamics2D('FiniteDifference2D',Lattice,BCs,fIn,fEq,LS2D,'Right');  
%     fIn = BoundaryDynamics2D('FiniteDifference2D',Lattice,BCs,fIn,fEq,LS2D,'Up');
%     fIn = BoundaryDynamics2D('FiniteDifference2D',Lattice,BCs,fIn,fEq,LS2D,'Down');     
    %% 8. Collision Step
%     [fOut,fEq] = Collision2D(fIn,fEq,fOut,Lattice,LS2D);

    [fOutP,fEqP] = Collision2D(fInP,fEqP,fOutP,LatticeP,LS2D);
    [fOutQ,fEqQ] = Collision2D(fInQ,fEqQ,fOutQ,LatticeQ,LS2D);
    [fOutC,fEqC] = Collision2D(fInC,fEqC,fOutC,LatticeC,LS2D);
    %% 9. Streaming Step

    fOutC(:,1,(LatticeC.Ny+1)/2:LatticeC.Ny) = fOutP(:,LatticeP.Nx,:);
    fOutC(:,1:(LatticeC.Nx+1)/2,LatticeC.Ny) = fOutQ(:,:,1);
    [fInC,LatticeC] = Streaming2D(fInC,fOutC,LatticeC,LS2D,BCsC,'C');        
    
    [fInP,LatticeP] = Streaming2D(fInP,fOutP,LatticeP,LS2D,BCsP,'P');
    
    [fInQ,LatticeQ] = Streaming2D(fInQ,fOutQ,LatticeQ,LS2D,BCsQ,'Q');


    fInC([2,6,9],1,(LatticeC.Ny+1)/2:LatticeC.Ny) = fInP([2,6,9],LatticeP.Nx,:);  
    fInP(LS2D.opp([2,6,9]),LatticeP.Nx,:) = fInC(LS2D.opp([2,6,9]),1,(LatticeC.Ny+1)/2:LatticeC.Ny);
    fInC([5,8,9],1:(LatticeC.Nx+1)/2,LatticeC.Ny) = fInQ([5,8,9],:,1);    
    fInQ(LS2D.opp([5,8,9]),:,1) = fInC(LS2D.opp([5,8,9]),1:(LatticeC.Nx+1)/2,LatticeC.Ny);
    

%     fOutC(:,1,(LatticeC.Ny+3)/2:LatticeC.Ny-1) = fOutP(:,LatticeP.Nx,2:end-1);
%     fOutC(:,2:(LatticeC.Nx-1)/2,LatticeC.Ny-1) = fOutQ(:,2:end-1,1);
%     [fInC,LatticeC] = Streaming2D(fInC,fOutC,LatticeC,LS2D,BCsC,'C');    

%     fOutP([4,7,8],LatticeP.Nx,2:end-1) = fOutC([4,7,8],1,(LatticeC.Ny+3)/2:LatticeC.Ny-1); 
%     [fInP,LatticeP] = Streaming2D(fInP,fOutP,LatticeP,LS2D,BCsP,'P');
%     
%     fOutQ([3,6,7],2:end-1,1) = fOutC([3,6,7],2:(LatticeC.Nx-1)/2,LatticeC.Ny-1);    
%     [fInQ,LatticeQ] = Streaming2D(fInQ,fOutQ,LatticeQ,LS2D,BCsQ,'Q');
%     
%     fInC([2,6,9],1,(LatticeC.Ny+3)/2:LatticeC.Ny-1) = fInP([2,6,9],LatticeP.Nx,2:end-1);  
%     fInP(LS2D.opp([2,6,9]),LatticeP.Nx,2:end-1) = fInC(LS2D.opp([2,6,9]),1,(LatticeC.Ny+3)/2:LatticeC.Ny-1);
%     fInC([5,8,9],2:(LatticeC.Nx-1)/2,LatticeC.Ny) = fInQ([5,8,9],2:end-1,1);    
%     fInQ(LS2D.opp([5,8,9]),2:end-1,1) = fInC(LS2D.opp([5,8,9]),2:(LatticeC.Nx-1)/2,LatticeC.Ny);

    %% 10. Visualization
    % available variables, {pressure,velocity}
    if (mod(cycle,PlotStep)==1)
        Visualization2D(LatticeP,CharP,LatticeQ,CharQ,LatticeC,CharC,'velocity')
%         Visualization2D(LatticeP,CharP,LatticeQ,CharQ,LatticeC,CharC,'pressure')

        cycle
    end
    %%
    cycle = cycle + 1;
end 
