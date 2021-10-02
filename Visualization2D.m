function Visualization2D(LatticeP,CharP,LatticeQ,CharQ,LatticeC,CharC,Variable)
%% Visualization  
switch Variable
    case 'pressure'
        imagesc(LatticeP.gridx,LatticeP.gridy,flipud(reshape((LatticeP.Rho-CharP.DimensionlessRho)/3*CharP.DimensionlessUMax^2*CharP.PhysicalUMax^2,LatticeP.Nx,LatticeP.Ny)'));
        hold on 
        imagesc(LatticeQ.gridx,LatticeQ.gridy,flipud(reshape((LatticeQ.Rho-CharQ.DimensionlessRho)/3*CharQ.DimensionlessUMax^2*CharQ.PhysicalUMax^2,LatticeQ.Nx,LatticeQ.Ny)'));
        hold on 
        imagesc(LatticeC.gridx,LatticeC.gridy,flipud(reshape((LatticeC.Rho-CharC.DimensionlessRho)/3*CharC.DimensionlessUMax^2*CharC.PhysicalUMax^2,LatticeC.Nx,LatticeC.Ny)'));
        axis equal off; 
        colorbar
        drawnow 
    case 'velocity'
        imagesc(LatticeC.gridx,LatticeC.gridy,flipud(reshape(sqrt((LatticeC.Ux.^2+LatticeC.Uy.^2))*CharC.DimensionlessUMax*CharC.PhysicalUMax,LatticeC.Nx,LatticeC.Ny)'));
        hold on
        imagesc(LatticeP.gridx,LatticeP.gridy,flipud(reshape(sqrt((LatticeP.Ux.^2+LatticeP.Uy.^2))*CharP.DimensionlessUMax*CharP.PhysicalUMax,LatticeP.Nx,LatticeP.Ny)')); 
        hold on
        imagesc(LatticeQ.gridx,LatticeQ.gridy,flipud(reshape(sqrt((LatticeQ.Ux.^2+LatticeQ.Uy.^2))*CharQ.DimensionlessUMax*CharQ.PhysicalUMax,LatticeQ.Nx,LatticeQ.Ny)')); 
        axis equal off; 
        colorbar
        drawnow
end        

