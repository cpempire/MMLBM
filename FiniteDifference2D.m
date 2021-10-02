function fIn = FiniteDifference2D(Lattice,BCs,fIn,fEq,LS2D,Location)
%% Regularized Boundary Dynamics
in = 1;  out = Lattice.Nx;       % in and out boundary x-direction indices
up = Lattice.Ny; down = 1;       % up and down wall boundary y-direction indices

switch Location
    case 'Boundary'
         wall = BCs.B;
    case 'Left'
         wall = BCs.L;          
    case 'Right'
         wall = BCs.R;
    case 'Up'
         wall = BCs.U;
    case 'Down'
         wall = BCs.D;
    otherwise 
        disp('please choose Boundary, Left, Right,Up,Down ')
end 

% equilibrium distribution
for i=1:9
    cu = 3*(LS2D.cx(i)*Lattice.Ux(:,wall)+LS2D.cy(i)*Lattice.Uy(:,wall));
    fEq(i,wall) = Lattice.Rho(:,wall)*LS2D.w(i) .* ...
        ( 1  + cu + 1/2*(cu.*cu) - 3/2*(Lattice.Ux(:,wall).^2+Lattice.Uy(:,wall).^2) );    
end

S = zeros(2,2,length(wall));

switch Location 
    case 'Boundary'
        % Down
        S(1,1,2:Lattice.Nx-1) = (Lattice.Ux(:,3:end,down) - Lattice.Ux(:,1:end-2,down))/2;
        S(1,1,1)              = (-3*Lattice.Ux(:,in,down)+4*Lattice.Ux(:,in+1,down)-Lattice.Ux(:,in+2,down))/2;
        S(1,1,Lattice.Nx)     = (-3*Lattice.Ux(:,out,down)+4*Lattice.Ux(:,out-1,down)-Lattice.Ux(:,out-2,down))/2;
        S(1,2,1:Lattice.Nx)   = (-3*Lattice.Ux(:,:,down)+4*Lattice.Ux(:,:,down+1)-Lattice.Ux(:,:,down+2))/2;
        
        S(2,1,2:Lattice.Nx-1) = (Lattice.Uy(:,3:end,down) - Lattice.Uy(:,1:end-2,down))/2;
        S(2,1,1)              = (-3*Lattice.Uy(:,in,down)+4*Lattice.Uy(:,in+1,down)-Lattice.Uy(:,in+2,down))/2;
        S(2,1,Lattice.Nx)     = (-3*Lattice.Uy(:,out,down)+4*Lattice.Uy(:,out-1,down)-Lattice.Uy(:,out-2,down))/2;
        S(2,2,1:Lattice.Nx)   = (-3*Lattice.Uy(:,:,down)+4*Lattice.Uy(:,:,down+1)-Lattice.Uy(:,:,down+2))/2;
        
        % Left
        
        S(1,2,Lattice.Nx+1:2:Lattice.Nx+1+2*(Lattice.Ny-3)) = (Lattice.Ux(:,in,3:end) - Lattice.Ux(:,in,1:end-2))/2;
        S(1,1,Lattice.Nx+1:2:Lattice.Nx+1+2*(Lattice.Ny-3)) = (-3*Lattice.Ux(:,in,2:end-1)+4*Lattice.Ux(:,in+1,2:end-1)-Lattice.Ux(:,in+2,2:end-1))/2;
        
        S(2,2,Lattice.Nx+1:2:Lattice.Nx+1+2*(Lattice.Ny-3)) = (Lattice.Uy(:,in,3:end) - Lattice.Uy(:,in,1:end-2))/2;
        S(2,1,Lattice.Nx+1:2:Lattice.Nx+1+2*(Lattice.Ny-3)) = (-3*Lattice.Uy(:,in,2:end-1)+4*Lattice.Uy(:,in+1,2:end-1)-Lattice.Uy(:,in+2,2:end-1))/2;         

        % Right
        
        S(1,2,Lattice.Nx+2:2:Lattice.Nx+2+2*(Lattice.Ny-3)) = (Lattice.Ux(:,out,3:end) - Lattice.Ux(:,out,1:end-2))/2;
        S(1,1,Lattice.Nx+2:2:Lattice.Nx+2+2*(Lattice.Ny-3)) = (-3*Lattice.Ux(:,out,2:end-1)+4*Lattice.Ux(:,out-1,2:end-1)-Lattice.Ux(:,out-2,2:end-1))/2;
        
        S(2,2,Lattice.Nx+2:2:Lattice.Nx+2+2*(Lattice.Ny-3)) = (Lattice.Uy(:,out,3:end) - Lattice.Uy(:,out,1:end-2))/2;
        S(2,1,Lattice.Nx+2:2:Lattice.Nx+2+2*(Lattice.Ny-3)) = (-3*Lattice.Uy(:,out,2:end-1)+4*Lattice.Uy(:,out-1,2:end-1)-Lattice.Uy(:,out-2,2:end-1))/2;   
        
        % Up
        S(1,1,end-Lattice.Nx+2:end-1) = (Lattice.Ux(:,3:end,up) - Lattice.Ux(:,1:end-2,up))/2;
        S(1,1,end-Lattice.Nx+1)       = (-3*Lattice.Ux(:,in,up)+4*Lattice.Ux(:,in+1,up)-Lattice.Ux(:,in+2,up))/2;
        S(1,1,end)                    = (-3*Lattice.Ux(:,out,up)+4*Lattice.Ux(:,out-1,up)-Lattice.Ux(:,out-2,up))/2;
        S(1,2,end-Lattice.Nx+1:end)   = (-3*Lattice.Ux(:,:,up)+4*Lattice.Ux(:,:,up-1)-Lattice.Ux(:,:,up-2))/2;
        
        S(2,1,end-Lattice.Nx+2:end-1) = (Lattice.Uy(:,3:end,up) - Lattice.Uy(:,1:end-2,up))/2;
        S(2,1,end-Lattice.Nx+1)       = (-3*Lattice.Uy(:,in,up)+4*Lattice.Uy(:,in+1,up)-Lattice.Uy(:,in+2,up))/2;
        S(2,1,end)                    = (-3*Lattice.Uy(:,out,up)+4*Lattice.Uy(:,out-1,up)-Lattice.Uy(:,out-2,up))/2;
        S(2,2,end-Lattice.Nx+1:end)   = (-3*Lattice.Uy(:,:,up)+4*Lattice.Uy(:,:,up-1)-Lattice.Uy(:,:,up-2))/2; 
             
    case 'Left'

        S(1,2,2:end-1) = (Lattice.Ux(:,in,3:end) - Lattice.Ux(:,in,1:end-2))/2;
        S(1,2,1)       = (-3*Lattice.Ux(:,in,down)+4*Lattice.Ux(:,in,down+1)-Lattice.Ux(:,in,down+2))/2;
        S(1,2,end)     = (-3*Lattice.Ux(:,in,up)+4*Lattice.Ux(:,in,up-1)-Lattice.Ux(:,in,up-2))/2;
        S(1,1,:)       = (-3*Lattice.Ux(:,in,:)+4*Lattice.Ux(:,in+1,:)-Lattice.Ux(:,in+2,:))/2;
        
        S(2,2,2:end-1) = (Lattice.Uy(:,in,3:end) - Lattice.Uy(:,in,1:end-2))/2;
        S(2,2,1)       = (-3*Lattice.Uy(:,in,down)+4*Lattice.Uy(:,in,down+1)-Lattice.Uy(:,in,down+2))/2;
        S(2,2,end)     = (-3*Lattice.Uy(:,in,up)+4*Lattice.Uy(:,in,up-1)-Lattice.Uy(:,in,up-2))/2;
        S(2,1,:)       = (-3*Lattice.Uy(:,in,:)+4*Lattice.Uy(:,in+1,:)-Lattice.Uy(:,in+2,:))/2; 
      
    case 'Right'
        
        S(1,2,2:end-1) = (Lattice.Ux(:,out,3:end) - Lattice.Ux(:,out,1:end-2))/2;
        S(1,2,1)       = (-3*Lattice.Ux(:,out,down)+4*Lattice.Ux(:,out,down+1)-Lattice.Ux(:,out,down+2))/2;
        S(1,2,end)     = (-3*Lattice.Ux(:,out,up)+4*Lattice.Ux(:,out,up-1)-Lattice.Ux(:,out,up-2))/2;
        S(1,1,:)       = (-3*Lattice.Ux(:,out,:)+4*Lattice.Ux(:,out-1,:)-Lattice.Ux(:,out-2,:))/2;
        
        S(2,2,2:end-1) = (Lattice.Uy(:,out,3:end) - Lattice.Uy(:,out,1:end-2))/2;
        S(2,2,1)       = (-3*Lattice.Uy(:,out,down)+4*Lattice.Uy(:,out,down+1)-Lattice.Uy(:,out,down+2))/2;
        S(2,2,end)     = (-3*Lattice.Uy(:,out,up)+4*Lattice.Uy(:,out,up-1)-Lattice.Uy(:,out,up-2))/2;
        S(2,1,:)       = (-3*Lattice.Uy(:,out,:)+4*Lattice.Uy(:,out-1,:)-Lattice.Uy(:,out-2,:))/2; 
        
    case 'Up'
        
        S(1,1,2:end-1) = (Lattice.Ux(:,3:end,up) - Lattice.Ux(:,1:end-2,up))/2;
        S(1,1,1)       = (-3*Lattice.Ux(:,in,up)+4*Lattice.Ux(:,in+1,up)-Lattice.Ux(:,in+2,up))/2;
        S(1,1,end)     = (-3*Lattice.Ux(:,out,up)+4*Lattice.Ux(:,out-1,up)-Lattice.Ux(:,out-2,up))/2;
        S(1,2,:)       = (-3*Lattice.Ux(:,:,up)+4*Lattice.Ux(:,:,up-1)-Lattice.Ux(:,:,up-2))/2;
        
        S(2,1,2:end-1) = (Lattice.Uy(:,3:end,up) - Lattice.Uy(:,1:end-2,up))/2;
        S(2,1,1)       = (-3*Lattice.Uy(:,in,up)+4*Lattice.Uy(:,in+1,up)-Lattice.Uy(:,in+2,up))/2;
        S(2,1,end)     = (-3*Lattice.Uy(:,out,up)+4*Lattice.Uy(:,out-1,up)-Lattice.Uy(:,out-2,up))/2;
        S(2,2,:)       = (-3*Lattice.Uy(:,:,up)+4*Lattice.Uy(:,:,up-1)-Lattice.Uy(:,:,up-2))/2;        
        
    case 'Down'
        
        S(1,1,2:end-1) = (Lattice.Ux(:,3:end,down) - Lattice.Ux(:,1:end-2,down))/2;
        S(1,1,1)       = (-3*Lattice.Ux(:,in,down)+4*Lattice.Ux(:,in+1,down)-Lattice.Ux(:,in+2,down))/2;
        S(1,1,end)     = (-3*Lattice.Ux(:,out,down)+4*Lattice.Ux(:,out-1,down)-Lattice.Ux(:,out-2,down))/2;
        S(1,2,:)       = (-3*Lattice.Ux(:,:,down)+4*Lattice.Ux(:,:,down+1)-Lattice.Ux(:,:,down+2))/2;
        
        S(2,1,2:end-1) = (Lattice.Uy(:,3:end,down) - Lattice.Uy(:,1:end-2,down))/2;
        S(2,1,1)       = (-3*Lattice.Uy(:,in,down)+4*Lattice.Uy(:,in+1,down)-Lattice.Uy(:,in+2,down))/2;
        S(2,1,end)     = (-3*Lattice.Uy(:,out,down)+4*Lattice.Uy(:,out-1,down)-Lattice.Uy(:,out-2,down))/2;
        S(2,2,:)       = (-3*Lattice.Uy(:,:,down)+4*Lattice.Uy(:,:,down+1)-Lattice.Uy(:,:,down+2))/2;
        
    otherwise 
         disp('please choose Boundary, Left, Right,Up,Down ')
end

TempNEq= (LS2D.wQ*reshape(S,4,length(wall)));
for i = 1:length(wall)
fIn(:,wall(i)) = fEq(:,wall(i)) - ...
                           3*Lattice.Tau*Lattice.Rho(:,wall(i)).*TempNEq(:,i);
end
