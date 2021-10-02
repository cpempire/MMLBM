function Lattice = BoundaryConditions2D(Lattice,fIn,BCs,Location)
%% Zou/He Boundary Dynamics
in = 1;  out = Lattice.Nx;       % in and out boundary x-direction indices
col = 2:Lattice.Ny-1;            % in and out boundary y-direction indices
up = Lattice.Ny; down = 1;       % up and down wall boundary y-direction indices
row = 2:Lattice.Nx-1;            % up and down wall boundary x-direction indices


switch Location
    case 'Boundary'
        switch BCs.Left
            case 'velocity'
                Lattice.Ux(:,in,:)     = BCs.UxLeftStable;%BCs.UxLeft(in,Lattice.gridy,Lattice.UMax);
                Lattice.Uy(:,in,:)     = BCs.UyLeftStable;%BCs.UyLeft(in,Lattice.gridy,Lattice.UMax);
                Lattice.Rho(:,in,col)  = 1./(1-Lattice.Ux(:,in,col)).*( ...
                    sum(fIn([1,3,5],in,col)) + 2*sum(fIn([4,7,8],in,col)) ); 
                Lattice.Rho(:,in,up)   = 1/2*(Lattice.Rho(:,in,up-1)+Lattice.Rho(:,in+1,up));
                Lattice.Rho(:,in,down) = 1/2*(Lattice.Rho(:,in,down+1)+Lattice.Rho(:,in+1,down));
            case 'pressure'
                Lattice.Rho(:,in,:)    = BCs.RhoLeftStable;%BCs.RhoLeft(out,Lattice.gridy);
                Lattice.Ux(:,in,col)   = 1-1./Lattice.Rho(:,in,col).*( ...
                    sum(fIn([1,3,5],in,col)) + 2*sum(fIn([4,7,8],in,col)) );                 
                Lattice.Uy(:,in,:)     = BCs.UyLeftStable;%BCs.UyLeft(out,Lattice.gridy,Lattice.UMax);   
            otherwise
                disp('please impose known velocity or pressure')
        end        
        
        switch BCs.Right
            case 'velocity'
                Lattice.Ux(:,out,:)     = BCs.UxRightStable;%BCs.UxRight(out,Lattice.gridy,Lattice.UMax);
                Lattice.Uy(:,out,:)     = BCs.UyRightStable;%BCs.UyRight(out,Lattice.gridy,Lattice.UMax);
                Lattice.Rho(:,out,col)  = 1./(1+Lattice.Ux(:,out,col)).*( ...
                    sum(fIn([1,3,5],out,col)) + 2*sum(fIn([2,6,9],out,col)) );                 
                Lattice.Rho(:,out,up)   = 1/2*(Lattice.Rho(:,in,up-1)+Lattice.Rho(:,in+1,up));
                Lattice.Rho(:,out,down) = 1/2*(Lattice.Rho(:,in,down+1)+Lattice.Rho(:,in+1,down));             
            case 'pressure'
                Lattice.Rho(:,out,:)  = BCs.RhoRightStable;%BCs.RhoRight(out,Lattice.gridy);
                Lattice.Ux(:,out,col) = -1 + 1 ./ (Lattice.Rho(:,out,col)) .* ( ...
                    sum(fIn([1,3,5],out,col)) + 2*sum(fIn([2,6,9],out,col)) );
                Lattice.Uy(:,out,:)   = BCs.UyRightStable;%BCs.UyRight(out,Lattice.gridy,Lattice.UMax);      
            otherwise
                disp('please impose known velocity or pressure')
        end    
        
        switch BCs.Up
            case 'velocity'
            Lattice.Ux(:,:,up)     = BCs.UxUpStable;%BCs.UxUp(Lattice.gridx,up,Lattice.UMax);
            Lattice.Uy(:,:,up)     = BCs.UyUpStable;%BCs.UyUp(Lattice.gridx,up,Lattice.UMax);
            Lattice.Rho(:,row,up)  = 1./(1+Lattice.Uy(:,row,up)).*( ...
                sum(fIn([1,2,4],row,up)) + 2*sum(fIn([3,6,7],row,up)) );                 
            Lattice.Rho(:,in,up)   = 1/2*(Lattice.Rho(:,in,up-1)+Lattice.Rho(:,in+1,up));
            Lattice.Rho(:,out,up)  = 1/2*(Lattice.Rho(:,out,up-1)+Lattice.Rho(:,out-1,up));                      
            case 'pressure'
                Lattice.Rho(:,:,up)    = BCs.RhoUpStable;%BCs.RhoUp(out,Lattice.gridy);
                Lattice.Uy(:,row,up)   = -1 + 1./(Lattice.Rho(:,row,up)).*( ...
                    sum(fIn([1,2,4],row,up)) + 2*sum(fIn([3,6,7],row,up)) ); 
                Lattice.Ux(:,:,up)     = BCs.UxUpStable;%BCs.UxUp(out,Lattice.gridy,Lattice.UMax);     
            otherwise 
                disp('please impose known velocity or pressure')
        end
        
        switch BCs.Down
            case 'velocity'
                Lattice.Ux(:,:,down)     = BCs.UxDownStable;%BCs.UxUp(Lattice.gridx,down,Lattice.UMax);
                Lattice.Uy(:,:,down)     = BCs.UyDownStable;%BCs.UyUp(Lattice.gridx,down,Lattice.UMax);
                Lattice.Rho(:,row,down)  = 1./(1-Lattice.Uy(:,row,down)).*( ...
                    sum(fIn([1,2,4],row,down)) + 2*sum(fIn([5,8,9],row,down)) );                 
                Lattice.Rho(:,in,down)   = 1/2*(Lattice.Rho(:,in,down+1)+Lattice.Rho(:,in+1,down));
                Lattice.Rho(:,out,down)  = 1/2*(Lattice.Rho(:,out,down+1)+Lattice.Rho(:,out-1,down));      
            case 'pressure'
                Lattice.Rho(:,:,down)    = BCs.RhoDownStable;%BCs.RhoDown(out,Lattice.gridy);
                Lattice.Uy(:,row,down)   = 1 - 1./(Lattice.Rho(:,row,down)).*( ...
                    sum(fIn([1,2,4],row,down)) + 2*sum(fIn([5,8,9],row,down)) );    
                Lattice.Ux(:,:,down)     = BCs.UxDownStable;%BCs.UxDown(out,Lattice.gridy,Lattice.UMax);
            otherwise 
                disp('please impose known velocity or pressure')
        end        
        
    case 'Left'
        
        switch BCs.Left
            case 'velocity'
                Lattice.Ux(:,in,:)     = BCs.UxLeftStable;%BCs.UxLeft(in,Lattice.gridy,Lattice.UMax);
                Lattice.Uy(:,in,:)     = BCs.UyLeftStable;%BCs.UyLeft(in,Lattice.gridy,Lattice.UMax);
                Lattice.Rho(:,in,col)  = 1./(1-Lattice.Ux(:,in,col)).*( ...
                    sum(fIn([1,3,5],in,col)) + 2*sum(fIn([4,7,8],in,col)) ); 
                Lattice.Rho(:,in,up)   = 1/2*(Lattice.Rho(:,in,up-1)+Lattice.Rho(:,in+1,up));
                Lattice.Rho(:,in,down) = 1/2*(Lattice.Rho(:,in,down+1)+Lattice.Rho(:,in+1,down));
            case 'pressure'
                Lattice.Rho(:,in,:)    = BCs.RhoLeftStable;%BCs.RhoLeft(out,Lattice.gridy);
                Lattice.Ux(:,in,col)   = 1-1./Lattice.Rho(:,in,col).*( ...
                    sum(fIn([1,3,5],in,col)) + 2*sum(fIn([4,7,8],in,col)) );                 
                Lattice.Uy(:,in,:)     = BCs.UyLeftStable;%BCs.UyLeft(out,Lattice.gridy,Lattice.UMax);   
            otherwise
                disp('please impose known velocity or pressure')
        end
       
    case 'Right'
        switch BCs.Right
            case 'velocity'
                Lattice.Ux(:,out,:)     = BCs.UxRightStable;%BCs.UxRight(out,Lattice.gridy,Lattice.UMax);
                Lattice.Uy(:,out,:)     = BCs.UyRightStable;%BCs.UyRight(out,Lattice.gridy,Lattice.UMax);
                Lattice.Rho(:,out,col)  = 1./(1+Lattice.Ux(:,out,col)).*( ...
                    sum(fIn([1,3,5],out,col)) + 2*sum(fIn([2,6,9],out,col)) );                 
                Lattice.Rho(:,out,up)   = 1/2*(Lattice.Rho(:,in,up-1)+Lattice.Rho(:,in+1,up));
                Lattice.Rho(:,out,down) = 1/2*(Lattice.Rho(:,in,down+1)+Lattice.Rho(:,in+1,down));             
            case 'pressure'
                Lattice.Rho(:,out,:)  = BCs.RhoRightStable;%BCs.RhoRight(out,Lattice.gridy);
                Lattice.Ux(:,out,col) = -1 + 1 ./ (Lattice.Rho(:,out,col)) .* ( ...
                    sum(fIn([1,3,5],out,col)) + 2*sum(fIn([2,6,9],out,col)) );
                Lattice.Uy(:,out,:)   = BCs.UyRightStable;%BCs.UyRight(out,Lattice.gridy,Lattice.UMax);      
            otherwise
                disp('please impose known velocity or pressure')
        end

    case 'Up'
        switch BCs.Up
            case 'velocity'
            Lattice.Ux(:,:,up)     = BCs.UxUpStable;%BCs.UxUp(Lattice.gridx,up,Lattice.UMax);
            Lattice.Uy(:,:,up)     = BCs.UyUpStable;%BCs.UyUp(Lattice.gridx,up,Lattice.UMax);
            Lattice.Rho(:,row,up)  = 1./(1+Lattice.Uy(:,row,up)).*( ...
                sum(fIn([1,2,4],row,up)) + 2*sum(fIn([3,6,7],row,up)) );                 
            Lattice.Rho(:,in,up)   = 1/2*(Lattice.Rho(:,in,up-1)+Lattice.Rho(:,in+1,up));
            Lattice.Rho(:,out,up)  = 1/2*(Lattice.Rho(:,out,up-1)+Lattice.Rho(:,out-1,up));                      
            case 'pressure'
                Lattice.Rho(:,:,up)    = BCs.RhoUpStable;%BCs.RhoUp(out,Lattice.gridy);
                Lattice.Uy(:,row,up)   = -1 + 1./(Lattice.Rho(:,row,up)).*( ...
                    sum(fIn([1,2,4],row,up)) + 2*sum(fIn([3,6,7],row,up)) ); 
%                 Lattice.Uy(:,in,up)    = 1/2*(Lattice.Uy(:,in,up-1)+Lattice.Uy(:,in+1,up));
%                 Lattice.Uy(:,out,up)   = 1/2*(Lattice.Uy(:,out,up-1)+Lattice.Uy(:,out-1,up));
                Lattice.Ux(:,:,up)     = BCs.UxUpStable;%BCs.UxUp(out,Lattice.gridy,Lattice.UMax);     
            otherwise 
                disp('please impose known velocity or pressure')
        end
    
    case 'Down'
        switch BCs.Down
            case 'velocity'
                Lattice.Ux(:,:,down)     = BCs.UxDownStable;%BCs.UxUp(Lattice.gridx,down,Lattice.UMax);
                Lattice.Uy(:,:,down)     = BCs.UyDownStable;%BCs.UyUp(Lattice.gridx,down,Lattice.UMax);
                Lattice.Rho(:,row,down)  = 1./(1-Lattice.Uy(:,row,down)).*( ...
                    sum(fIn([1,2,4],row,down)) + 2*sum(fIn([5,8,9],row,down)) );                 
                Lattice.Rho(:,in,down)   = 1/2*(Lattice.Rho(:,in,down+1)+Lattice.Rho(:,in+1,down));
                Lattice.Rho(:,out,down)  = 1/2*(Lattice.Rho(:,out,down+1)+Lattice.Rho(:,out-1,down));      
            case 'pressure'
                Lattice.Rho(:,:,down)    = BCs.RhoDownStable;%BCs.RhoDown(out,Lattice.gridy);
                Lattice.Uy(:,row,down)   = 1 - 1./(Lattice.Rho(:,row,down)).*( ...
                    sum(fIn([1,2,4],row,down)) + 2*sum(fIn([5,8,9],row,down)) );    
                Lattice.Ux(:,:,down)     = BCs.UxDownStable;%BCs.UxDown(out,Lattice.gridy,Lattice.UMax);
            otherwise 
                disp('please impose known velocity or pressure')
        end
    case 'Curve'
        switch BCs.Curve
            case 'velocity'
                Lattice.Ux(:,Lattice.CurveBoundary) = 0;
                Lattice.Uy(:,Lattice.CurveBoundary) = 0;
                
            case 'pressure'
                disp('possible?')  
            otherwise 
                disp('please impose known velocity or pressure')                
        end
    otherwise 
        disp('please impose boundary conditions on the four edges')
end        
