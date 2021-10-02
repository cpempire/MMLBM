function fIn = HalfWay2D(Lattice,fIn,LS2D,Location)
%% HalfWay2D Boundary Dynamics
in = 1;  out = Lattice.Nx;       % in and out boundary x-direction indices
col = 2:Lattice.Ny-1;            % in and out boundary y-direction indices
up = Lattice.Ny; down = 1;       % up and down wall boundary y-direction indices
row = 2:Lattice.Nx-1;            % up and down wall boundary x-direction indices

switch Location
    case 'Boundary'
        fIn(2,in,col) = fIn(4,in,col) + 2/3*Lattice.Rho(:,in,col).*Lattice.Ux(:,in,col); 
        fIn(6,in,col) = fIn(8,in,col)  + 1/6*Lattice.Rho(:,in,col).*Lattice.Uy(:,in,col) ...
                                        + 1/6*Lattice.Rho(:,in,col).*Lattice.Ux(:,in,col); 
        fIn(9,in,col) = fIn(7,in,col) - 1/6*Lattice.Rho(:,in,col).*Lattice.Uy(:,in,col) ...
                                        + 1/6*Lattice.Rho(:,in,col).*Lattice.Ux(:,in,col);          
        
        fIn(4,out,col) = fIn(2,out,col) - 2/3*Lattice.Rho(:,out,col).*Lattice.Ux(:,out,col); 
        fIn(8,out,col) = fIn(6,out,col)  - 1/6*Lattice.Rho(:,out,col).*Lattice.Uy(:,out,col) ...
                                          - 1/6*Lattice.Rho(:,out,col).*Lattice.Ux(:,out,col); 
        fIn(7,out,col) = fIn(9,out,col) + 1/6*Lattice.Rho(:,out,col).*Lattice.Uy(:,out,col) ...
                                          - 1/6*Lattice.Rho(:,out,col).*Lattice.Ux(:,out,col);         
        
        fIn(5,row,up) = fIn(3,row,up) - 2/3*Lattice.Rho(:,row,up).*Lattice.Uy(:,row,up); 
        fIn(8,row,up) = fIn(6,row,up)  - 1/6*Lattice.Rho(:,row,up).*Lattice.Ux(:,row,up) ...
                                        - 1/6*Lattice.Rho(:,row,up).*Lattice.Uy(:,row,up); 
        fIn(9,row,up) = fIn(7,row,up) + 1/6*Lattice.Rho(:,row,up).*Lattice.Ux(:,row,up) ...
                                        - 1/6*Lattice.Rho(:,row,up).*Lattice.Uy(:,row,up);
                                    
                                    
        fIn(3,row,down) = fIn(5,row,down) + 2/3*Lattice.Rho(:,row,down).*Lattice.Uy(:,row,down); 
        fIn(6,row,down) = fIn(8,row,down) + 1/6*Lattice.Rho(:,row,down).*Lattice.Ux(:,row,down) ...
                                        + 1/6*Lattice.Rho(:,row,down).*Lattice.Uy(:,row,down);  
        fIn(7,row,down) = fIn(9,row,down) - 1/6*Lattice.Rho(:,row,down).*Lattice.Ux(:,row,down) ...
                                        + 1/6*Lattice.Rho(:,row,down).*Lattice.Uy(:,row,down);                    
                                    
        % corners scheme 2
        fIn(:,in,up) = fIn(:,in+1,up);
        fIn(:,in,down) = fIn(:,in+1,down);                    
        fIn(:,out,up) = fIn(:,out-1,up);                     
        fIn(:,out,down) = fIn(:,out-1,down);
               
        fIn(:,in,up)  = fIn(LS2D.opp,in,up);
        fIn([7,9],in,up) = 1./2*(Lattice.Rho(:,in,up) - ...
                                    sum(fIn([1,2,3,4,5,6,8],in,up)));

        fIn(:,in,down) = fIn(LS2D.opp,in,down);
        fIn([6,8],in,down) = 1./2*(Lattice.Rho(:,in,down) - ...
                                    sum(fIn([1,2,3,4,5,7,9],in,down)));  
                                
        fIn(:,out,up) = fIn(LS2D.opp,out,up);
        fIn([6,8],out,up) = 1./2*(Lattice.Rho(:,out,up) - ...
                                    sum(fIn([1,2,3,4,5,7,9],out,up)));

        fIn(:,out,down) = fIn(LS2D.opp,out,down);
        fIn([7,9],out,down) = 1./2*(Lattice.Rho(:,out,down) - ...
                                    sum(fIn([1,2,3,4,5,6,8],out,down)));   
                                
    
    case 'Left' 

        fIn(2,in,col) = fIn(4,in,col) + 2/3*Lattice.Rho(:,in,col).*Lattice.Ux(:,in,col); 
        fIn(6,in,col) = fIn(8,in,col)  + 1/6*Lattice.Rho(:,in,col).*Lattice.Uy(:,in,col) ...
                                        + 1/6*Lattice.Rho(:,in,col).*Lattice.Ux(:,in,col); 
        fIn(9,in,col) = fIn(7,in,col) - 1/6*Lattice.Rho(:,in,col).*Lattice.Uy(:,in,col) ...
                                        + 1/6*Lattice.Rho(:,in,col).*Lattice.Ux(:,in,col);   
                                    
        fIn(:,in,up)  = fIn(LS2D.opp,in,up);
        fIn([7,9],in,up) = 1./2*(Lattice.Rho(:,in,up) - ...
                                    sum(fIn([1,2,3,4,5,6,8],in,up)));

        fIn(:,in,down) = fIn(LS2D.opp,in,down);
        fIn([6,8],in,down) = 1./2*(Lattice.Rho(:,in,down) - ...
                                    sum(fIn([1,2,3,4,5,7,9],in,down)));                                    

    case 'Right' 

        fIn(4,out,col) = fIn(2,out,col) - 2/3*Lattice.Rho(:,out,col).*Lattice.Ux(:,out,col); 
        fIn(8,out,col) = fIn(6,out,col)  - 1/6*Lattice.Rho(:,out,col).*Lattice.Uy(:,out,col) ...
                                          - 1/6*Lattice.Rho(:,out,col).*Lattice.Ux(:,out,col); 
        fIn(7,out,col) = fIn(9,out,col) + 1/6*Lattice.Rho(:,out,col).*Lattice.Uy(:,out,col) ...
                                          - 1/6*Lattice.Rho(:,out,col).*Lattice.Ux(:,out,col);            
                                      
        fIn(:,out,up) = fIn(LS2D.opp,out,up);
        fIn([6,8],out,up) = 1./2*(Lattice.Rho(:,out,up) - ...
                                    sum(fIn([1,2,3,4,5,7,9],out,up)));

        fIn(:,out,down) = fIn(LS2D.opp,out,down);
        fIn([7,9],out,down) = 1./2*(Lattice.Rho(:,out,down) - ...
                                    sum(fIn([1,2,3,4,5,6,8],out,down)));                                          
                                      

    case 'Up'

        fIn(5,row,up) = fIn(3,row,up) - 2/3*Lattice.Rho(:,row,up).*Lattice.Uy(:,row,up); 
        fIn(8,row,up) = fIn(6,row,up)  - 1/6*Lattice.Rho(:,row,up).*Lattice.Ux(:,row,up) ...
                                        - 1/6*Lattice.Rho(:,row,up).*Lattice.Uy(:,row,up); 
        fIn(9,row,up) = fIn(7,row,up) + 1/6*Lattice.Rho(:,row,up).*Lattice.Ux(:,row,up) ...
                                        - 1/6*Lattice.Rho(:,row,up).*Lattice.Uy(:,row,up);
                                    
        fIn(:,in,up) = fIn(LS2D.opp,in,up);
        fIn([7,9],in,up) = 1./2*(Lattice.Rho(:,in,up) - ...
                                    sum(fIn([1,2,3,4,5,6,8],in,up)));

        fIn(:,out,up) = fIn(LS2D.opp,out,up);
        fIn([6,8],out,up) = 1./2*(Lattice.Rho(:,out,up) - ...
                                    sum(fIn([1,2,3,4,5,7,9],out,up)));
    case 'Down'

        fIn(3,row,down) = fIn(5,row,down) + 2/3*Lattice.Rho(:,row,down).*Lattice.Uy(:,row,down); 
        fIn(6,row,down) = fIn(8,row,down) + 1/6*Lattice.Rho(:,row,down).*Lattice.Ux(:,row,down) ...
                                        + 1/6*Lattice.Rho(:,row,down).*Lattice.Uy(:,row,down);  
        fIn(7,row,down) = fIn(9,row,down) - 1/6*Lattice.Rho(:,row,down).*Lattice.Ux(:,row,down) ...
                                        + 1/6*Lattice.Rho(:,row,down).*Lattice.Uy(:,row,down);     
                                    
        fIn(:,in,down) = fIn(LS2D.opp,in,down);
        fIn([6,8],in,down) = 1./2*(Lattice.Rho(:,in,down) - ...
                                    sum(fIn([1,2,3,4,5,7,9],in,down)));

        fIn(:,out,down) = fIn(LS2D.opp,out,down);
        fIn([7,9],out,down) = 1./2*(Lattice.Rho(:,out,down) - ...
                                    sum(fIn([1,2,3,4,5,6,8],out,down)));                              
    otherwise 
        disp('please impose boundary conditions on the four edges')
end        
