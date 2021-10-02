function fIn = Regularized2D(Lattice,BCs,fIn,fEq,LS2D,Location)
%% Regularized Boundary Dynamics
in = 1;  out = Lattice.Nx;       % in and out boundary x-direction indices
col = 2:Lattice.Ny-1;            % in and out boundary y-direction indices
up = Lattice.Ny; down = 1;       % up and down wall boundary y-direction indices
row = 2:Lattice.Nx-1;            % up and down wall boundary x-direction indices

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
    case 'Curve'
         wall = Lattice.CurveBoundary;
    otherwise 
        disp('please choose Boundary, Left, Right,Up,Down ')
end 

% equilibrium distribution
for i=1:9
    cu = 3*(LS2D.cx(i)*Lattice.Ux(:,wall)+LS2D.cy(i)*Lattice.Uy(:,wall));
    fEq(i,wall) = Lattice.Rho(:,wall)*LS2D.w(i) .* ...
        ( 1  + cu + 1/2*(cu.*cu) - 3/2*(Lattice.Ux(:,wall).^2+Lattice.Uy(:,wall).^2) );    
end

% unknown distribution function on the boundary 
% obtained by non-equilibrium bounceback rule
% Left and Right
switch Location 
    case 'Boundary'
        % Left and Right
        fIn([2,6,9],in,col) = fEq([2,6,9],in,col) + ...
          fIn(LS2D.opp([2,6,9]),in, col) - fEq(LS2D.opp([2,6,9]),in, col);

        fIn(LS2D.opp([2,6,9]),out,col) = fEq(LS2D.opp([2,6,9]),out,col) + ...
          fIn([2,6,9],out, col) - fEq([2,6,9],out, col);  
        % Up and Down 
        fIn([5,8,9],row,up) = fEq([5,8,9],row,up) + ...
          fIn(LS2D.opp([5,8,9]),row, up) - fEq(LS2D.opp([5,8,9]),row, up);

        fIn(LS2D.opp([5,8,9]),row,down) = fEq(LS2D.opp([5,8,9]),row,down) + ...
          fIn([5,8,9],row, down) - fEq([5,8,9],row, down);
        % Corner Left-Up and Right-Down
        fIn([2,5,9],in,up) = fEq([2,5,9],in,up) + ...
          fIn(LS2D.opp([2,5,9]),in,up) - fEq(LS2D.opp([2,5,9]),in,up);

        fIn(LS2D.opp([2,5,9]),out,down) = fEq(LS2D.opp([2,5,9]),out,down) + ...
          fIn([2,5,9],out,down) - fEq([2,5,9],out,down);
        % Corner Left-Down and Right-Up
        fIn([2,3,6],in,down) = fEq([2,3,6],in,down) + ...
          fIn(LS2D.opp([2,3,6]),in,down) - fEq(LS2D.opp([2,3,6]),in,down);

        fIn(LS2D.opp([2,3,6]),out,up) = fEq(LS2D.opp([2,3,6]),out,up) + ...
          fIn([2,3,6],out,up) - fEq([2,3,6],out,up);   

    case 'Left'
        fIn([2,6,9],in,col) = fEq([2,6,9],in,col) + ...
          fIn(LS2D.opp([2,6,9]),in, col) - fEq(LS2D.opp([2,6,9]),in, col);        
        
        fIn([2,5,9],in,up) = fEq([2,5,9],in,up) + ...
          fIn(LS2D.opp([2,5,9]),in,up) - fEq(LS2D.opp([2,5,9]),in,up);        
        fIn([2,3,6],in,down) = fEq([2,3,6],in,down) + ...
          fIn(LS2D.opp([2,3,6]),in,down) - fEq(LS2D.opp([2,3,6]),in,down);
      
    case 'Right'
        fIn(LS2D.opp([2,6,9]),out,col) = fEq(LS2D.opp([2,6,9]),out,col) + ...
          fIn([2,6,9],out, col) - fEq([2,6,9],out, col);          
        
        fIn(LS2D.opp([2,5,9]),out,down) = fEq(LS2D.opp([2,5,9]),out,down) + ...
          fIn([2,5,9],out,down) - fEq([2,5,9],out,down);        
        fIn(LS2D.opp([2,3,6]),out,up) = fEq(LS2D.opp([2,3,6]),out,up) + ...
          fIn([2,3,6],out,up) - fEq([2,3,6],out,up);         
    case 'Up'
        fIn([5,8,9],row,up) = fEq([5,8,9],row,up) + ...
          fIn(LS2D.opp([5,8,9]),row, up) - fEq(LS2D.opp([5,8,9]),row, up);        
        
        fIn([2,5,9],in,up) = fEq([2,5,9],in,up) + ...
          fIn(LS2D.opp([2,5,9]),in,up) - fEq(LS2D.opp([2,5,9]),in,up);        
        fIn(LS2D.opp([2,3,6]),out,up) = fEq(LS2D.opp([2,3,6]),out,up) + ...
          fIn([2,3,6],out,up) - fEq([2,3,6],out,up);          
    case 'Down'
        fIn(LS2D.opp([5,8,9]),row,down) = fEq(LS2D.opp([5,8,9]),row,down) + ...
          fIn([5,8,9],row, down) - fEq([5,8,9],row, down);        
        
        fIn(LS2D.opp([2,5,9]),out,down) = fEq(LS2D.opp([2,5,9]),out,down) + ...
          fIn([2,5,9],out,down) - fEq([2,5,9],out,down);        
        fIn([2,3,6],in,down) = fEq([2,3,6],in,down) + ...
          fIn(LS2D.opp([2,3,6]),in,down) - fEq(LS2D.opp([2,3,6]),in,down);        
    case 'Curve'
        for i = 1:length(Lattice.CurveBoundary)
           direction = find(BCs.C(1:9,Lattice.CurveBoundary(i))==1);
           fIn(LS2D.opp(direction),Lattice.CurveBoundary(i)) = ...
               fEq(LS2D.opp(direction),Lattice.CurveBoundary(i)) ...
               + fIn(direction,Lattice.CurveBoundary(i)) - fEq(direction,Lattice.CurveBoundary(i)); 
        end
    otherwise 
         disp('please choose Boundary, Left, Right,Up,Down ')
end 

% non-equilibrium distribution function on the boundary 
% second order tensor P corresponding to non-equilibrium part
% and reconstruct distribution as fIn = fEq + fnEq on the boundary 
P = zeros(length(wall),2,2);
for i = 1:length(wall)
  P(i,:,:) = reshape(LS2D.Cxy*fIn(:,wall(i)),2,2) - ... 
      1/3*Lattice.Rho(:,wall(i))*eye(2) - ...
      Lattice.Rho(:,wall(i))*[Lattice.Ux(:,wall(i));Lattice.Uy(:,wall(i))]*...
      [Lattice.Ux(:,wall(i)),Lattice.Uy(:,wall(i))]; 
  fIn(:,wall(i)) = fEq(:,wall(i)) + ...
      9/2*LS2D.w'.*(LS2D.Q*reshape(P(i,:,:),4,1)); % new distribution
end
