function fIn = BoundaryDynamics2D(BoundaryDynamics2D,Lattice,BCs,fIn,fEq,LS2D,Location)
%% Boundary Dynamics 
switch BoundaryDynamics2D
    case 'ZouHe2D'
    %% ZouHe Boundary Dynamics 
    % Location can be chosen as 'Boundary' for all the boundary nodes
    % 'Left', 'Right','Up','Down' for their relavent edge nodes 
    switch Location
        case 'Boundary'
            fIn = ZouHe2D(Lattice,fIn,fEq,BCs,LS2D,'Boundary');
        case 'Left'
            fIn = ZouHe2D(Lattice,fIn,fEq,BCs,LS2D,'Left'); 
        case 'Right'
            fIn = ZouHe2D(Lattice,fIn,fEq,BCs,LS2D,'Right');  
        case 'Up'
            fIn = ZouHe2D(Lattice,fIn,fEq,BCs,LS2D,'Up'); 
        case 'Down'
            fIn = ZouHe2D(Lattice,fIn,fEq,BCs,LS2D,'Down');
        case 'Curve'
            fIn = ZouHe2D(Lattice,fIn,fEq,BCs,LS2D,'Curve');
        otherwise 
            disp('please choose Boundary, Left, Right,Up,Down ')
    end
    case 'HalfWay2D'
    %% HalfWay Boundary Dynamics 
    % Location can be chosen as 'Boundary' for all the boundary nodes
    % 'Left', 'Right','Up','Down' for their relavent edge nodes 
    switch Location
        case 'Boundary'
            fIn = HalfWay2D(Lattice,fIn,LS2D,'Boundary');
        case 'Left'
            fIn = HalfWay2D(Lattice,fIn,LS2D,'Left'); 
        case 'Right'
            fIn = HalfWay2D(Lattice,fIn,LS2D,'Right');  
        case 'Up'
            fIn = HalfWay2D(Lattice,fIn,LS2D,'Up'); 
        case 'Down'
            fIn = HalfWay2D(Lattice,fIn,LS2D,'Down');
        otherwise 
            disp('please choose Boundary, Left, Right,Up,Down ')
    end    
    case 'Regularized2D'
    %%  Regularized Boundary Dynamics
        % BCs.B is a vector of the position of all the boundary nodes, while
        % BCs.L,BCs.R,BCs.U,BCs.D are those for Left,Right,Up,Down parts
    switch Location
        case 'Boundary'
             fIn = Regularized2D(Lattice,BCs,fIn,fEq,LS2D,'Boundary');
        case 'Left'
             fIn = Regularized2D(Lattice,BCs,fIn,fEq,LS2D,'Left');
        case 'Right'
             fIn = Regularized2D(Lattice,BCs,fIn,fEq,LS2D,'Right');
        case 'Up'
             fIn = Regularized2D(Lattice,BCs,fIn,fEq,LS2D,'Up');
        case 'Down'
             fIn = Regularized2D(Lattice,BCs,fIn,fEq,LS2D,'Down');
        case 'Curve'
             fIn = Regularized2D(Lattice,BCs,fIn,fEq,LS2D,'Curve');
        otherwise 
            disp('please choose Boundary, Left, Right,Up,Down ')
    end 
    case 'FiniteDifference2D'
    %% Finite Difference Dynamics 
    % Location can be chosen as 'Boundary' for all the boundary nodes
    % 'Left', 'Right','Up','Down' for their relavent edge nodes 
    switch Location
        case 'Boundary'
            fIn = FiniteDifference2D(Lattice,BCs,fIn,fEq,LS2D,'Boundary');
        case 'Left'
            fIn = FiniteDifference2D(Lattice,BCs,fIn,fEq,LS2D,'Left'); 
        case 'Right'
            fIn = FiniteDifference2D(Lattice,BCs,fIn,fEq,LS2D,'Right');  
        case 'Up'
            fIn = FiniteDifference2D(Lattice,BCs,fIn,fEq,LS2D,'Up'); 
        case 'Down'
            fIn = FiniteDifference2D(Lattice,BCs,fIn,fEq,LS2D,'Down');
        otherwise 
            disp('please choose Boundary, Left, Right,Up,Down ')
    end    
    otherwise
        disp('please choose ZouHe2D,HalfWay2D, Regularized2D or FiniteDifference2D')     
end
    