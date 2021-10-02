function [Lattice,BCs,Char] = Discretization2D(NumCell,Order,Block,Dimensionless,BCs,Char,LS2D)
%% Discretization of the Dimensionless System to Lattice System 
Lattice.dx        = 1/NumCell; % spatial resolution
Lattice.gridx     = Block(1,1):Lattice.dx:Block(1,2);   % vector of position in x-direction
Lattice.gridy     = Block(2,1):Lattice.dx:Block(2,2);   % vector of position in y-direction

[Lattice.y,Lattice.x]     = meshgrid(Lattice.gridy,Lattice.gridx); % matrix of position for (x,y)
[Lattice.Nx,Lattice.Ny]   = size(Lattice.x);    % size of the matrix

Lattice.Lx        = Dimensionless.Lx;
Lattice.Ly        = Dimensionless.Ly;

switch Order
    case 2
        Lattice.dt        = Lattice.dx^2;  % temporal resolution, second order 
    case 1
        Lattice.dt        = Lattice.dx;  % temporal resolution, first order 
    case 0
        Lattice.dt        = 1e-6;
    otherwise
        disp('change temporal resulution by oneself here in Discreitzation')
end     

Lattice.Rho       = Dimensionless.Rho;
Lattice.Nu        = Dimensionless.Nu*Lattice.dt/Lattice.dx^2; % lattice viscosity 
Lattice.Tau       = 3*Lattice.Nu + 0.5;  % single relaxation time
Lattice.UMax      = Lattice.dt/Lattice.dx*Dimensionless.UMax;
Lattice.TimeMax   = Dimensionless.TimeMax/Lattice.dt;

Lattice.Fx        = Dimensionless.Fx(Lattice.x,Lattice.y,1,Lattice.UMax);
Lattice.Fy        = Dimensionless.Fy(Lattice.x,Lattice.y,1,Lattice.UMax);


Char.DimensionlessRho        = Lattice.Rho;
Char.DimensionlessLength     = Lattice.dx;
Char.DimensionlessTime       = Lattice.dt;
Char.DimensionlessUMax       = Lattice.dx/Lattice.dt;

% boundary indices for B = all boundary nodes, {L,R,U,D} = Left,Right,Up,Down boundary nodes 
BCs.B  =  find(Lattice.x==Block(1,1)|Lattice.x==Block(1,2)|Lattice.y==Block(2,1)|Lattice.y==Block(2,2));
BCs.L  =  find(Lattice.x==Block(1,1));
BCs.R  =  find(Lattice.x==Block(1,2));
BCs.U  =  find(Lattice.y==Block(2,1));      
BCs.D  =  find(Lattice.y==Block(2,2));



switch BCs.Left           
    case 'velocity'
        BCs.UxLeftStable  = BCs.UxLeft(:,Lattice.gridy,Lattice.UMax); 
        BCs.UyLeftStable  = BCs.UyLeft(:,Lattice.gridy,Lattice.UMax); 
    case 'pressure'
        BCs.RhoLeftStable = BCs.RhoLeft(:,Lattice.gridy,Lattice.UMax); 
        BCs.UyLeftStable  = BCs.UyLeft(:,Lattice.gridy,Lattice.UMax); 
end

switch BCs.Right           
    case 'velocity'
        BCs.UxRightStable  = BCs.UxRight(:,Lattice.gridy,Lattice.UMax); 
        BCs.UyRightStable  = BCs.UyRight(:,Lattice.gridy,Lattice.UMax); 
    case 'pressure'
        BCs.RhoRightStable = BCs.RhoRight(:,Lattice.gridy,Lattice.UMax); 
        BCs.UyRightStable  = BCs.UyRight(:,Lattice.gridy,Lattice.UMax); 
end

switch BCs.Up           
    case 'velocity'
        BCs.UxUpStable  = BCs.UxUp(Lattice.gridx,:,Lattice.UMax); 
        BCs.UyUpStable  = BCs.UyUp(Lattice.gridx,:,Lattice.UMax); 
    case 'pressure'
        BCs.RhoUpStable = BCs.RhoUp(Lattice.gridx,:,Lattice.UMax); 
        BCs.UxUpStable  = BCs.UxUp(Lattice.gridx,:,Lattice.UMax); 
end

switch BCs.Down           
    case 'velocity'
        BCs.UxDownStable  = BCs.UxDown(Lattice.gridx,:,Lattice.UMax); 
        BCs.UyDownStable  = BCs.UyDown(Lattice.gridx,:,Lattice.UMax); 
    case 'pressure'
        BCs.RhoDownStable = BCs.RhoDown(Lattice.gridx,:,Lattice.UMax); 
        BCs.UxDownStable  = BCs.UxDown(Lattice.gridx,:,Lattice.UMax); 
end

% switch BCs.Curve           
%     case 'velocity'
%         BCs.UxCurveStable  = BCs.UxCurve(:,Lattice.Boundary,Lattice.UMax); 
%         BCs.UyCurveStable  = BCs.UyCurve(:,Lattice.Boundary,Lattice.UMax); 
%     case 'pressure'
%         BCs.RhoLeftStable = BCs.RhoLeft(:,Lattice.gridy,Lattice.UMax); 
%         BCs.UyLeftStable  = BCs.UyLeft(:,Lattice.gridy,Lattice.UMax); 
% end


BCs.C  =  zeros(10,Lattice.Nx,Lattice.Ny);
Lattice.Region = zeros(Lattice.Nx,Lattice.Ny);

Lattice.Region(find(((Lattice.x-(Block(1,1)+Block(1,2))/2).^2+(Lattice.y-(Block(2,1)+Block(2,2))/2).^2<1) ...
                          |((Lattice.x<=(Block(1,1)+Block(1,2))/2)&(Lattice.y>=(Block(2,1)+Block(2,2))/2)))) = 1;
                      
Lattice.Obstacle = find(Lattice.Region==0);                      


for i = 1:Lattice.Nx
    for j = 1:Lattice.Ny
        if Lattice.Region(i,j)==1 
            for k = 1:9
                m = LS2D.cx(k);
                n = LS2D.cy(k);
                if i+m<=Lattice.Nx&&i+m>=1&&j+n<=Lattice.Ny&&j+n>=1
                    if Lattice.Region(i+m,j+n)==0
                        delta = 1/2;
                        for s = 1:10
                            xx = Lattice.x(i,j) + delta*LS2D.cx(k)*Lattice.dx;
                            yy = Lattice.y(i,j) + delta*LS2D.cy(k)*Lattice.dx;
                            pn = (xx-(Block(1,1)+Block(1,2))/2).^2+(yy-(Block(2,1)+Block(2,2))/2).^2>1;
                            delta = delta +1/2^(s+1) - pn*1/2^s;
                        end
                        BCs.C(k,i,j) = delta;
                        BCs.C(10,i,j) = BCs.C(10,i,j)+Lattice.Region(i+m,j+n);
                    end
                    BCs.C(10,i,j) = BCs.C(10,i,j)+Lattice.Region(i+m,j+n);
                else
%                     BCs.C(k,i,j) = 1;
                    BCs.C(10,i,j)=10;
                end
            end
            if BCs.C(10,i,j)<9
                BCs.C(10,i,j)=1;
            elseif  BCs.C(10,i,j)>10
                BCs.C(10,i,j)=2;
            end
        end   
    end
end
Lattice.CurveBoundary = find(BCs.C(10,:,:)==1);
Lattice.Boundary = find(reshape(BCs.C(10,:,:),Lattice.Nx,Lattice.Ny)~=9&Lattice.Region==1);
Lattice.InnerRegion = Lattice.Region;
Lattice.InnerRegion(Lattice.Boundary)=0;



% for i = 1:Lattice.Nx
%     for j = 1:Lattice.Ny
%         if Lattice.Region(i,j)==1 
%             for k = 1:9
%                 m = LS2D.cx(k);
%                 n = LS2D.cy(k);
%                 if i+m<=Lattice.Nx&&i+m>=1&&j+n<=Lattice.Ny&&j+n>=1
%                     if Lattice.Region(i+m,j+n)==0
%                         BCs.C(k,i,j)=1;
%                     end
%                     BCs.C(10,i,j) = BCs.C(10,i,j)+Lattice.Region(i+m,j+n);
%                 else
%                     BCs.C(k,i,j)=1;
%                 end
%             end
%             if BCs.C(10,i,j)<9
%                 BCs.C(10,i,j)=1;
%             else 
%                 BCs.C(10,i,j)=0;
%             end
%         end   
%     end
% end
% Lattice.Boundary = find(BCs.C(10,:,:)==1);
% Lattice.CurveBoundary = find(reshape(BCs.C(10,:,:),Lattice.Nx,Lattice.Ny)==1&Lattice.x~=Block(1,1)&Lattice.y~=Block(2,2));
% Lattice.InnerRegion = Lattice.Region;
% Lattice.InnerRegion(Lattice.Boundary)=0;




