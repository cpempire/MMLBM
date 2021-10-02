function fIn = Corner2D(Lattice,fIn,LS2D)
in = 1;  out = Lattice.Nx;       % in and out boundary x-direction indices
up = Lattice.Ny; down = 1;       % up and down wall boundary y-direction indices
%% Corners Left-Up Left-Down Right-Up Right-Down

% corners scheme 0
fIn(:,in,up) = fIn(LS2D.opp,in,up);
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
                                                
% % corners scheme 1
% fIn(:,in,up) = fIn(:,in+1,up-1);
% 
% fIn(:,in,down) = fIn(:,in+1,down+1);
%                         
% fIn(:,out,up) = fIn(:,out-1,up-1);
%                         
% fIn(:,out,down) = fIn(:,out-1,down+1);

% % corners scheme 2
% fIn(:,in,up) = fIn(:,in+1,up);
% 
% fIn(:,in,down) = fIn(:,in+1,down);
%                         
% fIn(:,out,up) = fIn(:,out-1,up);
%                         
% fIn(:,out,down) = fIn(:,out-1,down);
% 
% % corners scheme 3
% fIn(:,in,up) = 1/2*(fIn(:,in+1,up)+fIn(:,in,up-1));
% 
% fIn(:,in,down) = 1/2*(fIn(:,in+1,down)+fIn(:,in,down+1));
%                         
% fIn(:,out,up) = 1/2*(fIn(:,out-1,up)+fIn(:,out,up-1));
%                         
% fIn(:,out,down) = 1/2*(fIn(:,out-1,down)+fIn(:,out,down+1));

% % corners scheme 4
% fIn(:,in,up) = 1/(4+sqrt(2))*(2*fIn(:,in,up-1)+2*fIn(:,in+1,up)+sqrt(2)*fIn(:,in+1,up-1));
% 
% fIn(:,in,down) = 1/(4+sqrt(2))*(2*fIn(:,in,down+1)+2*fIn(:,in+1,down)+sqrt(2)*fIn(:,in+1,down+1));
%                         
% fIn(:,out,up) = 1/(4+sqrt(2))*(2*fIn(:,out,up-1)+2*fIn(:,out-1,up)+sqrt(2)*fIn(:,out-1,up-1));
%                         
% fIn(:,out,down) = 1/(4+sqrt(2))*(2*fIn(:,out,down+1)+2*fIn(:,out-1,down)+sqrt(2)*fIn(:,out-1,down+1));
                        