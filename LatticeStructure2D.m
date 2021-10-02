function LS2D = LatticeStructure2D
%% Lattice Structure D2Q9
%  7  3  6
%   \ | /
% 4 - 1 - 2 
%   / | \
%  8  5  9
LS2D.w   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];% quadrature weights
LS2D.opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7]; % opposite indices
LS2D.cx  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];% lattice velocity -x
LS2D.cy  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];% lattice velocity -y

% for regularized boundary dynamics, also for the non-equilibrium
% distribution related computation for second order moment tensors. 
LS2D.Q = [];
LS2D.Cxy = [];
LS2D.wQ  =[];
for i = 1:9
    LS2D.Cxy = [LS2D.Cxy,reshape([LS2D.cx(i);LS2D.cy(i)]*[LS2D.cx(i),LS2D.cy(i)],4,1)];
    LS2D.Q   = [LS2D.Q;reshape([LS2D.cx(i);LS2D.cy(i)]*[LS2D.cx(i),LS2D.cy(i)]-1/3*eye(2),1,4)];   
    LS2D.wQ  = [LS2D.wQ;LS2D.w(i)*reshape([LS2D.cx(i);LS2D.cy(i)]*[LS2D.cx(i),LS2D.cy(i)]-1/3*eye(2),1,4)]; 
end


