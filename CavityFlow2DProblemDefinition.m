function [Dimensionless, BCs,ICs, Char] = CavityFlow2DProblemDefinition
%%%%Problem can be difined both in physical or dimensionless system %%%%%%
%% Problem Definition in Physical System

% Physical Domain and Time Interval
% The domain of Cavity flow is a recktangle in 2 dimenstion looks like 
%                Up
%       -------------------
%       |                 |
% Left  |                 |  Right
%       |                 |
%       |-----------------|
%              Down
% The physical system is in CGS(cm-s-g) system
Physical.Lx      = 1;   % length in x-direction  [cm]
Physical.Ly      = 1;   % length in y-direction  [cm]
Physical.TimeMax = 10;  % stopping time [s]

Physical.Rho     = 1;   % mass density [g/cm^3] 
Physical.UMax    = 1;   % maximum velocity [cm/t]
Physical.Nu      = 0.1; % shear viscosity [g/(cm*s)] 

% Characteristic Constants and Reynolds Number 
Char.PhysicalRho        = Physical.Rho;
Char.PhysicalUMax       = Physical.UMax;
Char.PhysicalLength     = Physical.Ly;
Char.PhysicalTime       = Physical.Ly/Physical.UMax;

% Physical Boundary Conditions and Initial Conditions 
% this part can be defined by function file, for velocity and pressure on
% the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCs.UxLeft    = inline('Physical.UMax*zeros(size(y))','x','y','Physical.UMax');
BCs.UyLeft    = inline('Physical.UMax*zeros(size(y))','x','y','Physical.UMax');
BCs.UxRight   = inline('Physical.UMax*zeros(size(y))','x','y','Physical.UMax');
BCs.UyRight   = inline('Physical.UMax*zeros(size(y))','x','y','Physical.UMax');
BCs.UxUp      = inline('Physical.UMax*ones(size(x))','x','y','Physical.UMax');
BCs.UyUp      = inline('Physical.UMax*zeros(size(x))','x','y','Physical.UMax');
BCs.UxDown    = inline('Physical.UMax*zeros(size(x))','x','y','Physical.UMax');
BCs.UyDown    = inline('Physical.UMax*zeros(size(x))','x','y','Physical.UMax');

% type of boundary condition
BCs.Left      = 'velocity';
BCs.Right     = 'velocity';
BCs.Up        = 'velocity';
BCs.Down      = 'velocity';

ICs.Ux        = inline('Physical.UMax*zeros(size(x))','x','y','Physical.UMax');
ICs.Uy        = inline('Physical.UMax*zeros(size(x))','x','y','Physical.UMax');
%% Problem Definition in Dimensionless System
% This part is generated automatically provided everything is given above
% otherwise, one can specify the problem directly starting from here 

% Dimensionless Domain and Time Interval
Dimensionless.Lx     = Physical.Lx/Char.PhysicalLength;
Dimensionless.Ly     = Physical.Ly/Char.PhysicalLength;   
Dimensionless.TimeMax= Physical.TimeMax/Char.PhysicalTime;

Dimensionless.Rho    = Physical.Rho/Char.PhysicalRho;
Dimensionless.UMax   = Physical.UMax/Char.PhysicalUMax;
Dimensionless.Nu     = Physical.Nu*Char.PhysicalTime^2/Char.PhysicalLength;
Dimensionless.Re     = 1/Dimensionless.Nu;


Dimensionless.Fx     = inline('1e-2*Dimensionless.UMax*t*zeros(size(x))','x','y','t','Dimensionless.UMax');
Dimensionless.Fy     = inline('-1e-2*Dimensionless.UMax*t*zeros(size(y))','x','y','t','Dimensionless.UMax');

% Dimensionless Boundary Conditions and Initial Conditions 
% this part can be defined either as constants or functions for velocity 
% and pressure on the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCs.UxLeft = inline('Dimensionless.UMax*zeros(size(y))','x','y','Dimensionless.UMax');
BCs.UyLeft = inline('Dimensionless.UMax*zeros(size(y))','x','y','Dimensionless.UMax');
BCs.UxRight= inline('Dimensionless.UMax*zeros(size(y))','x','y','Dimensionless.UMax');
BCs.UyRight= inline('Dimensionless.UMax*zeros(size(y))','x','y','Dimensionless.UMax');
BCs.UxUp   = inline('Dimensionless.UMax*ones(size(x))','x','y','Dimensionless.UMax');
BCs.UyUp   = inline('Dimensionless.UMax*zeros(size(x))','x','y','Dimensionless.UMax');
BCs.UxDown = inline('Dimensionless.UMax*zeros(size(x))','x','y','Dimensionless.UMax');
BCs.UyDown = inline('Dimensionless.UMax*zeros(size(x))','x','y','Dimensionless.UMax');

% type of boundary condition
BCs.Left      = 'velocity';
BCs.Right     = 'velocity';
BCs.Up        = 'velocity';
BCs.Down      = 'velocity';


ICs.Ux        = inline('Dimensionless.UMax*zeros(size(x))','x','y','Dimensionless.UMax');
ICs.Uy        = inline('Dimensionless.UMax*zeros(size(x))','x','y','Dimensionless.UMax');
ICs.dPdx      = 0;
