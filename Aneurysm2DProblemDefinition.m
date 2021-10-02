function [DimensionlessP,DimensionlessQ,DimensionlessC, BCsP, BCsQ,BCsC, ICsP,ICsQ,ICsC,CharP,CharQ,CharC] = Aneurysm2DProblemDefinition
%%%%Problem can be difined both in physical or dimensionless system %%%%%%
%% Problem Definition in Physical System

% Physical Domain and Time Interval
% The domain of PoiseuilleCavity flow is two recktangles in 2 dimenstion look like 
% |----------------UpP---------------------------|
% | -->           -->                        --> |
% |LeftP                                   RightP|
% | -->           -->                        --> |  
% |----------    UpC\DownP    -------------------|
%           |                 |
%     LeftC |                 | RightC
%           |                 |
%           |-----------------|
%                 DownC


% The physical system is in CGS(cm-s-g) system
PhysicalP.Lx      = 3;   % length in x-direction  [cm]
PhysicalP.Ly      = 1;   % length in y-direction  [cm]
PhysicalP.TimeMax = 1000;  % stopping time [s]

PhysicalP.Rho     = 1;   % mass density [g/cm^3] 
PhysicalP.UMax    = 1;   % maximum velocity [cm/t]
PhysicalP.Nu      = 0.1; % shear viscosity [g/(cm*s)] 

% Characteristic Constants and Reynolds Number 
CharP.PhysicalRho        = PhysicalP.Rho;
CharP.PhysicalUMax       = PhysicalP.UMax;
CharP.PhysicalLength     = PhysicalP.Ly;
CharP.PhysicalTime       = PhysicalP.Ly/PhysicalP.UMax;

% Physical Boundary Conditions and Initial Conditions 
% this part can be defined by function file, for velocity and pressure on
% the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCsP.UxLeft    = inline('4*Physical.UMax*(y-y.*y)','x','y','Physical.UMax');
BCsP.UyLeft    = inline('PhysicalP.UMax*zeros(size(y))','x','y','PhysicalP.UMax');
BCsP.UxRight   = inline('PhysicalP.UMax*zeros(size(y))','x','y','PhysicalP.UMax');
BCsP.UyRight   = inline('PhysicalP.UMax*zeros(size(y))','x','y','PhysicalP.UMax');
BCsP.UxUp      = inline('PhysicalP.UMax*zeros(size(x))','x','y','PhysicalP.UMax');
BCsP.UyUp      = inline('PhysicalP.UMax*zeros(size(x))','x','y','PhysicalP.UMax');
BCsP.UxDown    = inline('PhysicalP.UMax*zeros(size(x))','x','y','PhysicalP.UMax');
BCsP.UyDown    = inline('PhysicalP.UMax*zeros(size(x))','x','y','PhysicalP.UMax');

% type of boundary condition
BCsP.Left      = 'velocity';
BCsP.Right     = 'pressure';
BCsP.Up        = 'velocity';
BCsP.Down      = 'velocity';
BCsP.Curve     = 'velocity';

ICsP.Ux        = inline('4*PhysicalP.UMax*(y-y.*y)','x','y','PhysicalP.UMax');
ICsP.Uy        = inline('PhysicalP.UMax*zeros(size(x))','x','y','PhysicalP.UMax');





% The physical system is in CGS(cm-s-g) system
PhysicalQ.Lx      = 1;   % length in x-direction  [cm]
PhysicalQ.Ly      = 3;   % length in y-direction  [cm]
PhysicalQ.TimeMax = 10;  % stopping time [s]

PhysicalQ.Rho     = 1;   % mass density [g/cm^3] 
PhysicalQ.UMax    = 1;   % maximum velocity [cm/t]
PhysicalQ.Nu      = 0.1; % shear viscosity [g/(cm*s)] 

% Characteristic Constants and Reynolds Number 
CharQ.PhysicalRho        = PhysicalQ.Rho;
CharQ.PhysicalUMax       = PhysicalQ.UMax;
CharQ.PhysicalLength     = PhysicalQ.Lx;
CharQ.PhysicalTime       = PhysicalQ.Ly/PhysicalQ.UMax;

% Physical Boundary Conditions and Initial Conditions 
% this part can be defined by function file, for velocity and pressure on
% the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCsQ.UxLeft    = inline('PhysicalQ.UMax*zeros(size(y))','x','y','PhysicalQ.UMax');
BCsQ.UyLeft    = inline('PhysicalQ.UMax*zeros(size(y))','x','y','PhysicalQ.UMax');
BCsQ.UxRight   = inline('PhysicalQ.UMax*zeros(size(y))','x','y','PhysicalQ.UMax');
BCsQ.UyRight   = inline('PhysicalQ.UMax*zeros(size(y))','x','y','PhysicalQ.UMax');
BCsQ.UxUp      = inline('PhysicalQ.UMax*zeros(size(x))','x','y','PhysicalQ.UMax');
BCsQ.UyUp      = inline('PhysicalQ.UMax*zeros(size(x))','x','y','PhysicalQ.UMax');
BCsQ.UxDown    = inline('PhysicalQ.UMax*zeros(size(x))','x','y','PhysicalQ.UMax');
BCsQ.UyDown    = inline('PhysicalQ.UMax*zeros(size(x))','x','y','PhysicalQ.UMax');

% type of boundary condition
BCsQ.Left      = 'velocity';
BCsQ.Right     = 'velocity';
BCsQ.Up        = 'pressure';
BCsQ.Down      = 'velocity';
BCsQ.Curve     = 'velocity';

ICsQ.Ux        = inline('PhysicalQ.UMax*zeros(size(x))','x','y','PhysicalQ.UMax');
ICsQ.Uy        = inline('PhysicalQ.UMax*zeros(size(x))','x','y','PhysicalQ.UMax');






% The physical system is in CGS(cm-s-g) system
PhysicalC.Lx      = 1;   % length in x-direction  [cm]
PhysicalC.Ly      = 1;   % length in y-direction  [cm]
PhysicalC.TimeMax = 10;  % stopping time [s]

PhysicalC.Rho     = 1;   % mass density [g/cm^3] 
PhysicalC.UMax    = 1;   % maximum velocity [cm/t]
PhysicalC.Nu      = 0.1; % shear viscosity [g/(cm*s)] 

% Characteristic Constants and Reynolds Number 
CharC.PhysicalRho        = PhysicalC.Rho;
CharC.PhysicalUMax       = PhysicalC.UMax;
CharC.PhysicalLength     = PhysicalC.Ly;
CharC.PhysicalTime       = PhysicalC.Ly/PhysicalC.UMax;

% Physical Boundary Conditions and Initial Conditions 
% this part can be defined by function file, for velocity and pressure on
% the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCsC.UxLeft    = inline('PhysicalC.UMax*zeros(size(y))','x','y','PhysicalC.UMax');
BCsC.UyLeft    = inline('PhysicalC.UMax*zeros(size(y))','x','y','PhysicalC.UMax');
BCsC.UxRight   = inline('PhysicalC.UMax*zeros(size(y))','x','y','PhysicalC.UMax');
BCsC.UyRight   = inline('PhysicalC.UMax*zeros(size(y))','x','y','PhysicalC.UMax');
BCsC.UxUp      = inline('PhysicalC.UMax*zeros(size(x))','x','y','PhysicalC.UMax');
BCsC.UyUp      = inline('PhysicalC.UMax*zeros(size(x))','x','y','PhysicalC.UMax');
BCsC.UxDown    = inline('PhysicalC.UMax*zeros(size(x))','x','y','PhysicalC.UMax');
BCsC.UyDown    = inline('PhysicalC.UMax*zeros(size(x))','x','y','PhysicalC.UMax');
% BCsC.PresRight = inline('zeros(size(y))','x','y','Physical.UMax');

% type of boundary condition
BCsC.Left      = 'velocity';
BCsC.Right     = 'velocity';
BCsC.Up        = 'velocity';
BCsC.Down      = 'velocity';
BCsC.Curve     = 'velocity';


ICsC.Ux        = inline('PhysicalC.UMax*zeros(size(x))','x','y','PhysicalC.UMax');
ICsC.Uy        = inline('PhysicalC.UMax*zeros(size(x))','x','y','PhysicalC.UMax');









%% Problem Definition in Dimensionless System
% This part is generated automatically provided everything is given above
% otherwise, one can specify the problem directly starting from here 

% Dimensionless Domain and Time Interval
DimensionlessP.Lx     = PhysicalP.Lx/CharP.PhysicalLength;
DimensionlessP.Ly     = PhysicalP.Ly/CharP.PhysicalLength;   
DimensionlessP.TimeMax= PhysicalP.TimeMax/CharP.PhysicalTime;

DimensionlessP.Rho    = PhysicalP.Rho/CharP.PhysicalRho;
DimensionlessP.UMax   = PhysicalP.UMax/CharP.PhysicalUMax;
DimensionlessP.Nu     = PhysicalP.Nu*CharP.PhysicalTime^2/CharP.PhysicalLength;
DimensionlessP.Re     = 1/DimensionlessP.Nu;


DimensionlessP.Fx     = inline('1e-2*DimensionlessP.UMax*t*zeros(size(x))','x','y','t','DimensionlessP.UMax');
DimensionlessP.Fy     = inline('-1e-2*DimensionlessP.UMax*t*zeros(size(y))','x','y','t','DimensionlessP.UMax');

% Dimensionless Boundary Conditions and Initial Conditions 
% this part can be defined either as constants or functions for velocity 
% and pressure on the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCsP.UxLeft = inline('4*DimensionlessP.UMax*(y-y.*y)','x','y','DimensionlessP.UMax');
BCsP.UyLeft = inline('DimensionlessP.UMax*zeros(size(y))','x','y','DimensionlessP.UMax');
BCsP.UxRight= inline('DimensionlessP.UMax*zeros(size(y))','x','y','DimensionlessP.UMax');
BCsP.UyRight= inline('DimensionlessP.UMax*zeros(size(y))','x','y','DimensionlessP.UMax');
BCsP.UxUp   = inline('DimensionlessP.UMax*zeros(size(x))','x','y','DimensionlessP.UMax');
BCsP.UyUp   = inline('DimensionlessP.UMax*zeros(size(x))','x','y','DimensionlessP.UMax');
BCsP.UxDown = inline('DimensionlessP.UMax*zeros(size(x))','x','y','DimensionlessP.UMax');
BCsP.UyDown = inline('DimensionlessP.UMax*zeros(size(x))','x','y','DimensionlessP.UMax');

BCsP.RhoRight  = inline('ones(size(y))','x','y','DimensionlessP.UMax');

% type of boundary condition
BCsP.Left      = 'velocity';
BCsP.Right     = 'pressure';
BCsP.Up        = 'velocity';
BCsP.Down      = 'velocity';


ICsP.Ux        = inline('4*DimensionlessP.UMax*(y-y.*y)','x','y','DimensionlessP.UMax');
ICsP.Uy        = inline('DimensionlessP.UMax*zeros(size(x))','x','y','DimensionlessP.UMax');
ICsP.dPdx      = -8/DimensionlessP.Re;






% Dimensionless Domain and Time Interval
DimensionlessQ.Lx     = PhysicalQ.Lx/CharQ.PhysicalLength;
DimensionlessQ.Ly     = PhysicalQ.Ly/CharQ.PhysicalLength;   
DimensionlessQ.TimeMax= PhysicalQ.TimeMax/CharQ.PhysicalTime;

DimensionlessQ.Rho    = PhysicalQ.Rho/CharQ.PhysicalRho;
DimensionlessQ.UMax   = PhysicalQ.UMax/CharQ.PhysicalUMax;
DimensionlessQ.Nu     = PhysicalQ.Nu*CharQ.PhysicalTime^2/CharQ.PhysicalLength;
DimensionlessQ.Re     = 1/DimensionlessQ.Nu;


DimensionlessQ.Fx     = inline('1e-2*DimensionlessQ.UMax*t*zeros(size(x))','x','y','t','DimensionlessQ.UMax');
DimensionlessQ.Fy     = inline('-1e-2*DimensionlessQ.UMax*t*zeros(size(y))','x','y','t','DimensionlessQ.UMax');

% Dimensionless Boundary Conditions and Initial Conditions 
% this part can be defined either as constants or functions for velocity 
% and pressure on the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCsQ.UxLeft = inline('DimensionlessQ.UMax*zeros(size(y))','x','y','DimensionlessQ.UMax');
BCsQ.UyLeft = inline('DimensionlessQ.UMax*zeros(size(y))','x','y','DimensionlessQ.UMax');
BCsQ.UxRight= inline('DimensionlessQ.UMax*zeros(size(y))','x','y','DimensionlessQ.UMax');
BCsQ.UyRight= inline('DimensionlessQ.UMax*zeros(size(y))','x','y','DimensionlessQ.UMax');
BCsQ.UxUp   = inline('DimensionlessQ.UMax*zeros(size(x))','x','y','DimensionlessQ.UMax');
BCsQ.UyUp   = inline('DimensionlessQ.UMax*zeros(size(x))','x','y','DimensionlessQ.UMax');
BCsQ.UxDown = inline('DimensionlessQ.UMax*zeros(size(x))','x','y','DimensionlessQ.UMax');
BCsQ.UyDown = inline('DimensionlessQ.UMax*zeros(size(x))','x','y','DimensionlessQ.UMax');

BCsQ.RhoUp  = inline('ones(size(x))','x','y','DimensionlessQ.UMax');

% type of boundary condition
BCsQ.Left      = 'velocity';
BCsQ.Right     = 'velocity';
BCsQ.Up        = 'pressure';
BCsQ.Down      = 'velocity';


ICsQ.Ux        = inline('DimensionlessQ.UMax*zeros(size(x))','x','y','DimensionlessQ.UMax');
ICsQ.Uy        = inline('DimensionlessQ.UMax*zeros(size(x))','x','y','DimensionlessQ.UMax');
ICsQ.dPdx      = 0;









% Dimensionless Domain and Time Interval
DimensionlessC.Lx     = PhysicalC.Lx/CharC.PhysicalLength;
DimensionlessC.Ly     = PhysicalC.Ly/CharC.PhysicalLength;   
DimensionlessC.TimeMax= PhysicalC.TimeMax/CharC.PhysicalTime;

DimensionlessC.Rho    = PhysicalC.Rho/CharC.PhysicalRho;
DimensionlessC.UMax   = PhysicalC.UMax/CharC.PhysicalUMax;
DimensionlessC.Nu     = PhysicalC.Nu*CharC.PhysicalTime^2/CharC.PhysicalLength;
DimensionlessC.Re     = 1/DimensionlessC.Nu;

DimensionlessC.Fx     = inline('1e-3*DimensionlessC.UMax*t*zeros(size(x))','x','y','t','DimensionlessC.UMax');
DimensionlessC.Fy     = inline('-1e-3*DimensionlessC.UMax*t*zeros(size(y))','x','y','t','DimensionlessC.UMax');


% DimensionlessC Boundary Conditions and Initial Conditions 
% this part can be defined either as constants or functions for velocity 
% and pressure on the boundary Left-Right, Up-Down, Front-Back(for 3D)
BCsC.UxCurve = inline('DimensionlessC.UMax*zeros(size(y))','x','y','DimensionlessC.UMax');
BCsC.UyCurve = inline('DimensionlessC.UMax*zeros(size(y))','x','y','DimensionlessC.UMax');



% type of boundary condition
BCsC.Curve      = 'velocity';

ICsC.Ux        = inline('DimensionlessC.UMax*zeros(size(x))','x','y','DimensionlessC.UMax');
ICsC.Uy        = inline('DimensionlessC.UMax*zeros(size(x))','x','y','DimensionlessC.UMax');
ICsC.dPdx      = 0;



