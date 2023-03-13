% Simulation of a simple EMI gasketed enclosure made of a clamshell metal box

close all
clear
clc
confirm_recursive_rmdir(0)

%%%%%%%  openEMS Setup %%%%%%%
physical_constants;
unit = 1e-3;        % all length in mm
f0 = 1e9;           % center frequency
fc = 2e9;           % 20 dB corner frequency

% plot setup
plot_points = 1001; % number of points
fl = 1e6;           % lowest freq
fh = f0+fc;         % highest freq

% air properties
air.epsR   = 1;      % relative permeability
air.kappa  = 1e-9;   % air conductivity
air.mue  = 1;

% EMI gasket properties
gasket.epsR   = 4;        % relative permeability (neoprene)
gasket.kappa  = 3e2;      % gasket conductivity (carbon)
gasket.width  = 200;      % mm
gasket.length = 200;      % mm
gasket.z_thickness = 0.5; % mm
gasket.xy_thickness = 3;  % mm
gasket.cells = 4;       % number of FDTD cells to capture gasket thickness

% metal enclosure properties
enc.width  = 200;         % mm
enc.length = 200;         % mm
enc.height = 60;          % mm
enc.wall_thickness = 10;  % mm

% ground bond wire properties
gbw.thickness = 1;        % mm
gbw.x = 0;                % mm, dead center
gbw.y = 0;                % mm, dead center

% soure feeding
feed.R = 50;     %feed resistance, ohms

% size of the simulation box
SimBox = [enc.width*4 enc.length*4 enc.height*6];

%% Setup FDTD Parameter & Excitation Function

FDTD = InitFDTD( 'NrTs', 60000 );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/3 SimBox(3)*2/3];

%%%%%%% Drawing of EM elements with materials %%%%%%%
CSX = AddMaterial( CSX, 'air' );
CSX = SetMaterialProperty( CSX, 'air', 'Epsilon', air.epsR, 'Mue', air.mue , 'Kappa', air.kappa );

% Create top cover by building a solid block then applying an air block
% Set the top cover on TOP of the gasket, where the gasket z base = 0
CSX = AddMetal( CSX, 'top_cover' ); % create a perfect electric conductor (PEC)
start = [-enc.width/2 -enc.length/2 gasket.z_thickness];
stop  = [ enc.width/2  enc.length/2 enc.height];
CSX = AddBox(CSX, 'top_cover', 1, start, stop); % add a box-primitive to the metal top cover
start = [-enc.width/2+enc.wall_thickness -enc.length/2+enc.wall_thickness gasket.z_thickness];
stop  = [ enc.width/2-enc.wall_thickness enc.length/2-enc.wall_thickness enc.height-enc.wall_thickness];
CSX = AddBox(CSX, 'air', 2, start, stop); % subtract out the metal by applying higher priority 'air'

% Create Conductive EMI Gasket
% Gasket z base = 0
CSX = AddMaterial( CSX, 'gasket' );
CSX = SetMaterialProperty( CSX, 'gasket', 'Epsilon', gasket.epsR, 'Kappa', gasket.kappa );
start = [-gasket.width/2 -gasket.length/2 0];
stop  = [gasket.width/2  gasket.length/2 gasket.z_thickness];
CSX = AddBox( CSX, 'gasket', 1, start, stop ); % add a box-primitive to the conductive sheet of EMI gasket
start = [-gasket.width/2+gasket.xy_thickness -gasket.length/2+gasket.xy_thickness 0];
stop  = [gasket.width/2-gasket.xy_thickness  gasket.length/2-gasket.xy_thickness gasket.z_thickness];
CSX = AddBox( CSX, 'air', 2, start, stop ); % subtract out the conductive material by applying higher priority 'air'
% add extra cells to discretize the gasket thickness
mesh.z = [linspace(0, gasket.z_thickness, gasket.cells + 1) mesh.z];

% Create bottom cover by building a solid block then applying an air block
% Set the bottom cover below the gasket, where the gasket z base = 0
CSX = AddMetal( CSX, 'bottom_cover' ); % create a perfect electric conductor (PEC)
start = [-enc.width/2 -enc.length/2 0];
stop  = [ enc.width/2  enc.length/2 -enc.height];
CSX = AddBox(CSX, 'bottom_cover', 1, start, stop); % add a box-primitive to the metal top cover
start = [-enc.width/2+enc.wall_thickness -enc.length/2+enc.wall_thickness gasket.z_thickness];
stop  = [ enc.width/2-enc.wall_thickness enc.length/2-enc.wall_thickness -enc.height+enc.wall_thickness];
CSX = AddBox(CSX, 'air', 2, start, stop); % subtract out the metal by applying higher priority 'air'

% Create Ground bond wire that connects both shells
# don't forget to get the priority higher than the air boxes
CSX = AddMetal( CSX, 'gnd_bond_wire' ); % create a perfect electric conductor (PEC)
start = [gbw.x-gbw.thickness/2 gbw.y-gbw.thickness/2 enc.height-enc.wall_thickness];
stop  = [ gbw.x+gbw.thickness/2 gbw.y+gbw.thickness/2 -enc.height+enc.wall_thickness];
CSX = AddBox(CSX, 'gnd_bond_wire', 3, start, stop);

% Apply the Excitation & Resist as a Current Source
% Connect the TOP and BOTTOM covers through a vertical stub at the corner
% of the gasket
start = [-gasket.width/2+1 -gasket.length/2+1 0];
stop  = [-gasket.width/2+1 -gasket.length/2+1 gasket.z_thickness];
[CSX port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], true);

% Finalize the Mesh
% detect all edges except of the patch
mesh = DetectEdges(CSX, mesh);
% generate a smooth mesh with max. cell size: lambda_min / 20
mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20);
CSX = DefineRectGrid(CSX, unit, mesh);

%%%%%%% Prepare and Run Simulation %%%%%%%
Sim_Path = 'tmp_EMI_gasket';
Sim_CSX = 'EMI_gasket.xml';
% create an empty working directory
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%%%%%%% Postprocessing & Plots %%%%%%%

freq = linspace(fl, fh, plot_points);
port = calcPort(port, Sim_Path, freq);

% plot feed point impedance
Zin = port.uf.tot ./ port.if.tot;
figure
loglog( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
loglog( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
s11 = port.uf.ref ./ port.uf.inc;
figure
loglog( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow
