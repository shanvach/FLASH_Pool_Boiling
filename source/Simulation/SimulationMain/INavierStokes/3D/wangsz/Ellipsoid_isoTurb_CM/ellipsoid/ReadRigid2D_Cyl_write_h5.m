% Script to write Rigid Body h5 file, to be used by Splash.
% Name of file for simulation: sm_body.1.h5
% The Rigid is 3D.
% -------------------------------------------------------------------------
close all
clear all
clc

%% Define variables:
IAXIS = 1;
JAXIS = 2;
KAXIS = 3;
ITHETA= 3;

NDIM  = 3; % Dimensions

SM_FALSE = 0;
SM_TRUE  = 1;

% Rigid Body constants:
BODYTYPE_RIGID  = 1;

RB_ANNSPHERE = 55;
RB_ANNDISC   = 56;
RB_ANNRBC    = 57;

RB_EULER321 = 1321;

% Kinematics Constants:       
SM_PK_FIXED    =   0;
SM_PK_HARMONIC = 102;
SM_PK_CONSTVEL = 103;

DOFS_per_node   =  9;  % Degrees of freedom per node: 
                       % x,y,ang3,wz

ix = 1; ex = 3;
ia = 4; ea = 6;
iw = 7; ew = 9;

% Free body flag: Set to 1 for all free dofs:
FREE_FLAG=1;

%% Define variables:
NumBods = 1; 

np_el   = 3; % Points per wet surface element: 2 point segment in 2D.

% Output File names:
basedir_out = ['./'];

writefile = SM_TRUE; %SM_FALSE; 

% Prescribed Kinematics Per Body: 
%kinemflag1=['rt_x'];
%kinemflag1=['blck'];
%kinemflag1=['elst']; ELSTON_CASE='a1';
kinemflag1=['free'];

%% Ellipsoid Mesh File:
nX = 3;   % Number of Ellipsoids in the x direction
nY = 3;   % Number of Ellipsoids in the y direction
nZ = 3;   % Number of Ellipsoids in the z direction

Da  =      1; % The major diameter of ellipsoid
Db  =    0.5;
Dc  =    0.5;
Ra  =    Da/2; % The major diameter of ellipsoid
Rb  =    Db/2;
Rc  =    Dc/2;
L  =      1; % Unit depth.

NumBods = nX*nY*nZ;

% Build kinemflag:
for ibd=1:NumBods
  for ilet=1:4
    kinemflag(ibd,ilet) = kinemflag1(ilet);
  end
end

% Read the points from the file
npt = 2429;
nel = 4854;

load 'ellipsoil-025.points';
pos = ellipsoil_025;
load 'ellipsoil-025.elements';
ele = ellipsoil_025;

nnodes = npt + 1;
XYZ = zeros(nnodes,3);
ws_IEN = zeros(nel,3);

% nodes
inod = 1;
for it = 1:npt
    inod = inod + 1;
    XYZ(inod,1) = pos(it,2);
    XYZ(inod,2) = pos(it,3);
    XYZ(inod,3) = pos(it,4);
end

% Elems
iel   = 0;
for it=1:nel
   iel = iel+1;  
   % Segment  
   ws_IEN(iel,1) = ele(it,2)+1; 
   ws_IEN(iel,2) = ele(it,3)+1; 
   ws_IEN(iel,3) = ele(it,4)+1; 
end

% Plot figure:
figure
hold on
trimesh(ws_IEN, XYZ(:,1), XYZ(:,2), XYZ(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
axis equal


figure; hold on
% These params are in case we want to write a set of ellipsoids:
D1 = 2*Da; % Ellipsoids center to center distance in x
D2 = 2*Da; % Ellipsoids center to center distance in y
D3 = 2*Da; % Ellipsoids center to center distance in z

nEllp(IAXIS) = nX;
nEllp(JAXIS) = nY;
nEllp(KAXIS) = nZ;

xstart = -2;
ystart = -2;
zstart = -2;

delx = D1;
dely = D2;
delz = D3;
    
    
ibd = 0;
for kellp=1:nEllp(KAXIS)
zpos = zstart + (kellp-1)*delz; 
XYZ(1,3) = zpos;

for jellp=1:nEllp(JAXIS)
 
   ypos = ystart + (jellp-1)*dely; 
   
   XYZ(1,2) = ypos;
      
   for iellp = 1:nEllp(IAXIS)

    xpos = xstart + (iellp-1)*delx;
    
    XYZ(1,1) = xpos;
    
    ibd = ibd + 1;

    %% Rigid body properties:
     grav_vec     = [0 0 0];
     grav_flag    =  0;
    
    
    %% Inertia properties:
    massDensity = 1.1; % No mass, all prescribed problem.
    volume      = pi*4/3*Ra*Rb*Rc;
    mass        = (massDensity)*volume;
    Ixx         = 1/5*mass*(Rb*Rb + Rc*Rc);
    Iyy         = 1/5*mass*(Ra*Ra + Rc*Rc);
    Izz         = 1/5*mass*(Ra*Ra + Rb*Rb);
    I_body      = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
              
    kx          = 0.0;
    ky          = 0.0;
    kz          = 0.0;
    stiff  = [kx ky kz];
    
    % Analytical body type:
    flag_forceinside = SM_FALSE; %SM_TRUE; % Force inside.
    annbody_type     = RB_ANNDISC;
    annbody_nparam   = 3;
    annbody_param    = [Da]; % Cyl direction is always KAXIS. L = 1.

    %% Restraints:
    % Restrained Surface:
    nrsurf      =  1;
    fix_nodes_A = [2:nnodes]; % Restrained nodes in global numbering
    nfix = length(fix_nodes_A);
    
    %%%% UP TO HERE !!!
    if(FREE_FLAG)
    % No Restrained DOFS of Node 1:
    restnode      = 1;
    nrestcoords   = 0; % No restraints.
    maxrestparams = 1;
    restnodes     = 0;
    restdofs      = 0;
    resttype      = 0; 
    nparam        = 0; 
    param         = 0; 
    
    else
       
    % Change restraint on z displacement 
    break; % not set for constraint parameter currently
    end
    
    %% Export Solid Model to file:
    % Values to write:
    TRANSFORMATION  = RB_EULER321;  % Use 2D Rotation matrix
    nnp             =      nnodes;  % Total number of nodes
    if (ibd ==1)
    nel_load        =         nel;  % nel surface elements
    nel             =           1;  % 1 rigid body element
    end
    ned             =     NDIM;
    max_eltype      =       15;  % Max nodes per element type (that is one node elem)
    max_eltype_load =        2;  % Max nodes per surface element (02 is triangle)
    eltype          =       15*ones(1,nel); % Element type : One node Rigid Body
    eltype_load     =        2*ones(1,nel_load); % Surface Element type : triangle
    kinematics_idx   =       1;  %
    
    x = XYZ(:,IAXIS); % Positions of all nodes : nnodes
    y = XYZ(:,JAXIS);
    z = XYZ(:,KAXIS); % Positions of all nodes : nnodes
    
    IEN = [1];  % One rigid body with node 1 as connectivity
    
    
    if (writefile == SM_TRUE)
    
    % Write hdf5 file:
    % Mesh:
    hfilename = ['sm_body.' num2str(ibd,'%5.5d') '.h5'];
    hfile=[basedir_out hfilename];
    fprintf(1,'     Body %d hdf5 file started...\n',ibd);
    hdf5write(hfile,'BodyType',int32(BODYTYPE_RIGID))
    hdf5write(hfile,'mesh/x',x,'WriteMode','append')
    hdf5write(hfile,'mesh/y',y,'WriteMode','append')
    hdf5write(hfile,'mesh/z',z,'WriteMode','append')
    hdf5write(hfile,'mesh/nnp',int32(nnp),'WriteMode','append')
    hdf5write(hfile,'mesh/ned',int32(ned),'WriteMode','append')
    
    % Body:
    hdf5write(hfile,'body/IEN',int32(IEN'),'WriteMode','append')
    hdf5write(hfile,'body/nel',int32(nel),'WriteMode','append')
    %hdf5write(hfile,'body/nee',int32(nee),'WriteMode','append')
    hdf5write(hfile,'body/max_eltype',int32(max_eltype),'WriteMode','append')
    hdf5write(hfile,'body/DOFS_per_node',int32(DOFS_per_node),'WriteMode','append')
    
    % Body Properties:
    hdf5write(hfile,'body/eltype',int32(eltype),'WriteMode','append')
    hdf5write(hfile,'body/Mass',double(mass),'WriteMode','append')
    hdf5write(hfile,'body/Volume',double(volume),'WriteMode','append')
    hdf5write(hfile,'body/I_body',double(I_body),'WriteMode','append')
    hdf5write(hfile,'body/trmatrix',int32(TRANSFORMATION),'WriteMode','append')
    hdf5write(hfile,'body/Stiffness',double(stiff),'WriteMode','append')
    hdf5write(hfile,'body/gravity',double(grav_vec),'WriteMode','append')
    hdf5write(hfile,'body/gravity_flag',int32(grav_flag),'WriteMode','append')
    hdf5write(hfile,'body/flag_forceinside',int32(flag_forceinside),'WriteMode','append')
    hdf5write(hfile,'body/annbody_type',int32(annbody_type),'WriteMode','append')    
    hdf5write(hfile,'body/annbody_nparam',int32(annbody_nparam),'WriteMode','append')    
    hdf5write(hfile,'body/annbody_param',double(annbody_param),'WriteMode','append')
    
    % WetSurface:
    hdf5write(hfile,'WetSurface/nel',int32(nel_load),'WriteMode','append')
    hdf5write(hfile,'WetSurface/max_eltype',int32(max_eltype_load),'WriteMode','append')
    hdf5write(hfile,'WetSurface/IEN',int32(ws_IEN'),'WriteMode','append')
    hdf5write(hfile,'WetSurface/eltype',int32(eltype_load),'WriteMode','append')
    
    % RestSurface:
    if( isempty(fix_nodes_A) )
        hdf5write(hfile,'RestSurface/fix_list',int32(0),'WriteMode','append')
    else
        hdf5write(hfile,'RestSurface/fix_list',int32(fix_nodes_A),'WriteMode','append')
    end
    hdf5write(hfile,'RestSurface/nrestsurf',int32(nrsurf),'WriteMode','append')
    hdf5write(hfile,'RestSurface/nfix',int32(nfix),'WriteMode','append')
    hdf5write(hfile,'RestSurface/kinematics_idx',int32(kinematics_idx),'WriteMode','append')
    
    
    % Write Restrained nodes, just node 1:
    hdf5write(hfile,'RestNodes/nrestcoords'  ,int32(nrestcoords),'WriteMode','append')
    hdf5write(hfile,'RestNodes/maxrestparams',int32(maxrestparams),'WriteMode','append')
    hdf5write(hfile,'RestNodes/nodes'        ,int32(restnodes),'WriteMode','append')
    hdf5write(hfile,'RestNodes/restdof'      ,int32(restdofs),'WriteMode','append')
    hdf5write(hfile,'RestNodes/restype'      ,int32(resttype),'WriteMode','append')
    hdf5write(hfile,'RestNodes/nparam'       ,int32(nparam),'WriteMode','append')
    hdf5write(hfile,'RestNodes/param'        ,double(param),'WriteMode','append')
    
    fprintf(1,'   Body %d hdf5 file done.\n',ibd);

    
    end
    
    trimesh(ws_IEN, XYZ(1:nnodes,1)+xpos, XYZ(1:nnodes,2)+ypos, XYZ(1:nnodes,3)+zpos)

    XYZ(1,:)
    
  end
end
end

xlabel('x')
ylabel('y')
axis equal

return
