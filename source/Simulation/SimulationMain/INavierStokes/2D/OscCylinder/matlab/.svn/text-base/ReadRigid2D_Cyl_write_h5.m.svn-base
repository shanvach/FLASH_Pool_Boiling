% Script to write Rigid Body h5 file, to be used by Splash.
% Name of file for simulation: sm_body.1.h5
% The Rigid is 2D.
% -------------------------------------------------------------------------
close all
clear all
clc

%% Define variables:
IAXIS = 1;
JAXIS = 2;
KAXIS = 3;
ITHETA= 3;

NDIM  = 2; % Dimensions

SM_FALSE = 0;
SM_TRUE  = 1;

% Rigid Body constants:
BODYTYPE_RIGID  = 1;

RB_ANNSPHERE = 55;
RB_ANNDISC   = 56;
RB_ANNRBC    = 57;

RB_TWODIM = 1002;

% Kinematics Constants:       
SM_PK_FIXED    =   0;
SM_PK_HARMONIC = 102;
SM_PK_CONSTVEL = 103;

DOFS_per_node   =  4;  % Degrees of freedom per node: 
                       % x,y,ang3,wz

ix = 1; ex = 2;
ia = 3; ea = 3;
iw = 4; ew = 4;

% Free body flag: Set to 1 for all free dofs:
FREE_FLAG=0;

%% Define variables:
NumBods = 1; 

np_el   = 2; % Points per wet surface element: 2 point segment in 2D.

% Output File names:
basedir_out = ['/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/source/Simulation/SimulationMain/INavierStokes/2D/OscCylinder/'];

writefile = SM_TRUE; %SM_FALSE; 

% Prescribed Kinematics Per Body: 
%kinemflag1=['rt_x'];
%kinemflag1=['blck'];
kinemflag1=['elst']; ELSTON_CASE='a1';

%% Cylinder Mesh File:
nCyl_1 = 1;   % Number of Cylinders in the x direction
nCyl_2 = 1;   % Number of Cylinders in the y direction

D  =      1; % Cylinder Diameter
L  =      1; % Unit depth.

n_sections_theta = 256; % Number of sections in circumf direction.


NumBods = nCyl_1*nCyl_2;

% Build kinemflag:
for ibd=1:NumBods
  for ilet=1:4
    kinemflag(ibd,ilet) = kinemflag1(ilet);
  end
end

dtheta = pi*D/n_sections_theta;

% Now make the mesh
% Radius, circumference and npoints in cylinder:
r      =                 D/2;
circmf =              2*pi*r;
npt    = ceil(circmf/dtheta);

% Check if number is odd -> make even:
if mod(npt,2) ~= 0
    npt=npt+1;
end

% Delta theta:
dth=2*pi/npt;

% 1D arrays:
theta = linspace(0,2*pi-dth,npt);
X_circ = r*cos(theta);
Y_circ = r*sin(theta);

% Nodes and Elements:
nnodes = npt+1;
nel    = npt;
XY     = zeros(nnodes,2);
ELEM   = zeros(nel,2);

% Nodes
inod = 1;
for it=1:npt
  inod = inod+1;
  XY(inod,:) = [X_circ(it) Y_circ(it)];  
end


% Elems
iel   = 0;
for it=1:npt
      
   iel = iel+1;  
   % Segment  
   ws_IEN(iel,:) = [it+1 it+2]; 
    
end
ws_IEN(end,2) = ws_IEN(1,1);


% Plot figure:
figure
hold on
plot(XY(:,IAXIS),XY(:,JAXIS),'xb')
for iel=1:nel
   plot([XY(ws_IEN(iel,1),IAXIS) XY(ws_IEN(iel,2),IAXIS)], ...
        [XY(ws_IEN(iel,1),JAXIS) XY(ws_IEN(iel,2),JAXIS)],'-k')
end
xlabel('x')
ylabel('y')
axis equal


figure; hold on
% These params are in case we want to write a set of cylinders:
D1 = 1.5*D; % Cylinders center to center distance in x
D2 = 1.5*D; % Cylinders center to center distance in y

nCyl(IAXIS) = nCyl_1;
nCyl(JAXIS) = nCyl_2;

xstart = -(nCyl_1-1)/2*D1;
ystart = -(nCyl_2-1)/2*D2;

delx = D1;
dely = D2;
    
    
ibd = 0;
for jcyl=1:nCyl(JAXIS)
 
   ypos = ystart + (jcyl-1)*dely;   
      
   for icyl = 1:nCyl(IAXIS)

    xpos = xstart + (icyl-1)*delx;
       
    ibd = ibd + 1;

    %% Rigid body properties:
    grav_vec     = [0 -1.];
    grav_flag    =       0; % 0 = No gravity, 1 with gravity
    
    
    %% Inertia properties:
    massDensity = 1.; % No mass, all prescribed problem.
    volume      = pi*(D/2)^2 * L;
    mass        = (massDensity)*volume;
    IL          = mass*(D/2)^2/2;    
    I_body      = [IL 0; 0 IL];
              
    kx          = (pi)^2*mass;
    ky          = kx;
    ktheta      = (2*pi)^2*IL;
    stiff  = [kx ky ktheta];
    
    % Analytical body type:
    flag_forceinside = SM_FALSE; %SM_TRUE; % Force inside.
    annbody_type     = RB_ANNDISC;
    annbody_nparam   = 3;
    annbody_param    = [KAXIS D L]; % Cyl direction is always KAXIS. L = 1.

    %% Restraints:
    % Restrained Surface:
    nrsurf      =  1;
    fis_nodes_A = [];
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
    
    % Restrained DOFS of Node 1:
    restnode      = 1;
    nrestcoords   = 4; % All 3 + 1 dofs of node 1 are restrained.
    maxrestparams = 4;
    restnodes     = restnode*ones(1,nrestcoords);
    restdofs      = [1 2 3 4];
    resttype      = SM_PK_FIXED*ones(1,nrestcoords);
    nparam        = 1*ones(1,nrestcoords);
    param         = zeros(maxrestparams,nrestcoords); % All dofs set to zero
        
    % Change restraint on z displacement 
    if (strcmp(kinemflag(1,:),'fixd'))
        
        % Locations in x,y,z
        param(1,IAXIS) = xpos;
        param(1,JAXIS) = ypos;
        param(1,KAXIS) = zpos;

    elseif (strcmp(kinemflag(ibd,:),'rt_x'))
        
        % Locations in x,y,z
        param(1,IAXIS) = xpos;
        param(1,JAXIS) = ypos;
        param(1,KAXIS) = zpos;        
        
        % Rotate angle phi with constant velocity:
        resttype(ITHETA) = SM_PK_CONSTVEL; nparam(ITHETA) = 2;
        phio = 0; phido=-5; param(1:2,ITHETA) = [phio phido]; %SO OMEGA = 5 -> alpha=OMG*D/(2*Uoo) = 2.5        

    elseif (strcmp(kinemflag(ibd,:),'blck'))
        
        % Locations in x,y,z
        param(1,IAXIS) = xpos;
        param(1,JAXIS) = ypos;
        
        % do sin in z:
        dof_to_harm = JAXIS;
        resttype(dof_to_harm)    = SM_PK_HARMONIC;
        % Properties are such that z(t) = At*sin(2*pi*f*t+phase) + fixed_coord
        Ao=0.5; T=pi; f=1/T; phase = 0; fixed_coord=ypos;
        nparam(dof_to_harm) = 4;
        param(1:4,dof_to_harm) = [Ao f phase fixed_coord]';

        % Rotate angle phi with constant velocity:
        dof_to_harm = ITHETA;
        resttype(dof_to_harm)    = SM_PK_HARMONIC;
        % Properties are such that phi(t) = Ath*sin(2*pi*f*t+phase) + fixed_coord
        Ao=1; T=pi; f=1/T; phase = pi; fixed_coord=0.;
        nparam(dof_to_harm) = 4;
        param(1:4,dof_to_harm) = [Ao f phase fixed_coord]';
        
        
    elseif (strcmp(kinemflag(ibd,:),'elst'))
        
        % Locations in x,y,z
        param(1,IAXIS)  = xpos;
        param(1,ITHETA) = 0.;
        
        % do sin in y:
        dof_to_harm = JAXIS;
        resttype(dof_to_harm)    = SM_PK_HARMONIC;
        % Properties are such that z(t) = At*sin(2*pi*f*t+phase) + fixed_coord

        % Case (a1) Table 2. Elston, Blackburn, Sheridan JFM 2006.
        if (strcmp(ELSTON_CASE,'a1'))
          KC=8;
          beta=12.5;
        elseif (strcmp(ELSTON_CASE,'a3'))
          KC=2.5;
          beta=100;
        end
        
        Re=KC*beta
        Ao=KC*D/(2*pi)
        nu=1/Re
        f=nu*beta/D^2
        T=1/f;
        
        disp(['2pi*f*A=' num2str(2*pi*f*Ao)])
        
        phase = 0; fixed_coord=ypos;
        nparam(dof_to_harm) = 4;
        param(1:4,dof_to_harm) = [Ao f phase fixed_coord]';        
        
        
    end

    end
    
    %% Export Solid Model to file:
    % Values to write:
    TRANSFORMATION  =   RB_TWODIM;  % Use 2D Rotation matrix
    nnp             =      nnodes;  % Total number of nodes
    if (ibd ==1)
    nel_load        =         nel;  % nel surface elements
    nel             =           1;  % 1 rigid body element
    end
    ned             =     NDIM;
    max_eltype      =       15;  % Max nodes per element type (that is one node elem)
    max_eltype_load =        1;  % Max nodes per surface element (01 is 2 node segment)
    eltype          =       15*ones(1,nel); % Element type : One node Rigid Body
    eltype_load     =        1*ones(1,nel_load); % Surface Element type : Line
    kinematics_idx   =       1;  %
    
    x = XY(:,IAXIS); % Positions of all nodes : nnodes
    y = XY(:,JAXIS);
    
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
    
    plot(XY(:,IAXIS)+xpos,XY(:,JAXIS)+ypos,'xb')
    for iel=1:nel_load
      plot([XY(ws_IEN(iel,1),IAXIS)+xpos XY(ws_IEN(iel,2),IAXIS)+xpos], ...
           [XY(ws_IEN(iel,1),JAXIS)+ypos XY(ws_IEN(iel,2),JAXIS)+ypos],'-k')
    end

    
  end
end

xlabel('x')
ylabel('y')
axis equal

return
