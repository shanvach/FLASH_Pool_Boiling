% Script to write Rigid Body h5 file, to be used by Splash.
% Name of file for simulation: sm_body.1.h5
% The Rigid is oriented using a 321 Euler angle sequence.
% -------------------------------------------------------------------------
close all
clear all
clc

%% Define variables:
IAXIS = 1;
JAXIS = 2;
KAXIS = 3;
IPHI  = 4;
NDIM  = 3; % Dimensions

SM_FALSE = 0;
SM_TRUE  = 1;

% Rigid Body constants:
BODYTYPE_RIGID  = 1;

RB_IDENTITY = 1000;
RB_EULER321 = 1321;
RB_QUATERNN = 1111;

RB_ANNSPHERE = 55;
RB_ANNDISC   = 56;
RB_ANNRBC    = 57;


% Kinematics Constants:       
SM_PK_FIXED    =   0;
SM_PK_HARMONIC = 102;
SM_PK_CONSTVEL = 103;

DOFS_per_node   =  9;  % Degrees of freedom per node: 
                       % x,y,z,ang1,ang2,ang3,wx,wy,wz

ix = 1; ex = 3;
ia = 4; ea = 6;
iw = 7; ew = 9;

%% Define variables:
NumBods = 1; 

np_el   = 3; % Points per wet surface element: 3 - Triangle

% Output File names:
basedir_out = ['/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/source/Simulation/SimulationMain/INavierStokes/3D/IB_fallingSphere/matlab/'];

%hfilename=['sm_bodySphere.1.h5'; ...
%           'sm_bodySphere.2.h5'];
hfilename=['sm_body.00001.h5'];

% Prescribed Kinematics Per Body: 
kinemflag=['    '];

%% Read Mesh File:

% Base directory:
basedir='/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/source/Simulation/SimulationMain/INavierStokes/3D/IB_fallingSphere/matlab/';

% Input File name:
readflg= 'gambit'; %'dat'; 
filename='sphere06.neu';

% Read File
if (strcmp(readflg,'dat')) % text file of nnodes, nelem and XYZ followed
                           % by triangle connectivities.
    [fid]=fopen([basedir filename],'r');
    [nnodes] = str2num(fgetl(fid));
    [nel]    = str2num(fgetl(fid));
    XYZ   = zeros(nnodes+1,NDIM);
    ws_IEN= zeros(nel,np_el);
    % Set the first node to {0,0,0} origin of reference frame:
    for i=2:nnodes+1
        XYZ(i,:) = str2num(fgetl(fid));
    end
    for i=1:nel
        ws_IEN(i,:) = str2num(fgetl(fid)) + [1 1 1];
    end
    fclose(fid);
elseif(strcmp(readflg,'gambit'))
    [XYZ,ws_IEN,nnodes,nel]=readsurf_gambit([basedir filename]);
end

% Plot figure:
figure
hold on
trimesh(ws_IEN,XYZ(1:nnodes+1,IAXIS),XYZ(1:nnodes+1,JAXIS),XYZ(1:nnodes+1,KAXIS))
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

% Test normals out:
for iel=1:nel
   x12 = XYZ(ws_IEN(iel,2),:) - XYZ(ws_IEN(iel,1),:);
   x13 = XYZ(ws_IEN(iel,3),:) - XYZ(ws_IEN(iel,1),:);
   vcr = cross(x12,x13);
   xcen = 1/3*(XYZ(ws_IEN(iel,1),:)+XYZ(ws_IEN(iel,2),:)+XYZ(ws_IEN(iel,3),:));
   sign_normal(iel) = sign(dot(xcen,vcr));
   if sign_normal(iel) < 0
       aux=ws_IEN(iel,3); ws_IEN(iel,3)=ws_IEN(iel,2); ws_IEN(iel,2)=aux;
   end
end

for iel=1:nel
   x12 = XYZ(ws_IEN(iel,2),:) - XYZ(ws_IEN(iel,1),:);
   x13 = XYZ(ws_IEN(iel,3),:) - XYZ(ws_IEN(iel,1),:);
   vcr = cross(x12,x13);
   xcen = 1/3*(XYZ(ws_IEN(iel,1),:)+XYZ(ws_IEN(iel,2),:)+XYZ(ws_IEN(iel,3),:));
   sign_normal2(iel) = sign(dot(xcen,vcr));
end


%return

for ibd=1:NumBods

    %% Rigid body properties:
    grav_vec     = [0 0 0];
    grav_flag    =  0;
    
    
    %% Inertia properties:
    rho = 2.56;
    rhof= 1.00;
    
%     % Case ellipse:
%     a_ell = 1/5;
%     b_ell = 1/2;
%     c_ell = 1/8;
%     Vol = 4/3*pi*a_ell*b_ell*c_ell;
%     mass = rho*Vol;
% 
%     Ixx = 1/5*mass*(b_ell^2+c_ell^2); Ixy = 0.0; Ixz = 0.0;
%     Iyx = Ixy; Iyy = 1/5*mass*(a_ell^2+c_ell^2); Iyz = 0.0;
%     Izx = Ixz; Izy = Iyz; Izz = 1/5*mass*(a_ell^2+b_ell^2);
% 
%     I_body = [Ixx Ixy Ixz; Iyx Iyy Iyz; Izx Izy Izz];

    
    massDensity = rho;
    % Case sphere:
    D           = 1.0;  % Sphere Diameter to compute sphere mass.
    volume      = 4/3*pi*(D/2)^3;
    mass        = (massDensity)*volume;
    Ixx         = 2/5*mass*(D/2)^2;
    Iyy         = Ixx;
    Izz         = Ixx;
    I_body = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
    
    kx          = 0.;
    ky          = kx;
    kz          = kx;
    stiff  = [kx ky kz];
    
    % Analytical body type:
    flag_forceinside = SM_FALSE; %TRUE; % Force inside.
    annbody_type     = RB_ANNSPHERE;
    annbody_nparam   = 1;
    annbody_param    = [D];

    %% Restraints:
    % Restrained Surface:
    nrsurf      =  1;
    fis_nodes_A = [];
    fix_nodes_A = [2:nnodes+1]; % Restrained nodes in global numbering
    nfix = length(fix_nodes_A);
    
    % Restrained DOFS of Node 1:
    % Restrained DOFS of Node 1:
    restnode      = 1;
    nrestcoords   = 9; % All six + 3 dofs of node 1 are restrained.
    maxrestparams = 2;
    restnodes     = restnode*ones(1,nrestcoords);
    restdofs      = [1 2 3 4 5 6 7 8 9];
    resttype      = SM_PK_FIXED*ones(1,nrestcoords);
    nparam        = 1*ones(1,nrestcoords);
    param         = zeros(maxrestparams,nrestcoords); % All dofs set to zero
        
    % Change restraint on z displacement 
    if (strcmp(kinemflag(ibd,:),'sinz'))
        % do sin in z:
        for dof_to_harm = 3:3
            resttype(dof_to_harm)    = SM_PK_HARMONIC;
            % Properties are such that phi(t) = Ao*sin(2*pi*f*t+phase) + fixed_coord
            % Ao=1/0.625; T=2*pi/0.625; f=1/T; phase = pi/2; fixed_coord=0.;
            Ao=0.5; T=0.5*pi; f=1/T; phase = pi/2; fixed_coord=0.;
            nparam(dof_to_harm) = 4;
            param(1:4,dof_to_harm) = [Ao f phase fixed_coord]';
        end

    elseif (strcmp(kinemflag(ibd,:),'restz'))
        
        disp('Rest z ...')
        % Z Direction
        resttype(KAXIS) = SM_PK_CONSTVEL; nparam(KAXIS) = 2;
        zo = 0.; wo=-1.; param(1:2,KAXIS) = [zo wo];        
        
    elseif (strcmp(kinemflag(ibd,:),'cv_1'))
        % First type of CONSTVEL 
        % Y direction
        resttype(JAXIS) = SM_PK_CONSTVEL; nparam(JAXIS) = 2;
        yo = 0.; vo=-1/sqrt(2); param(1:2,JAXIS) = [yo vo];
        % Z Direction
        resttype(KAXIS) = SM_PK_CONSTVEL; nparam(KAXIS) = 2;
        zo = 0.; wo=-1/sqrt(2); param(1:2,KAXIS) = [zo wo];
        
    elseif (strcmp(kinemflag(ibd,:),'cv_2'))
        % Second type of CONSTVEL 
        % Y direction
        resttype(JAXIS) = SM_PK_CONSTVEL; nparam(JAXIS) = 2;
        yo = -1.5; vo=1/sqrt(2); param(1:2,JAXIS) = [yo vo];
        % Z Direction
        resttype(KAXIS) = SM_PK_CONSTVEL; nparam(KAXIS) = 2;
        zo = 0.; wo=1/sqrt(2); param(1:2,KAXIS) = [zo wo];        

    elseif (strcmp(kinemflag(ibd,:),'rt_x'))
        % Rotate angle phi with constant velocity:
        resttype(IPHI) = SM_PK_CONSTVEL; nparam(IPHI) = 2;
        phio = 0; phido=2*pi; param(1:2,IPHI) = [phio phido];
        
    end

    %% Export Solid Model to file:
    % Values to write:
    TRANSFORMATION  = RB_EULER321;  % Use Euler 321
    nnp             =    nnodes+1;  % Total number of nodes
    if (ibd ==1)
    nel_load        =         nel;  % nel surface elements
    nel             =           1;  % 1 rigid body element
    end
    ned             =     NDIM;
    max_eltype      =       15;  % Max nodes per element type (that is one node elem)
    max_eltype_load =        2;  % Max nodes per surface element (02 is triangle)
    eltype          =       15*ones(1,nel); % Element type : One node Rigid Body
    eltype_load     =        2*ones(1,nel_load); % Surface Element type : Triangle
    kinematics_idx   =       1;  %
    
    x = XYZ(:,IAXIS); % Positions of all nodes : nnodes + 1
    y = XYZ(:,JAXIS);
    z = XYZ(:,KAXIS);
    
    IEN = [1];  % One rigid body with node 1 as connectivity
    
    
    % Write hdf5 file:
    % Mesh:
    hfile=[basedir_out hfilename(ibd,:)];
    fprintf(1,'     hdf5 file started...\n');
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
    hdf5write(hfile,'body/Mass',mass,'WriteMode','append')
    hdf5write(hfile,'body/Volume',volume,'WriteMode','append')
    hdf5write(hfile,'body/I_body',I_body,'WriteMode','append')
    hdf5write(hfile,'body/trmatrix',int32(TRANSFORMATION),'WriteMode','append')
    hdf5write(hfile,'body/Stiffness',stiff,'WriteMode','append')
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
    hdf5write(hfile,'RestNodes/param'        ,param,'WriteMode','append')
    
    fprintf(1,'   Body %d hdf5 file done.\n',ibd);

end

return