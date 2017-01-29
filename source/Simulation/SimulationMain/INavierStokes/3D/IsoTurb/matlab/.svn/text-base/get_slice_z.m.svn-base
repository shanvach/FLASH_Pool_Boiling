function [DATA]=get_slice_z(file)

IAXIS = 1;
JAXIS = 2;
KAXIS = 3;
MDIM  = 3;

U_IND     = 1;
V_IND     = 2;
W_IND     = 3;  
DUDX_IND  = 4;
DUDY_IND  = 5; 
DUDZ_IND  = 6; 
DVDX_IND  = 7; 
DVDY_IND  = 8; 
DVDZ_IND  = 9; 
DWDX_IND  =10; 
DWDY_IND  =11; 
DWDZ_IND  =12; 
PRES_IND  =13;
TVIS_IND  =14;


%---Define precision format.
iformat = 'int32';
rformat = 'float64';

%---Read file
fid = fopen(file, 'r');

% Read Time, and molecular viscosity:
temp = fread(fid,   1, iformat);
vars = fread(fid, [4], rformat);
temp = fread(fid,   1, iformat);

DATA.time = vars(1);
DATA.nu   = vars(2);
DATA.zc   = vars(3);
DATA.Zslci_original = vars(4);

% Read istep,blockCount,NX,NY,NZB:
temp = fread(fid, 1, iformat);
vars = fread(fid, [4], iformat);
temp = fread(fid, 1, iformat);

DATA.istep =      vars(1);
NX = vars(2); NY = vars(3); TotVars = vars(4); 

DATA.NX =            NX;
DATA.NY =            NY;
DATA.TotVars =    TotVars;

% Read xmin, xmax, ymin, ymax:
temp = fread(fid, 1, iformat);
vars = fread(fid, [4], rformat);
temp = fread(fid, 1, iformat);
DATA.xmin = vars(1);
DATA.xmax = vars(2);
DATA.ymin = vars(3);
DATA.ymax = vars(4);

DATA.del(IAXIS) = (DATA.xmax - DATA.xmin)/DATA.NX;
DATA.del(JAXIS) = (DATA.ymax - DATA.ymin)/DATA.NY;

DATA.xcell = zeros(DATA.NX,1);
DATA.xedge = zeros(DATA.NX+1,1);

DATA.ycell = zeros(DATA.NY,1);
DATA.yedge = zeros(DATA.NY+1,1);

% Generate coordinates:
for i=1:DATA.NX
    DATA.xedge(i,1)=DATA.xmin + (i-1)*DATA.del(IAXIS);
    DATA.xcell(i,1)=DATA.xedge(i,1) + DATA.del(IAXIS)/2; 
end
DATA.xedge(DATA.NX+1,1)=DATA.xmin + DATA.NX*DATA.del(IAXIS);

for j=1:DATA.NY
    DATA.yedge(j,1)=DATA.ymin + (j-1)*DATA.del(JAXIS);
    DATA.ycell(j,1)=DATA.yedge(j,1) + DATA.del(JAXIS)/2; 
end
DATA.yedge(DATA.NY+1,1)=DATA.ymin + DATA.NY*DATA.del(JAXIS);

% Load Variables:
temp = fread(fid, 1, iformat);
DATA.u = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.v = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.w = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dudx = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dudy = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dudz = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dvdx = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dvdy = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dvdz = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dwdx = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dwdy = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.dwdz = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.p = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);

temp = fread(fid, 1, iformat);
DATA.nut = fread(fid, [DATA.NX DATA.NY], rformat);
temp = fread(fid, 1, iformat);


fclose(fid);

return