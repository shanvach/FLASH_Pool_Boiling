close all
clear all
%clc
format long

basedir ='/Users/mvanella/Documents/ADAPTIVE/FLASH/test/INS_IB/XYZ_POISSON_SUPERLU/IOData/';

nprocs = 1;
nblocks= 1;
nxb    = 32; nyb    = 32; nzb    = 32;
n_glob = nxb*nyb*nzb*nblocks;

% Read A:
A=sparse(n_glob,n_glob);
for ip=1:nprocs
    filename=['matrix_' num2str(ip-1,'%4.4d') '.res']
    [VEC] = load([basedir filename]);    
    nnz = length(VEC(:,1));
    for ii=1:nnz
       i   = VEC(ii,1);
       j   = VEC(ii,2);
       val = VEC(ii,3);
       A(i,j)=val;
    end
end
spy(A)

% Read rhs:
rhs=zeros(n_glob,1);
for ip=1:nprocs
    filename=['rhs_' num2str(ip-1,'%4.4d') '.res']
    [VEC] = load([basedir filename]);    
    nnz = length(VEC(:,1));
    for ii=1:nnz
       i        = VEC(ii,1);
       rhs(i)   = VEC(ii,2);
    end
end

sol = A\rhs;

% Check against analytical res:
xob = 0; xeb = 1; Lx = (xeb-xob); dx = Lx/nxb;
yob = 0; yeb = 1; Ly = (yeb-yob); dy = Ly/nyb;
zob = 0; zeb = 1; Lz = (zeb-zob); dz = Lz/nzb;
wx = 1; wy = 1; wz = 1;
ii =0;
phiAnn = zeros(1,n_glob);
phiNum = zeros(1,n_glob);
for k=1:nzb
    zc = zob + (k-1+0.5)*dz; 
    for j=1:nyb
        yc = yob + (j-1+0.5)*dy;
        for i=1:nyb
            xc = xob + (i-1+0.5)*dx;
            ii = ii + 1;            
            phiAnn(ii) = sin(2*pi*wx*xc/Lx)*sin(2*pi*wy*yc/Ly)*sin(2*pi*wz*zc/Lz);      
            phiNum(ii) = sol(ii);
        end
    end
end
     
errinf = max(abs(phiAnn - phiNum))
