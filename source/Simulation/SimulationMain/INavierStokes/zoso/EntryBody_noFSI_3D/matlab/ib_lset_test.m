clear
close all
clc

basedir='/Users/Akash/Desktop/matlab/matlab_3D/';
filename= 'sphere_coarse';

[XYZ,ws_IEN,nnodes,nel,normal]=stlread([basedir filename]);

XYZ(:,3) = XYZ(:,3);

xcenter = zeros(nel,1);
ycenter = zeros(nel,1);
zcenter = zeros(nel,1);

nrm3 = zeros(nel,3);

for i=1:nel
    
    xcenter(i) = sum(XYZ(ws_IEN(i,:),1))/3;
    ycenter(i) = sum(XYZ(ws_IEN(i,:),2))/3;
    zcenter(i) = sum(XYZ(ws_IEN(i,:),3))/3;
                    
    PA   = [XYZ(ws_IEN(i,1),1) XYZ(ws_IEN(i,1),2) XYZ(ws_IEN(i,1),3)];                 
    PB   = [XYZ(ws_IEN(i,2),1) XYZ(ws_IEN(i,2),2) XYZ(ws_IEN(i,2),3)];                  
    PC   = [XYZ(ws_IEN(i,3),1) XYZ(ws_IEN(i,3),2) XYZ(ws_IEN(i,3),3)];
    
    nrm1 = cross(PB-PA,PC-PA);
    
    nrm3(i,:) = nrm1./norm(nrm1);
   
    
end

Ng = 20;
[Xg,Yg,Zg] = meshgrid(linspace(-1.5,1.5,Ng),linspace(-1.5,1.5,Ng),linspace(-1.5,1.5,Ng));

phi  = zeros(size(Xg));
dist = zeros(nel,1);

for k=1:Ng
    for j=1:Ng
        for i=1:Ng
                   
            count  = 0;
            
            for ielem=1:nel
               
                    P1   = [Xg(i,j,k) Yg(i,j,k) Zg(i,j,k)];
                    P0   = [xcenter(ielem) ycenter(ielem) zcenter(ielem)];
                    PA   = [XYZ(ws_IEN(ielem,1),1) XYZ(ws_IEN(ielem,1),2) XYZ(ws_IEN(ielem,1),3)];
                    PB   = [XYZ(ws_IEN(ielem,2),1) XYZ(ws_IEN(ielem,2),2) XYZ(ws_IEN(ielem,2),3)];
                    PC   = [XYZ(ws_IEN(ielem,3),1) XYZ(ws_IEN(ielem,3),2) XYZ(ws_IEN(ielem,3),3)];
                    
                    nrm  = [nrm3(ielem,1) nrm3(ielem,2) nrm3(ielem,3)];
                 
                    min_vec = [min([PA(1) PB(1) PC(1)]) min([PA(2) PB(2) PC(2)]) min([PA(3) PB(3) PC(3)])];
                    max_vec = [max([PA(1) PB(1) PC(1)]) max([PA(2) PB(2) PC(2)]) max([PA(3) PB(3) PC(3)])];
                                        
                    lx = [1 0 0]; 
                    ln = nrm;
                    
                    vec  = PA-P1; 
                    
                    vecA = PB-PA;
                    vecB = PC-PA;
                    dotD = dot(vecA,vecB)^2 - dot(vecA,vecA)*dot(vecB,vecB);
                    
                    if(abs(dot(lx,nrm)) < 1e-13)
                        du = 0.0;
                    else
                        du = dot(vec,nrm)/dot(lx,nrm);
                    end
                    
                    PI = P1 + du.*lx;
                    
                    vecW = PI-PA;                 
                    
                    da = (dot(vecA,vecB)*dot(vecW,vecB) - dot(vecB,vecB)*dot(vecW,vecA))/dotD;
                    db = (dot(vecA,vecB)*dot(vecW,vecA) - dot(vecA,vecA)*dot(vecW,vecB))/dotD;
                    
                    if(da>=0 && da<=1 && db>=0 && da+db<=1 && du > 0)
                        count = count+1;
                    end
                                       
                    dn = dot(vec,nrm)/dot(ln,nrm);
                    
                    PN = P1+dn.*ln;
                                                      
                    vecW = PN-PA;                 
                    
                    da = (dot(vecA,vecB)*dot(vecW,vecB) - dot(vecB,vecB)*dot(vecW,vecA))/dotD;
                    db = (dot(vecA,vecB)*dot(vecW,vecA) - dot(vecA,vecA)*dot(vecW,vecB))/dotD;
                                                            
                    if(da>=0 && da<=1 && db>=0 && da+db<=1)
                        dist(ielem) = norm(P1-PN);
                    else
                        dist(ielem) = min([norm(P1-P0),norm(P1-PA),norm(P1-PB),norm(P1-PC)]); 
                    end
                                   
            end
            
            if(mod(count,2) == 1)
                phi(i,j,k) = min(dist);
            else
                phi(i,j,k) = -min(dist);
            end
            
        end
    end
end

%%
close all

figure
hold on
trimesh(ws_IEN,XYZ(1:nnodes+1,1),XYZ(1:nnodes+1,2),XYZ(1:nnodes+1,3))
quiver3(xcenter,ycenter,zcenter,nrm3(:,1),nrm3(:,2),nrm3(:,3),'k')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

zplane = 3;

figure
hold on
%trimesh(ws_IEN,XYZ(1:nnodes+1,1),XYZ(1:nnodes+1,2),XYZ(1:nnodes+1,3))
contourslice(Xg,Yg,Zg,phi,[],0,0,20)
isosurface(Xg,Yg,Zg,phi,0.0)
colorbar
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
