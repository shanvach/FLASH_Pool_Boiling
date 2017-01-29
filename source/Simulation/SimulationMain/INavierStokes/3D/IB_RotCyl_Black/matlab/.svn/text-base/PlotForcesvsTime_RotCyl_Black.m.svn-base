% File to plot Forces in Z direction for Oscillating Cylinder:
% 
%
%         % Locations in x,y,z
%         param(1,IAXIS) = xpos;
%         param(1,JAXIS) = ypos;
%         
%         % do sin in z:
%         dof_to_harm = KAXIS;
%         resttype(dof_to_harm)    = SM_PK_HARMONIC;
%         % Properties are such that z(t) = At*sin(2*pi*f*t+phase) + fixed_coord
%         Ao=0.5; T=pi; f=1/T; phase = 0; fixed_coord=zpos;
%         nparam(dof_to_harm) = 4;
%         param(1:4,dof_to_harm) = [Ao f phase fixed_coord]';
% 
%         % Rotate angle phi with constant velocity:
%         dof_to_harm = IPHI;
%         resttype(dof_to_harm)    = SM_PK_HARMONIC;
%         % Properties are such that phi(t) = Ath*sin(2*pi*f*t+phase) + fixed_coord
%         Ao=1; T=pi; f=1/T; phase = pi; fixed_coord=0.;
%         nparam(dof_to_harm) = 4;
%         param(1:4,dof_to_harm) = [Ao f phase fixed_coord]';     
%
% -------------------------------------------------------------------------
close all
clear all
clc

% Strouhal Number:
rho = 1;
Uoo = 1;
D   = 1;
KC  = pi;
f   = 1/KC;
T   = 1/f; 
Re = 90*pi;

% Cylinder Lengths:
L_ug  = 0.0625;
L_amr = 0.625;

% Axes:
IAXIS =1;
JAXIS =2;
KAXIS =3;

% Basedirs for New computations:

% Uniform Grid and 3 Level AMR:
basedir_ug  ='/Users/mvanella/Dropbox/Documents/ADAPTIVE/CHRIS_OSCCYL/ug/'; 
basedir_amr ='/Users/mvanella/Dropbox/Documents/ADAPTIVE/CHRIS_OSCCYL/amr/';
% 
filename='force.000001.res';


% Load:
[vect]=load([basedir_ug filename],'r');
nstep_ug = vect(:,  1);
time_ug  = vect(:,  2);
Fpres_ug = vect(:,3:5)/L_ug;
Fvisc_ug = vect(:,6:8)/L_ug;

[vect2]=load([basedir_amr filename],'r');
nstep_amr = vect2(:,  1);
time_amr  = vect2(:,  2);
Fpres_amr = vect2(:,3:5)/L_amr;
Fvisc_amr = vect2(:,6:8)/L_amr;


% Normalization
q_A = 1/2*rho*Uoo^2*(pi*(D/2)^2);


istart =   10;
istep  =    1;

% Figures:
figure
hold on
plot(time_ug(istart:istep:end)/T,Fpres_ug(istart:istep:end,KAXIS)/q_A,'r')
plot(time_ug(istart:istep:end)/T,Fvisc_ug(istart:istep:end,KAXIS)/q_A,'b')

plot(time_amr(istart:istep:end)/T,Fpres_amr(istart:istep:end,KAXIS)/q_A,'--r')
plot(time_amr(istart:istep:end)/T,Fvisc_amr(istart:istep:end,KAXIS)/q_A,'--b')

ylabel('C_D','FontSize',14)
legend('UG Pressure','UG Viscous','AMR Pressure','AMR Viscous')
set(gca,'FontSize',14)
axis([0. 5. -10. 10.])
grid on
box on


return
