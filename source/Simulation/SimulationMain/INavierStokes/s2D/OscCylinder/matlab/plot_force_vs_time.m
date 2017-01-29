% Script for plot Cfy for oscillating cylinder in Elston et al. (JFM 2006) 
% cases of Table 2 on that Journal.
close all
clear all
clc

D = 1;
rho=1;

IND_TIME       = 2;
IND_FX_PRES    = 3;
IND_FY_PRES    = 4;
IND_FX_VISC    = 5;
IND_FY_VISC    = 6;

ELSTON_CASE = 'a1';

%% Computation with 64 cells in diameter.
if (strcmp(ELSTON_CASE,'a1'))
  KC=8;
  beta=12.5;
  Re=KC*beta;
  Ao=KC*D/(2*pi); 
  nu=1/Re;
  f=nu*beta/D^2; T=1/f;
  Ums = ((2*pi*f*Ao)/sqrt(2))^2;
  basedir='/Users/mvanella/Dropbox/Documents/ADAPTIVE/CHRIS_OSCCYL/OSCCYL_2D_ELSTON_A1/';

elseif (strcmp(ELSTON_CASE,'a3'))
  KC=2.5;
  beta=100;
  Re=KC*beta;
  Ao=KC*D/(2*pi); Ums = (Ao/sqrt(2))^2;
  nu=1/Re;
  f=nu*beta/D^2; T=1/f;
  Ums = ((2*pi*f*Ao)/sqrt(2))^2;
  basedir='';

end
filename_force='force.000001.res';

qaD = 1/2*rho*Ums*D;


[vect]=load([basedir filename_force],'r');
time    = vect(:,IND_TIME);
Cfx_pres = vect(:,IND_FX_PRES)/qaD;
Cfy_pres = vect(:,IND_FY_PRES)/qaD;
Cfx_visc = vect(:,IND_FX_VISC)/qaD;
Cfy_visc = vect(:,IND_FY_VISC)/qaD;

figure
subplot(1,2,1)
hold on
plot(time/T,Cfy_pres,'--r')
plot(time/T,Cfy_visc,'--b')
plot(time/T,Cfy_pres+Cfy_visc,'k')
ylabel('Cfy','FontSize',14)
xlabel('time/T','FontSize',14)
set(gca,'FontSize',12)
legend('pres','visc','total')
axis([0 3 -8 8])
grid on

subplot(1,2,2)
hold on
plot(time/T,Cfx_pres,'--r')
plot(time/T,Cfx_visc,'--b')
plot(time/T,Cfx_pres+Cfx_visc,'k')
ylabel('Cfx','FontSize',14)
xlabel('time/T','FontSize',14)
set(gca,'FontSize',12)
legend('pres','visc','total')
axis([0 3 -.8 .8])
grid on

%% Compare to Fine grid solution 96 cells in D:

if (strcmp(ELSTON_CASE,'a1'))
  basedir2='/Users/mvanella/Dropbox/Documents/ADAPTIVE/CHRIS_OSCCYL/OSCCYL_2D_ELSTON_A1_96D/';

elseif (strcmp(ELSTON_CASE,'a3'))
  basedir2='';

end

[vect2]=load([basedir2 filename_force],'r');
time2    = vect2(:,IND_TIME);
Cfx_pres2 = vect2(:,IND_FX_PRES)/qaD;
Cfy_pres2 = vect2(:,IND_FY_PRES)/qaD;
Cfx_visc2 = vect2(:,IND_FX_VISC)/qaD;
Cfy_visc2 = vect2(:,IND_FY_VISC)/qaD;

figure
hold on
plot(time/T,Cfy_pres,'--r')
plot(time/T,Cfy_visc,'--b')
plot(time/T,Cfy_pres+Cfy_visc,'--k')

legend('pres','visc','total')

plot(time2/T,Cfy_pres2,'r')
plot(time2/T,Cfy_visc2,'b')
plot(time2/T,Cfy_pres2+Cfy_visc2,'k')

ylabel('Cfy','FontSize',14)
xlabel('time/T','FontSize',14)
set(gca,'FontSize',12)
legend('pres','visc','total')
axis([0 3 -8 8])
grid on


