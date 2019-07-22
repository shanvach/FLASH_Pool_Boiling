close all
clear all
clc


IND_TIME = 2;
IND_ANG  = 3;
IND_X    = 3;
IND_Y    = 4;

basedir='/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/OSCCYL_2D_fineOPT/IOData/';
filename_ang='posvelacc_ang.000001.res';

[vect]=load([basedir filename_ang],'r');
time_ang=vect(:,IND_TIME);
ang     =vect(:,IND_ANG);
vel_ang =vect(:,IND_ANG+1);
acc_ang =vect(:,IND_ANG+2);

figure
subplot(3,1,1)
plot(time_ang,ang)
ylabel('ang')
grid on

subplot(3,1,2)
plot(time_ang,vel_ang)
ylabel('Vang')
grid on

subplot(3,1,3)
plot(time_ang,acc_ang)
ylabel('Aang')
grid on
xlabel('time')


basedir='/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/OSCCYL_2D_fineOPT/IOData/';
basedir2='/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/OSCCYL_2D_fineOPT/IOData/UDERIVS/';
filename_x='posvelacc_x.000001.res';

[vect]=load([basedir filename_x],'r');
time_x=vect(:,IND_TIME);
x     =vect(:,IND_Y);
vel_x =vect(:,IND_Y+2);
acc_x =vect(:,IND_Y+4);

[vect]=load([basedir2 filename_x],'r');
time_x2=vect(:,IND_TIME);
x2     =vect(:,IND_Y);
vel_x2 =vect(:,IND_Y+2);
acc_x2 =vect(:,IND_Y+4);

figure
subplot(3,1,1)
hold on
plot(time_x,x)
plot(time_x2,x2,'r')
ylabel('x')
grid on

subplot(3,1,2)
hold on
plot(time_x,vel_x)
plot(time_x2,vel_x2,'r')
ylabel('Vx')
grid on

subplot(3,1,3)
hold on
plot(time_x,acc_x)
plot(time_x2,acc_x2,'r')
ylabel('Ax')
grid on
xlabel('time')
