% Read Turbkinlam.res and plot:

%close all
clear all
clc

TIME=2;
TKE =3;

basedir='/Users/mvanella/Dropbox/Documents/ADAPTIVE/PARTICLES/ISOTURB_384/';


V1=load([basedir 'TURBKINLAM.res']);

end_stat1 = length(V1(:,1));
%end_stat2 = length(V2(:,1));


% Figure
figure
hold on
plot(V1(1:end_stat1,TIME),V1(1:end_stat1,TKE))
%plot(V2(1:end_stat2,TIME),V2(1:end_stat2,TKE),'--r')
xlabel('t','FontSize',14)
ylabel('tke','FontSize',14)
set(gca,'FontSize',12)
%legend('Euler','RK2')
grid on
box on
axis([0 250 0 2.5])