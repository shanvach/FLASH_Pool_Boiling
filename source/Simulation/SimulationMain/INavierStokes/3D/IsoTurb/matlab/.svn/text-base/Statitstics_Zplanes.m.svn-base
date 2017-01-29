% Statistics_multiphase:
% Compute statistics as a function of z for simulations of IsoTurb Flash
% case. Use two methods:
% 1. Read the statistics space averaged by planes and average in time, in
% function read_stats.m . The corresponding files are named :
% statsz.00istat.00iproc.bin
% 2. Read the Raw data on defined plane on z and compute statistics in
% space-time.
% -------------------------------------------------------------------------
close all
clear all
clc

%% Define parameters for stats:
% DNS 2pi in COLONIAL:
stats_start    =   35;
stats_end      =   35;
NumProcs       =  256;
stats_nslicesz =   16;
invRe          =0.001;
nu = invRe;
basedir='/Users/mvanella/Dropbox/Documents/ADAPTIVE/PARTICLES/ISOTURB_384/';



%% 1. Read stats already computed and averaged in time in IsoTurb/ins_turbstats_z.F90 
% Warning, these stats are only computed for uniform grid UG unit with procs in z.
%[DATAavg]=read_stats(basedir,nu,stats_start,stats_end,NumProcs);


%% 2. Read slice info and compute statistics. Files produced by
%     IsoTurb/ins_exportslicces_Z.F90
[DATAavg_slice] = slicez_stats(basedir,stats_start,stats_end,stats_nslicesz);


%%
% Plot time averaged velocities:
figure
subplot(1,3,1)
hold on
% for iproc = 1 : NumProcs
%    for lb = 1 : DATAavg(iproc).blockCount 
%      plot(DATAavg(iproc).uavg(:,lb),DATAavg(iproc).zc(:,lb),'b')
%    end
% end
for islice = 1 : stats_nslicesz
    plot(DATAavg_slice(islice).uavg,DATAavg_slice(islice).zc,'.k')
end
xlabel('uavg')
ylabel('z')
grid on
box on

subplot(1,3,2)
hold on
% for iproc = 1 : NumProcs
%    for lb = 1 : DATAavg(iproc).blockCount 
%      plot(DATAavg(iproc).vavg(:,lb),DATAavg(iproc).zc(:,lb),'b')
%    end
% end
for islice = 1 : stats_nslicesz
    plot(DATAavg_slice(islice).vavg,DATAavg_slice(islice).zc,'.k')
end
xlabel('vavg')
ylabel('z')
grid on
box on

subplot(1,3,3)
hold on
% for iproc = 1 : NumProcs
%    for lb = 1 : DATAavg(iproc).blockCount 
%      plot(DATAavg(iproc).q2(:,lb)/2,DATAavg(iproc).zc(:,lb),'b')
%    end
% end
for islice = 1 : stats_nslicesz
    plot(DATAavg_slice(islice).q2/2,DATAavg_slice(islice).zc,'.k')
end

xlabel('k')
ylabel('z')
grid on
box on

figure
subplot(1,3,1)
hold on
% for iproc = 1 : NumProcs
%    for lb = 1 : DATAavg(iproc).blockCount 
%      plot(sqrt(DATAavg(iproc).urms2(:,lb)),DATAavg(iproc).zc(:,lb),'b')
%      plot(sqrt(DATAavg(iproc).vrms2(:,lb)),DATAavg(iproc).zc(:,lb),'--r')
%      plot(sqrt(DATAavg(iproc).wrms2(:,lb)),DATAavg(iproc).zc(:,lb),':k')
% 
%      plot(sqrt(DATAavg(iproc).uu(:,lb)),DATAavg(iproc).zc(:,lb),'--k')
%      
%    end
% end
for islice = 1 : stats_nslicesz
    plot(1/3*(sqrt(DATAavg_slice(islice).urms2)+ ...
              sqrt(DATAavg_slice(islice).vrms2)+ ...
              sqrt(DATAavg_slice(islice).wrms2)),DATAavg_slice(islice).zc,'.k')
end

xlabel('urms')
ylabel('z')
grid on
box on

subplot(1,3,2)
hold on
% for iproc = 1 : NumProcs
%    for lb = 1 : DATAavg(iproc).blockCount 
%      plot(DATAavg(iproc).Relambda(:,lb),DATAavg(iproc).zc(:,lb),'b')
%    end
% end
for islice = 1 : stats_nslicesz
    plot(DATAavg_slice(islice).Relambda,DATAavg_slice(islice).zc,'.k')
end
xlabel('Re_\lambda')
ylabel('z')
grid on
box on

subplot(1,3,3)
hold on
% for iproc = 1 : NumProcs
%    for lb = 1 : DATAavg(iproc).blockCount 
%      plot(DATAavg(iproc).lambdaf(:,lb),DATAavg(iproc).zc(:,lb),'b')
%    end
% end
for islice = 1 : stats_nslicesz
    plot(DATAavg_slice(islice).lambdag,DATAavg_slice(islice).zc,'.k')
end
xlabel('Lambda_g')
ylabel('z')
grid on
box on


figure
for islice = 1 : stats_nslicesz
  hold off  
  plot(DATAavg_slice(islice).kx,DATAavg_slice(islice).Elong(1:floor(DATAavg_slice(islice).NX/2)+1))
  hold on
  plot(DATAavg_slice(islice).kx,DATAavg_slice(islice).kx.^(-5/3),'--k')
  title(['Slice=' num2str(islice) ', zc=' num2str(DATAavg_slice(islice).zc)])
  xlabel('k')
  ylabel('Euu')
  set(gca,'Xscale','log','Yscale','log','FontSize',12)
  axis([1 1000 10^-17 10^-1])
  pause
end

% Mean urms, lambdag, ReLambda:
urms=0; vrms=0; wrms=0;
lambdag=0;
ReLambda=0;
for islice = 1 : stats_nslicesz
    urms=urms+sqrt(DATAavg_slice(islice).urms2)/stats_nslicesz;
    vrms=vrms+sqrt(DATAavg_slice(islice).vrms2)/stats_nslicesz;
    wrms=wrms+sqrt(DATAavg_slice(islice).wrms2)/stats_nslicesz;
   
    lambdag=lambdag+DATAavg_slice(islice).lambdag/stats_nslicesz;
    ReLambda=ReLambda+DATAavg_slice(islice).Relambda/stats_nslicesz;
    
end
disp(['urms=' num2str(urms) ', vrms=' num2str(vrms) ', wrms=' num2str(wrms)])
disp(['lambdag=' num2str(lambdag) ', ReLambda=' num2str(ReLambda)])

save([basedir 'AVGDATA.mat'],'DATAavg_slice')


return