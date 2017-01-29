function [DATAavg] = slicez_stats(basedir,stats_start,stats_end,stats_nslicesz)
% Read and compute statistics from stat_slcz.00iwrite.0islice.bin
% 
%
% -------------------------------------------------------------------------

% close all
% clear all
% clc
% stats_start    =   1;
% stats_end      =   2;
% stats_nslicesz =   3;
% basedir='/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/ISOTURB_UG_DEBUG/IOData/';

IAXIS=1;
JAXIS=2;

nstats = stats_end - stats_start + 1; 
indexstat = [stats_start:stats_end];

%% Read data, by time ensemble and processor:
% for islice = 1 : stats_nslicesz
% 
%     for istat = 1 : nstats 
%         
%        file = [basedir 'stat_slcz.' num2str(indexstat(istat),'%4.4d') '.' ...
%                                     num2str(islice,'%3.3d') '.bin']; 
%         
%        [data]=get_slice_z(file);
%         
%        % Massage data:
%        DATA(istat,islice) = data;
%         
%     end 
%     disp(['Slice=' num2str(islice) ', zc=' num2str(DATA(istat,islice).zc)])
% end



%% Now compute statistics:
uavg = zeros(nstats,stats_nslicesz);
vavg = zeros(nstats,stats_nslicesz);
wavg = zeros(nstats,stats_nslicesz);
pavg = zeros(nstats,stats_nslicesz);

uu   = zeros(nstats,stats_nslicesz);
vv   = zeros(nstats,stats_nslicesz);
ww   = zeros(nstats,stats_nslicesz);
uv   = zeros(nstats,stats_nslicesz);
uw   = zeros(nstats,stats_nslicesz);
vw   = zeros(nstats,stats_nslicesz);

urms2= zeros(nstats,stats_nslicesz);
vrms2= zeros(nstats,stats_nslicesz);
wrms2= zeros(nstats,stats_nslicesz);

dudx = zeros(nstats,stats_nslicesz);
dudy = zeros(nstats,stats_nslicesz);
lambdaf = zeros(nstats,stats_nslicesz);
lambdag = zeros(nstats,stats_nslicesz);

nut     = zeros(nstats,stats_nslicesz);
eps_sgs = zeros(nstats,stats_nslicesz);
eps     = zeros(nstats,stats_nslicesz);
Relambda= zeros(nstats,stats_nslicesz);
L       = zeros(nstats,stats_nslicesz);
ReL     = zeros(nstats,stats_nslicesz);

% For each slice, loop through time:
for islice = 1 : stats_nslicesz
    
   tic
   % Declarations:
   istat=1;
   file = [basedir 'stat_slcz.' num2str(indexstat(istat),'%4.4d') '.' ...
                                num2str(islice,'%3.3d') '.bin']; 
        
   [data]=get_slice_z(file);   
   nptfftx=data.NX;
   nptffty=data.NY;    
   Puut_x = zeros(nstats,nptfftx);
   Pvvt_x = zeros(nstats,nptfftx);
   Pwwt_x = zeros(nstats,nptfftx);
   Puut_y = zeros(nstats,nptffty);
   Pvvt_y = zeros(nstats,nptffty);
   Pwwt_y = zeros(nstats,nptffty);
   
   
   for istat = 1 : nstats 
    
       
       file = [basedir 'stat_slcz.' num2str(indexstat(istat),'%4.4d') '.' ...
                                    num2str(islice,'%3.3d') '.bin']; 
       [data]=get_slice_z(file);
       
       % Mean velocities:
       uavg(istat,islice) = mean(mean(data.u));
       vavg(istat,islice) = mean(mean(data.v));
       wavg(istat,islice) = mean(mean(data.w));
       pavg(istat,islice) = mean(mean(data.p));
       
       % Reynolds stresses:
       uu(istat,islice) = mean(mean(data.u .^2)) - uavg(istat,islice)^2;
       vv(istat,islice) = mean(mean(data.v .^2)) - vavg(istat,islice)^2;
       ww(istat,islice) = mean(mean(data.w .^2)) - wavg(istat,islice)^2;
       
       uv(istat,islice) = mean(mean(data.u .* data.v)) - ...
                          uavg(istat,islice)*vavg(istat,islice);
       uw(istat,islice) = mean(mean(data.u .* data.w)) - ...
                          uavg(istat,islice)*wavg(istat,islice);
       vw(istat,islice) = mean(mean(data.v .* data.w)) - ...
                          vavg(istat,islice)*wavg(istat,islice);
                      
       % Urms ans TKE:               
       % Substract mean values on the slice:
       data.u = data.u - uavg(istat,islice);
       data.v = data.v - vavg(istat,islice);               
       data.w = data.w - wavg(istat,islice);
       
       % urms2:
       urms2(istat,islice)  = mean(mean(data.u .^2));
       vrms2(istat,islice)  = mean(mean(data.v .^2));
       wrms2(istat,islice)  = mean(mean(data.w .^2));
       
       % q2:
       q2(istat,islice)     = mean(mean(data.u .^2 + ...
                                        data.v .^2 + ...
                                        data.w .^2 ));
                                    
       % Longitudinal and transverse derivatives:
       dudx(istat,islice)   = mean(mean(data.dudx .^2 + ...
                                        data.dvdy .^2 + ...
                                        data.dwdz .^2 ))/3;
                                    
       dudy(istat,islice)   = mean(mean(data.dudy .^2 + ...
                                        data.dudz .^2 + ...
                                        data.dvdx .^2 + ...
                                        data.dvdz .^2 + ...
                                        data.dwdx .^2 + ...
                                        data.dwdy .^2 ))/6;
                                    
       % Microscales:
       lambdaf(istat,islice)= sqrt(2/3*q2(istat,islice)/dudx(istat,islice));
       lambdag(istat,islice)= sqrt(2/3*q2(istat,islice)/dudy(istat,islice)); 
       
       % Molecular and SGS dissipation:
       % Sij:
       S11 = data.dudx;
       S12 = 0.5*(data.dudy + data.dvdx);
       S13 = 0.5*(data.dudz + data.dwdx);
       S22 = data.dvdy;
       S23 = 0.5*(data.dvdz + data.dwdy);
       S33 = data.dwdz;
       
       % SGS stresses:
       T11_sgs = -2 * data.nut .* S11;
       T12_sgs = -2 * data.nut .* S12;
       T13_sgs = -2 * data.nut .* S13;
       T22_sgs = -2 * data.nut .* S22;
       T23_sgs = -2 * data.nut .* S23;
       T33_sgs = -2 * data.nut .* S33;
       
       % Average sgs dissipation TijSij:
       eps_sgs(istat,islice) = mean(mean( T11_sgs.*S11 + T22_sgs.*S22 + T33_sgs.*S33 + ...
                                       2*(T12_sgs.*S12 + T13_sgs.*S13 + T23_sgs.*S23 )));
                                   
       % Molecular dissipation:
       eps(istat,islice) = data.nu * ...
                           (3*dudx(istat,islice) + 6*dudy(istat,islice));

       % Reynolds Number based on transverse microscale:
       Relambda(istat,islice) = 1/3*(sqrt(urms2(istat,islice)) + ...
                                     sqrt(vrms2(istat,islice)) + ...
                                     sqrt(wrms2(istat,islice)))* ...
                                     lambdag(istat,islice)/data.nu;
                           
                           
       % Turbulent viscosity:
       nut(istat,islice) = mean(mean(data.nut));
       
       % Length scale, using total dissipation:
       L(istat,islice)   = (1/2*q2(istat,islice))^(3/2) / (eps(istat,islice)+eps_sgs(istat,islice));
       
       % Reynolds number based on the length scale:
       ReL(istat,islice) = (1/2*q2(istat,islice))^(1/2) * L(istat,islice) / data.nu;
       
       % Velocity spectra:
%        nptfftx=data.NX;
%        nptffty=data.NY;
       
       % ffts in the x direction:
       fftuo = zeros(nptffty,nptfftx);
       fftvo = zeros(nptffty,nptfftx);
       fftwo = zeros(nptffty,nptfftx);
       
       for jj=1:nptffty
           
          uof = zeros(1,nptfftx);
          vof = zeros(1,nptfftx);
          wof = zeros(1,nptfftx);
          
          for ii=1:nptfftx
             uof(ii) = data.u(ii,jj);
             vof(ii) = data.v(ii,jj);
             wof(ii) = data.w(ii,jj);
          end
          
          fftuo(jj,:) = fft(uof,nptfftx);
          fftvo(jj,:) = fft(vof,nptfftx);
          fftwo(jj,:) = fft(wof,nptfftx);
          
       end

       % Power Spectra:
       Puu_x = fftuo(:,1:nptfftx) .* conj(fftuo(:,1:nptfftx)) / nptfftx^2;
       Pvv_x = fftvo(:,1:nptfftx) .* conj(fftvo(:,1:nptfftx)) / nptfftx^2;
       Pww_x = fftwo(:,1:nptfftx) .* conj(fftwo(:,1:nptfftx)) / nptfftx^2;
       
%        Puut_x(istat,:) = zeros(1,nptfftx);
%        Pvvt_x(istat,:) = zeros(1,nptfftx);
%        Pwwt_x(istat,:) = zeros(1,nptfftx);
       for ii=1:nptfftx
           Puut_x(istat,ii) = mean(Puu_x(:,ii));
           Pvvt_x(istat,ii) = mean(Pvv_x(:,ii));           
           Pwwt_x(istat,ii) = mean(Pww_x(:,ii));
       end
       
       % ffts in the y direction:
       fftuo = zeros(nptfftx,nptffty);
       fftvo = zeros(nptfftx,nptffty);
       fftwo = zeros(nptfftx,nptffty);   
       
       for ii=1:nptfftx
           
          uof = zeros(1,nptffty);
          vof = zeros(1,nptffty);
          wof = zeros(1,nptffty);
          
          for jj=1:nptffty
             uof(jj) = data.u(ii,jj);
             vof(jj) = data.v(ii,jj);
             wof(jj) = data.w(ii,jj);
          end
          
          fftuo(ii,:) = fft(uof,nptffty);
          fftvo(ii,:) = fft(vof,nptffty);
          fftwo(ii,:) = fft(wof,nptffty);
          
       end
       
       % Power Spectra:
       Puu_y = fftuo(:,1:nptffty) .* conj(fftuo(:,1:nptffty)) / nptffty^2;
       Pvv_y = fftvo(:,1:nptffty) .* conj(fftvo(:,1:nptffty)) / nptffty^2;
       Pww_y = fftwo(:,1:nptffty) .* conj(fftwo(:,1:nptffty)) / nptffty^2;
       
%        Puut_y(istat,:) = zeros(1,nptffty);
%        Pvvt_y(istat,:) = zeros(1,nptffty);
%        Pwwt_y(istat,:) = zeros(1,nptffty);
       for jj=1:nptffty
           Puut_x(istat,jj) = mean(Puu_y(:,jj));
           Pvvt_x(istat,jj) = mean(Pvv_y(:,jj));           
           Pwwt_x(istat,jj) = mean(Pww_y(:,jj));
       end       
       
       
   end 
      
   % Mean statistics in time:
   DATAavg(islice).uavg = mean(uavg(:,islice));
   DATAavg(islice).vavg = mean(vavg(:,islice));
   DATAavg(islice).wavg = mean(wavg(:,islice));
   DATAavg(islice).pavg = mean(pavg(:,islice));
   
   DATAavg(islice).uu   = mean(uu(:,islice));
   DATAavg(islice).vv   = mean(vv(:,islice));
   DATAavg(islice).ww   = mean(ww(:,islice));
   DATAavg(islice).uv   = mean(uv(:,islice));
   DATAavg(islice).uw   = mean(uw(:,islice));
   DATAavg(islice).vw   = mean(vw(:,islice));
   
   DATAavg(islice).urms2= mean(urms2(:,islice));
   DATAavg(islice).vrms2= mean(vrms2(:,islice));
   DATAavg(islice).wrms2= mean(wrms2(:,islice));
   
   DATAavg(islice).q2   = mean(q2(:,islice));
   
   DATAavg(islice).dudx = mean(dudx(:,islice));
   DATAavg(islice).dudy = mean(dudy(:,islice));
   
   DATAavg(islice).lambdaf = mean(lambdaf(:,islice));
   DATAavg(islice).lambdag = mean(lambdag(:,islice));
   
   DATAavg(islice).nut     = mean(nut(:,islice));
   DATAavg(islice).eps_sgs = mean(eps_sgs(:,islice));
   DATAavg(islice).eps     = mean(eps(:,islice));
   
   DATAavg(islice).Relambda= mean(Relambda(:,islice));
   
   DATAavg(islice).L       = mean(L(:,islice));
   DATAavg(islice).ReL     = mean(ReL(:,islice));
   
   % Spectra in x direction
   nptfftx=data.NX;
   nptffty=data.NY;   

   Euu_x = zeros(1,nptfftx);
   Evv_x = zeros(1,nptfftx);
   Eww_x = zeros(1,nptfftx);   
   for ii=1:nptfftx
      Euu_x(ii)   = mean(Puut_x(:,ii));
      Evv_x(ii)   = mean(Pvvt_x(:,ii));
      Eww_x(ii)   = mean(Pwwt_x(:,ii));
   end   
   
   DATAavg(islice).kx = 2*pi*[0:floor(nptfftx/2)]/(data.del(IAXIS)*nptfftx);
   
   % Spectra in y direction
   Euu_y = zeros(1,nptffty);
   Evv_y = zeros(1,nptffty);
   Eww_y = zeros(1,nptffty);      
   for jj=1:nptffty
      Euu_y(jj)   = mean(Puut_y(:,jj)); 
      Evv_y(jj)   = mean(Pvvt_y(:,jj)); 
      Eww_y(jj)   = mean(Pwwt_y(:,jj)); 
   end
   
   DATAavg(islice).ky = 2*pi*[0:floor(nptffty/2)]/(data.del(JAXIS)*nptffty);
   
   % This only works if NX=NY:
   DATAavg(islice).Ecross = 0.25*(Euu_y + Eww_x + Evv_x + Eww_y);
   DATAavg(islice).Elong  =  0.5*(Euu_x + Evv_y);
   
   DATAavg(islice).NX      = data.NX;
   DATAavg(islice).NY      = data.NY;
   
   DATAavg(islice).zc      = data.zc; % all zc are the same
   DATAavg(islice).nu      = data.nu; % all nu are the same
   
   disp(['Slice=' num2str(islice) ', Stats done.'])
   
   toc
   
end


% %% Figure:
% stats_nslicesz =16;
% 
% secpause = 1;
% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
% for istat = stats_start : stats_end
%     
%     for islice = 1 : stats_nslicesz
% 
%         subplot(1,4,1)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).u','EdgeColor','flat')
%         title(['U veloc for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         %colorbar
%         view([0 90])
%         
%         subplot(1,4,2)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dudx','EdgeColor','flat')
%         title(['dUdx for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         %colorbar
%         view([0 90])        
%         
%         subplot(1,4,3)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dudy','EdgeColor','flat')
%         title(['dUdy for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         %colorbar
%         view([0 90])        
% 
%         subplot(1,4,4)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dudz','EdgeColor','flat')
%         title(['dUdz for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         %colorbar
%         view([0 90])        
%         
%         pause(secpause)
% 
%     end
% end
% 
% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
% for istat = stats_start : stats_end
%     
%     for islice = 1 : stats_nslicesz
% 
%         subplot(1,4,1)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).v','EdgeColor','flat')
%         title(['V veloc for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])
%         %colormap('gray')
%         %shading interp
%    
%         subplot(1,4,2)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dvdx','EdgeColor','flat')
%         title(['dVdx for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])        
%         
%         subplot(1,4,3)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dvdy','EdgeColor','flat')
%         title(['dVdy for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])        
% 
%         subplot(1,4,4)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dvdz','EdgeColor','flat')
%         title(['dVdz for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])              
%         
%         pause(secpause)       
%         
%     end
% end
% 
% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
% for istat = stats_start : stats_end
%     
%     for islice = 1 : stats_nslicesz
%         
%         subplot(1,4,1)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).w','EdgeColor','flat')
%         title(['W veloc for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])        
%         
%         subplot(1,4,2)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dwdx','EdgeColor','flat')
%         title(['dWdx for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])        
%         
%         subplot(1,4,3)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dwdy','EdgeColor','flat')
%         title(['dWdy for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])        
% 
%         subplot(1,4,4)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).dwdz','EdgeColor','flat')
%         title(['dWdz for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])             
%         
%         pause(secpause)
% 
%     end
% end
% 
% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
% for istat = stats_start : stats_end
%     
%     for islice = 1 : stats_nslicesz
%         
%         subplot(1,2,1)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).p','EdgeColor','flat')
%         title(['Pressure for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         %colorbar
%         axis equal
%         view([0 90])
%               
%         subplot(1,2,2)
%         surf(DATA(istat,islice).xcell, ...
%              DATA(istat,islice).ycell,DATA(istat,islice).nut','EdgeColor','flat')
%         title(['Tvis for istat, islice=' num2str(istat) ', ' ...
%                num2str(islice)],'FontSize',14)
%         xlabel('x','FontSize',14)
%         ylabel('y','FontSize',14)
%         set(gca,'FontSize',12)
%         grid on
%         box on
%         axis equal
%         view([0 90])        
%                        
%         pause(secpause)
% 
%     end
% end

return