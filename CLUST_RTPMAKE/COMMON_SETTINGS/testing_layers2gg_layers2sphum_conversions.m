%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
check
% Your layer boundaries in Pa (convert from mb)

p_bounds = [1050, 1010, 980, 960] * 100;   % Pa

% dp for each layer
dp = diff(p_bounds);   % [40, 30, 20] * 100 Pa

So your round-trip is straightforward:
% Forward: ECMWF q → molecules/cm2 for each layer
% q_layer is the layer-averaged specific humidity [kg/kg]
% dp is the layer pressure thickness [Pa]

N = mixratio_to_molecules(q_layer, dp);   % molecules/cm2

% Inverse: back to mixing ratio
w_recovered = molecules_to_mixratio(N, dp);

% These should be identical to machine precision
fprintf('Max round-trip error: %.2e kg/kg\n', max(abs(w_recovered - q_layer ./ (1 - q_layer))));

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath0

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 49;   %% daytime Australia, getting veddy veddy high PBLH
  JOB = 48;   %% daytime Australia, getting veddy veddy high PBLH  
  JOB = 200;  %% california  
  JOB = 213;  %% S.Africa to Antartica
  JOB = 219;  %% E. Pacific, with very interesting SKT dipole
  JOB = 210;  %% India
end

%% this is TESTING molecules/cm2 --> gg (mix ratio) --> sp. hum. kg/kg [units of ERA5, ECMWF, MERRA2]

gran = JOB;

ecm_rtpIN  = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.rtp'];
%% retrieval starting with ECM
fileECM    = ['/home/sergio/nogit/sergio_temp_rtp_files/singlefootprintretrievals_ccast_hires_jpss1/2024/11/13/'];
  fileECM  = [fileECM 'interp_analysis_ecm_retr' num2str(gran,'%03d') '_cris_-1_iDET_4_iStemp_ColWV_21_iCenterFov_-1_iCO2_Yes_No_Switch_-1_singlelayerclouds.mat'];

[h,ha,xpdmat,pa] = rtpread(ecm_rtpIN);

figure(1); clf; scatter_coast(xpdmat.rlon,xpdmat.rlat,10,xpdmat.stemp); title('SKT')

% specific humidity [kg/kg] --> mixing ratio [kg/kg]
%   w_mmr = q_specific_hum ./ (1 - q_specific_hum);
% mixing ratio [kg/kg] --> specific humidity [kg/kg]
%   q_specific_hum = w_mmr./ (1 + w_mmr);

xpdmat.gg = xpdmat.gas_1 ./ (1 - xpdmat.gas_1);

loader = ['load ' fileECM];
eval(loader);

pecm = poemNew;
pecm.plays = plevs2plays(pecm.plevs);
pecm.gas_1 = poemNew.gas_1_OEMinitialization;  %% in molecules/cm2
pecm.ptemp = poemNew.ptemp_OEMinitialization;  %% in K
pecm.stemp = poemNew.stemp_OEMinitialization;  %% in K
pecm.mmw = mmwater_rtp(hoemNew,pecm);
figure(2); clf; scatter_coast(xpdmat.rlon,xpdmat.rlat,10,pecm.mmw); title('mmw')

%%%%%%%%%%%%%%%%%%%%%%%%%

[ggLAY,ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2gg(hoemNew,pecm,1:length(pecm.stemp),1,1);   %% this is my orig code 

%%%%%%%%%%%%%%%%%%%%%%%%%

dpjunk = zeros(101,length(poemNew.stemp));
  dpjunk(1:100,:) = diff(poemNew.plevs,1)*100;  %% change mb to Pa
  [ggLAY1,qqLAY1] = molecules_to_mixratio(pecm.gas_1,abs(dpjunk),poemNew.nlevs-1);  %% gg is mass mix ratio while qq is sp. humidity

dLjunk = zeros(101,length(poemNew.stemp));
  dLjunk(1:100,:) = abs(diff(poemNew.palts,1));      %% in meteres
  [ggLAY2,qqLAY2] = recover_q_from_forward(pecm.gas_1,poemNew.plays*100,poemNew.ptemp,abs(dLjunk),poemNew.nlevs-1);  %% gg is mass mix ratio while qq is sp. humidity

%%%%%%%%%%%%%%%%%%%%%%%%%

iMethod = 1;
iMethod = 2;
if iMethod == 1
  ggLAYX = ggLAY1;
  qqLAYX = qqLAY1;
else
  ggLAYX = ggLAY2;
  qqLAYX = qqLAY2;
end
[mm,nn] = size(ggLAY);

%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(pecm.nlevs == 92,1);
ind = 1;

mmx = poemNew.nlevs(ind)-1;

figure(3); clf

%% xpdmat.gas_1 is straight from ECM so is specific humidity q in g/g           w = q/(1-q)
loglog(ggLAY(1:mmx,ind),pavgLAY(1:mmx,ind),'bx-',ggLAYX(1:mmx,ind),pavgLAY(1:mmx,ind),'r');
  legend('gg ORIG code','gg NEW code'); title('mass mix ratio')
  set(gca,'ydir','reverse');
  ylim([100 1010]);
  ylim([800 1010]);  
  
loglog(ggLAY(1:mmx,ind),pavgLAY(1:mmx,ind),'bx-',ggLAYX(1:mmx,ind),pavgLAY(1:mmx,ind),'r',xpdmat.gg(:,ind),xpdmat.plevs(:,ind),'g*-');;
  legend('gg ORIG code','gg NEW code','ecm'); title('mass mix ratio')
  set(gca,'ydir','reverse');
  ylim([100 1010])
  ylim([800 1010])  
  
%loglog(ggLAY(1:mmx,ind),pavgLAY(1:mmx,ind),'b',ggLAYX(1:mmx,ind),pavgLAY(1:mmx,ind),'r',...
%       xpdmat.gas_1(:,ind)./(1-xpdmat.gas_1(:,ind)),xpdmat.plevs(:,ind),'g');; set(gca,'ydir','reverse'); ylim([100 1010])

%%%%%%%%%%%%%%%%%%%%%%%%%

ocean = 1:length(xpdmat.landfrac);
ocean = find(xpdmat.landfrac == 1);
ocean = find(xpdmat.landfrac == 0);

whos ocean

qqLAY = ggLAY./(1+ggLAY);

%%% since xpdmat.plevs is all over the place, you really cant average
% loglog(nanmean(ggLAY(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'bx-',nanmean(ggLAYX(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'ro-',...
%        nanmean(xpdmat.gas_1(:,ocean),2)./(1-nanmean(xpdmat.gas_1(:,ocean),2)),nanmean(xpdmat.plevs(:,ocean),2),'gs-');;
%   set(gca,'ydir','reverse');
%   ylim([800 1010])
% loglog(nanmean(ggLAY(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'bx-',nanmean(ggLAYX(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'ro-',...
%        nanmean(xpdmat.gg(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'gs-');;
%   set(gca,'ydir','reverse');
%   xtitle('Mass Mix ratio')
%   ylim([800 1010])
% 
% loglog(nanmean(qqLAY(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'bx-',nanmean(qqLAYX(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'ro-',...
%   nanmean(xpdmat.gas_1(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'gs-');;
%   set(gca,'ydir','reverse');
%   xtitle('Specific humidity')
%   ylim([800 1010])

%%%%%%%%%%%%%%%%%%%%%%%%%

xpdmat.interped_gg = nan(101,length(pecm.stemp));
xpdmat.interped_sh = nan(101,length(pecm.stemp));
for ii = 1 : length(pecm.stemp)
  gaa = xpdmat.gas_1(:,ii);
  paa = xpdmat.plevs(:,ii);
  N   = pecm.nlevs(ii)-1;
  Paa = pecm.plays(1:N,ii);
  xpdmat.interped_sh(1:N,ii) = interp1(log(paa),gaa,log(Paa),[],'extrap');
end  
xpdmat.interped_gg = xpdmat.interped_sh./(1 + xpdmat.interped_sh);

loglog(nanmean(ggLAY(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'bx-',nanmean(ggLAYX(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'ro-',...
  nanmean(xpdmat.interped_gg(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'gs-');;
  set(gca,'ydir','reverse');
  xlabel('Mass Mix Ratio')
  ylim([800 1010])

loglog(nanmean(qqLAY(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'bx-',nanmean(qqLAYX(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'ro-',...
  nanmean(xpdmat.interped_sh(1:mm,ocean),2),nanmean(pavgLAY(1:mm,ocean),2),'gs-');;
  set(gca,'ydir','reverse');
  xlabel('Specific humidity')
  ylim([800 1010])

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
plot(nanmean(qqLAY(1:mm,ocean)./xpdmat.interped_sh(1:mm,ocean),2) - 1 ,nanmean(pavgLAY(1:mm,ocean),2),'b',...
     nanmean(qqLAYX(1:mm,ocean)./xpdmat.interped_sh(1:mm,ocean),2) - 1,nanmean(pavgLAY(1:mm,ocean),2),'r',...
     nanstd(qqLAY(1:mm,ocean)./xpdmat.interped_sh(1:mm,ocean),[],2)   ,nanmean(pavgLAY(1:mm,ocean),2),'c--',...
     nanstd(qqLAYX(1:mm,ocean)./xpdmat.interped_sh(1:mm,ocean),[],2)  ,nanmean(pavgLAY(1:mm,ocean),2),'m--','linewidth',2)
  set(gca,'ydir','reverse');
   plotaxis2;
  xlabel('Sp Humidity : layers->levels/ECM TRUE'); legend('My code','Claude')
  ylim([700 1010])

figure(4)
plot(nanmean(ggLAY(1:mm,ocean)./xpdmat.interped_gg(1:mm,ocean),2) - 1 ,nanmean(pavgLAY(1:mm,ocean),2),'b',...
     nanmean(ggLAYX(1:mm,ocean)./xpdmat.interped_gg(1:mm,ocean),2) - 1,nanmean(pavgLAY(1:mm,ocean),2),'r',...
     nanstd(ggLAY(1:mm,ocean)./xpdmat.interped_gg(1:mm,ocean),[],2)   ,nanmean(pavgLAY(1:mm,ocean),2),'c--',...
     nanstd(ggLAYX(1:mm,ocean)./xpdmat.interped_gg(1:mm,ocean),[],2)  ,nanmean(pavgLAY(1:mm,ocean),2),'m--','linewidth',2)
  set(gca,'ydir','reverse');
   plotaxis2;
  xlabel('Mass Mix Ratio : layers->levels/ECM TRUE'); legend('My code','Claude')
  ylim([700 1010])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% print out stats
pecm.dL = nan(101,length(pecm.stemp));
pecm.dL(1:100,:) = abs(diff(pecm.palts,1));
pecm.dP = nan(101,length(pecm.stemp));
pecm.dP(1:100,:) = abs(diff(pecm.plevs,1));
pecm.Ntotal = pecm.plays*100 .* pecm.dL/8.31./pecm.ptemp * 6.023e23/1e4;

nbot = mean(pecm.nlevs(ocean)) - 1;
nbot = nbot - 2;

fprintf('Q_H2O    [mol/cm2]: %.6e\n', mean(pecm.gas_1(nbot,ocean)))
fprintf('N_total  [mol/cm2]: %.6e\n', mean(pecm.Ntotal(nbot,ocean)))
fprintf('X        [mol/mol]: %.6e\n', mean(pecm.gas_1(nbot,ocean))/mean(pecm.Ntotal(nbot,ocean)));
fprintf('P_avg    [Pa]:      %.2f\n', mean(pecm.plays(nbot,ocean))*100);
fprintf('T_avg    [K]:       %.2f\n', mean(pecm.ptemp(nbot,ocean)));
fprintf('L        [m]:       %.2f\n', mean(pecm.dL(nbot,ocean)));
fprintf('dp       [Pa]:      %.2f\n', mean(pecm.dP(nbot,ocean))*100);
fprintf('q_ecmwf  [kg/kg]:   %.6f\n', mean(xpdmat.interped_sh(nbot,ocean)));
fprintf('q_recov  [kg/kg]:   %.6f\n', mean(qqLAYX(nbot,ocean)));
fprintf('ratio:              %.4f\n', mean(qqLAYX(nbot,ocean))/mean(xpdmat.interped_sh(nbot,ocean)))

for iL = 1:100   % sample 100 ocean profiles
  % fprintf('Q_H2O    [mol/cm2]: %.6e\n', mean(pecm.gas_1(iL,ocean)))
  % fprintf('N_total  [mol/cm2]: %.6e\n', mean(pecm.Ntotal(iL,ocean)))
  % fprintf('X        [mol/mol]: %.6e\n', mean(pecm.gas_1(iL,ocean))/mean(pecm.Ntotal(iL,ocean)));
  % fprintf('P_avg    [Pa]:      %.2f\n', mean(pecm.plays(iL,ocean))*100);
  % fprintf('T_avg    [K]:       %.2f\n', mean(pecm.ptemp(iL,ocean)));
  % fprintf('L        [m]:       %.2f\n', mean(pecm.dL(iL,ocean)));
  % fprintf('dp       [Pa]:      %.2f\n', mean(pecm.dP(iL,ocean))*100);
  % fprintf('q_ecmwf  [kg/kg]:   %.6f\n', mean(xpdmat.interped_sh(iL,ocean)));
  % fprintf('q_recov  [kg/kg]:   %.6f\n', mean(qqLAYX(iL,ocean)));
  % fprintf('ratio:              %.4f\n', mean(qqLAYX(iL,ocean))/mean(xpdmat.interped_sh(iL,ocean)))
  ratio(iL)  = mean(qqLAYX(iL,ocean))/mean(xpdmat.interped_sh(iL,ocean));
  ratio2(iL) = mean(qqLAYX(iL,ocean)./xpdmat.interped_sh(iL,ocean));  
end
figure(5)
semilogy(ratio2(1:mm),nanmean(pavgLAY(1:mm,ocean),2))
set(gca,'ydir','reverse'); ylim([700 1010])

error('kshjskgjhskgjhsg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wawoo1 = nanmean(ggLAY(1:mm,ocean),2);
wawoo2 = nanmean(ggLAYX(1:mm,ocean),2);
wawooP = nanmean(pavgLAY(1:mm,ocean),2);
wawoo1 = interp1(log(wawooP),wawoo1,log(nanmean(xpdmat.plevs,2)),[],'extrap');
wawoo2 = interp1(log(wawooP),wawoo2,log(nanmean(xpdmat.plevs,2)),[],'extrap');

wawooE = nanmean(xpdmat.gas_1,2)./(1+nanmean(xpdmat.gas_1,2));
wawooPE = nanmean(xpdmat.plevs,2);

semilogy(wawoo1./wawooE - 1, wawooPE, wawoo2./wawooE - 1, wawooPE); set(gca,'ydir','reverse'); ylim([100 1010])
ylim([800 1020])

%%%%%%%%%%%%%%%%%%%%%%%%%

