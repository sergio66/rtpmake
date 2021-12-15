function [h,ha,p2x,pa] = quickmake_oneERA_file(h,ha,p,pa)

%[h,ha,p,pa] = rtpadd_era_data(h,ha,p,pa,cldfields); %%% add on era
[p,h] = fill_era(p,h);
    
p0 = p;

%[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
%p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
%p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
%[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);

addpath /asl/rtp_prod2/emis/
addpath /asl/rtp_prod2/util/
%[p,pa] = rtp_add_emis(p,pa);
    
%figure(1)
%scatter_coast(p.rlon,p.rlat,10,p.nemis); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params
if codeX == 0
  icestr = '_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end

run_sarta.ForceNewSlabs = +1;  %% force new slabs
[p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

[h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);

