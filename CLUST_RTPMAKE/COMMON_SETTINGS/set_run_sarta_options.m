run_sarta.clear = +1;
run_sarta.cloud = +1;

if iSlabCld_CumSumStrowORGeorge > 0
  run_sarta.cumsum = 9999;  %% strow pick, cloud at PEAK of wgt fcn
else
  run_sarta.cumsum = -1;  %% aumann pick, cloud at wgt mean of profile
end
run_sarta.klayers_code = klayers;

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
code1 = sartaCld;

if codeX == 0
  icestr = '_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end

run_sarta.sartaclear_code = run_sarta.sartacloud_code;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
