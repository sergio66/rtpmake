%% testing ~/MATLABCODE/matlib/clouds/sarta/driver_sarta_cloud_rtp.m

prof = profXYZ;
  prof = reset_cprtop_cloudOD(prof,run_sarta.cumsum/100,airslevels,airsheights);
  plot_hists_cprtop_cprtop2

disp('---> checking cprtop vs cprbot vs spres')
iNotOK = +1;
iFix = 0;
%% before used to give up at iFix == 10
while iFix < 12 & iNotOK > 0
  iFix = iFix + 1;
  [prof,iNotOK] = check_for_errors(prof,run_sarta.cfrac,iFix);  %% see possible pitfalls in clouds
    plot_hists_cprtop_cprtop2;
    disp('ret to continue'); pause
  fprintf(1,' did n=%2i try at checking clouds \n',iFix)
end
if iFix >= 12 & iNotOK > 0
  %disp('oops, could not fix cprtop vs cprbot vs spres'); %keyboard
  error('oops, could not fix cprtop vs cprbot vs spres')
end
