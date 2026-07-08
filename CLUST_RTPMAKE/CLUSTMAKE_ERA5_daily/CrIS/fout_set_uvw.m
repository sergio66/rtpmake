%% airs_l1c ---> cris_fsr

%{
if iInterp <= 0  
  if iERAorECM == 1
    fmatout = ['fsr_allfov_era5_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['uvw_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['uvw_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fmatout = ['fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['uvw_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['uvw_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
else
  if iERAorECM == 1
    fmatout = ['interp_analysis_fsr_allfov_era5_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['interp_analysis_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['interp_analysis_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fmatout = ['interp_analysis_fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['interp_analysis_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['interp_analysis_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fout_set

if iInterp <= 0  
  if iERAorECM == 1
    fmatout = ['fsr_allfov_era5_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['uvw_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['uvw_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  elseif iERAorECM == -1
    fmatout = ['fsr_allfov_ecm_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['uvw_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['uvw_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  end
else
  if iERAorECM == 1
    fmatout = ['interp_analysis_fsr_allfov_era5_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_cris_fsr_era5_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  elseif iERAorECM == -1
    fmatout = ['interp_analysis_fsr_allfov_ecm_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_cris_fsr_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  end
end
