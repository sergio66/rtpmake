%{
if iInterp <= 0  
  if iERAorECM == 1
    fmatout = ['fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['uvw_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['uvw_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fmatout = ['fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['uvw_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['uvw_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
else
  if iERAorECM == 1
    fmatout = ['interp_analysis_fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fmatout = ['interp_analysis_fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fmatout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fmatout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iInterp <= 0  
  if iERAorECM == 1
    fmatout = ['fsr_allfov_era_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['uvw_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['uvw_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  elseif iERAorECM == -1
    fmatout = ['fsr_allfov_ecm_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['uvw_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['uvw_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  end
else
  if iERAorECM == 1
    fmatout = ['interp_analysis_fsr_allfov_era_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  elseif iERAorECM == -1
    fmatout = ['interp_analysis_fsr_allfov_ecm_' num2str(gg,'%03d') '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.mat'];
    fmatout = ['interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.mat'];
  end
end
