%{
if iInterp <= 0  
  if iERAorECM == 1
    fout = ['fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fout = ['cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fout = ['fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fout = ['cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
else
  if iERAorECM == 1
    fout = ['interp_analysis_fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fout = ['interp_analysis_fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iInterp <= 0  
  if iERAorECM == 1
    fout = ['fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fout = ['cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fout = ['fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fout = ['cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
else
  if iERAorECM == 1
    fout = ['interp_analysis_fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fout = ['interp_analysis_fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
end
