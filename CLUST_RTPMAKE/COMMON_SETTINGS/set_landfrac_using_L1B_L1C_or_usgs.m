iUSGS = +1;  %% swtich using usgs      but we lost the database!
iUSGS = +2;  %% swtich using ETOPO2022 landfrac is good but salti for lakes and eg bahamas, S. Pacific is bad
iUSGS = +3;  %% swtich using GDEMM2024 landfrac is good and salti seems good
iUSGS = -1;  %% stick with L1B/L1C

if iUSGS == 1
  p.landfrac_fromL1B = p.landfrac;
  p.salti_fromL1B = p.salti;
  [salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
  p.landfrac = landfrac;
  p.salti    = salti;
elseif iUSGS == 2
  addpath /home/sergio/git/matlabcode/DEM_DigitalELeveationModel
  p.landfrac_fromL1B = p.landfrac;
  p.salti_fromL1B = p.salti;
  [salti,landfrac,gebco] = etopo_dem_and_imerg_lf(p.rlat,p.rlon);
  p.landfrac = landfrac;
  p.salti    = salti;
elseif iUSGS == 3
  addpath /home/sergio/git/matlabcode/DEM_DigitalELeveationModel
  p.landfrac_fromL1B = p.landfrac;
  p.salti_fromL1B = p.salti;
  [salti,landfrac,gebco] = gdemm_dem_and_imerg_lf(p.rlat,p.rlon);
  p.landfrac = landfrac;
  p.salti    = salti;
else
  disp('set_landfrac_using_L1B_L1C_or_usgs.m : we have lost usgs_dem and ETOPO2022 has issues in open ocean/inland lakes')
  disp('set_landfrac_using_L1B_L1C_or_usgs.m : just use what is in radiance files')
  if ~isfield(p,'salti')
    error('set_landfrac_using_L1B_L1C_or_usgs.m : p has no salti field')
  end
  if ~isfield(p,'landfrac')
    error('set_landfrac_using_L1B_L1C_or_usgs.m : p has no landfrac field')
  end
end
