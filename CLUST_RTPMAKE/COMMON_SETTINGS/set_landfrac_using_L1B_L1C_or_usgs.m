iUSGS = +1;  %% swtich, but we lost the database!
iUSGS = -1;  %% stick with L1B/L1C

if iUSGS > 0
  p.landfrac_fromL1B = p.landfrac;
  p.salti_fromL1B = p.salti;
  [salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
  p.landfrac = landfrac;
  p.salti    = salti;
end
