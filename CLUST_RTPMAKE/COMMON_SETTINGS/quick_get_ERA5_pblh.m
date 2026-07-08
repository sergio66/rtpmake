function  prof = quick_get_ERA5_pblh(profin,iTimeInterpolate);

if nargin == 1
  iTimeInterpolate = -1;  %% first time bracket UGH
  iTimeInterpolate =  0;  %% nearest in time
  iTimeInterpolate = +1;  %% interpolate across two 3-hour intervals
end

%% see eg /umbc/rs/pi_sergio/WorkDirDec2025/rtpmake/CLUST_RTPMAKE/CLUSTMAKE_ERA_or_ECM_CRIS_FSR/clustbatch_make_eraORecm_cloudrtp_sergio_sarta_YYMMDD_loopGG.m

prof = profin;

iAdd_pblh_nwp = -1;

[xyyyy,xmmmm,xdddd,xhhhh] = tai2utcSergio(nanmean(profin.rtime));
e_ymd = [num2str(xyyyy) num2str(xmmmm,'%02d') num2str(xdddd,'%02d')];
daysINyear = [31 28 31 30 31 30 31 31 30 31 30 31];
if mod(xyyyy,4) == 0
  daysINyear(2) = 29;
end  

ecm_8steps = 0:3:21;

fERA5A = ['/home/sergio/asl/isilonX/ERA5/' e_ymd(1:4) '/' e_ymd(5:6) '/era5_sfc_' e_ymd(1:4) e_ymd(5:6) e_ymd(7:8) '.nc'];
fERA5A = ['/umbc/rs/strow/asl/ERA5/' e_ymd(1:4) '/' e_ymd(5:6) '/era5_sfc_' e_ymd(1:4) e_ymd(5:6) e_ymd(7:8) '.nc'];

iHH5A = find(ecm_8steps <= xhhhh);
iHH5A = iHH5A(end);
iHH5A = min(iHH5A,8); %% only have 8 timesteps in the daily ERA5 file
wgtA = (xhhhh - ecm_8steps(iHH5A))/3;
wgtA = 1 - wgtA;   %% notice I do this

if wgtA > 1 | wgtA < 0
  wgtA
  error('quick_get_ERA5_pblh.m : wgtA outside 0,1')
end  

if xhhhh <= 21
  fERA5B = fERA5A;
  iHH5B = iHH5A + 1;
else
  meantime = profin.rtime;
  b_ymd = e_ymd;
  if xdddd < daysINyear(xmmmm)
    b_ymd = [num2str(xyyyy) num2str(xmmmm,'%02d') num2str(xdddd+1,'%02d')];
  elseif  xdddd == daysINyear(xmmmm) & xmmmm < 12
    b_ymd = [num2str(xyyyy) num2str(xmmmm+1,'%02d') num2str(1,'%02d')];
  elseif  xdddd == daysINyear(xmmmm) & xmmmm == 12
    b_ymd = [num2str(xyyyy+1) num2str(1,'%02d') num2str(1,'%02d')];
  end
  fERA5B = ['/home/sergio/asl/isilonX/ERA5/' b_ymd(1:4) '/' b_ymd(5:6) '/era5_sfc_' b_ymd(1:4) b_ymd(5:6) b_ymd(7:8) '.nc'];
  fERA5B = ['/umbc/rs/strow/asl/ERA5/' b_ymd(1:4) '/' b_ymd(5:6) '/era5_sfc_' b_ymd(1:4) b_ymd(5:6) b_ymd(7:8) '.nc'];
  iHH5B = 1;  
end
wgtB = 1 - wgtA;   %% notice I do this

fprintf(1,'mean hour of observation and brackets = %2i < %.2f < %2i    fERA5 files straddling this time for PBLH are \n %s and %s \n',ecm_8steps(iHH5A),xhhhh,ecm_8steps(iHH5B),fERA5A,fERA5B);
fprintf(1,'   iHH5A (ecm_8stepsA), wgtA = %2i %2i %.2f     iHH5B (ecm_8stepsB), wgtB = %2i %2i %.2f \n',iHH5A,ecm_8steps(iHH5A),wgtA,iHH5B,ecm_8steps(iHH5B),wgtB)

if iTimeInterpolate == 0
  %% use nearest in time
  if wgtA >= wgtB
    wgtA = 1.0;
    wgtB = 0.0;
  else
    wgtA = 0.0;
    wgtB = 1.0;
  end
elseif iTimeInterpolate < 0
  %% use first time bracket
  wgtA = 1.0;
  wgtB = 0.0;  
end

if exist(fERA5A) & exist(fERA5B)
  
  %% profin.rtime is in seconds; typical granule is 6 minutes long = 360 seconds
  %% each chunk here is 3 hours long = 180 minutes
  if (max(profin.rtime) - min(profin.rtime))/60 < 180
    iAdd_pblh_nwp = +1;      
    F5A = grib_interpolateERA5sfc(fERA5A,iHH5A);
    pblh_nwpA = F5A.pblh_nwp.ig(profin.rlat,profin.rlon)/1000; %%% in km
    scatter_coast(profin.rlon,profin.rlat,50,pblh_nwpA);   title('PBLH A  NWP (km)'); caxis([0 4])

    F5B = grib_interpolateERA5sfc(fERA5B,iHH5B);
    pblh_nwpB = F5B.pblh_nwp.ig(profin.rlat,profin.rlon)/1000; %%% in km
    scatter_coast(profin.rlon,profin.rlat,50,pblh_nwpB);   title('PBLH B NWP (km)'); caxis([0 4])

    pblh_nwp = pblh_nwpB * wgtB + pblh_nwpA * wgtA;
    %{
    addpath ~/git/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/
    ax = axis;
    read_ERA5_PBLH_daily_grib3
    figure(10); 
      simplemap(flipud(squeeze(pblhDataPBH(:,:,iHH5)))/1000);  colormap jet; colorbar; title('ERA5 PBH (km) MEAN'); caxis([0 4])
      axis(ax);
    %}
    
  end
end

if iAdd_pblh_nwp > 0
  prof.pblh_nwp = pblh_nwp;
  scatter_coast(profin.rlon,profin.rlat,50,pblh_nwp);   title('PBLH FINAL WGTA+WGTB NWP (km)'); caxis([0 4])  
end  
