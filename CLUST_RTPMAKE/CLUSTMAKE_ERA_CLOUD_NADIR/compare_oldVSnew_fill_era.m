addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/TIME
for hhx = 0:23
  fileOld = ['/asl/data/rtprod_airs/2002/09/01/sergio_nadir_cloudy_airs_l1b_era_sarta_baum_ice.2002.09.01.'           num2str(hhx,'%02d') '.rtp'];
  fileNew = ['/asl/data/rtprod_airs/2002/09/01/test_new_fillera_nadir_cloudy_airs_l1b_era_sarta_baum_ice.2002.09.01.' num2str(hhx,'%02d') '.rtp'];
  if exist(fileOld) & exist(fileNew)
    dirO = dir(fileOld);
    dirN = dir(fileNew);
    if dirO.bytes > 0 & dirN.bytes > 0
      [h0,ha0,p0,pa0] = rtpread(fileOld);
      [hN,haN,pN,paN] = rtpread(fileNew);
      if p0.rtime(1) > 1.4e9
        [yy,mm,dd,hh] = tai2utc(pN.rtime - offset1958_to_1993);
      else
        [yy,mm,dd,hh] = tai2utc(pN.rtime);
      end
      [sum(p0.rlon-pN.rlon) sum(p0.rlat-pN.rlat) sum(p0.rtime-pN.rtime)]
      figure(1); plot(hh,p0.stemp-pN.stemp); title(num2str(hhx,'%02d')); 
      figure(2); scatter_coast(p0.rlon,p0.rlat,20,p0.stemp-pN.stemp); title(num2str(hhx,'%02d'));
      disp('ret to continue'); pause
    end
  end
end
