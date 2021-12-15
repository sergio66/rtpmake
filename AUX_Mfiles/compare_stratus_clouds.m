all_waterprof = [];
all_watertop  = [];
all_wateramt  = [];
all_lat = [];
all_lon = [];
all_obs = [];
all_cal = [];

xall_waterprof = [];
xall_watertop  = [];
xall_wateramt  = [];
xall_lat = [];
xall_lon = [];
xall_obs = [];
xall_cal = [];

for mm = 3 : 3
  for dd = 10 : 10
    dir0 = ['//strowdataN/s3/robinso2/PUYEHUE/CLOUDS_6chans/TYPEIS_rCfracType=-1_x10/'];
    xthedir = dir([dir0 'eracld_type0_allfov_0.375cut_2011_03_10*.rtp']);
    xthedir = dir([dir0 'eracld_type0_allfov_0.425cut_2011_03_10*.rtp']);
    xthedir = dir([dir0 'eracld_type0_allfov_0.625cut_2011_03_10*.rtp']);
    thedir = dir([dir0 'eracld_type0_allfov_2011_03_10*.rtp']);

    for ii = 1 : length(thedir)
      figure(1); clf

      fname = [dir0 thedir(ii).name];
      [h,ha,p0,pa] = rtpread(fname);
      oo = find(p0.landfrac <= 0.01);
      if length(oo) > 0
        [h,p] = subset_rtp(h,p0,[],[],oo);
        p.clwc = p0.clwc(:,oo);
        p.ciwc = p0.ciwc(:,oo);
        p.cc = p0.cc(:,oo);

        figure(1);
        pcolor(1:length(p.stemp),p.plevs,double(log10(p.clwc))); 
        shading flat; colorbar; set(gca,'ydir','reverse'); hold all

        watertop = nan(size(p.stemp));  profswater = nan(size(p.clwc));  wateramt = watertop; 
        obs = watertop; cal = watertop;
        lons = nan(size(p.stemp)); lats = nan(size(p.stemp)); 

        water = find(p.ctype == 101);  
        watertop(water) = p.cprtop(water);  profswater(:,water) = p.clwc(:,water);
        i1 = length(water);  lons(water) = p.rlon(water); lats(water) = p.rlat(water);
        wateramt(water) = p.cngwat(water);
        obs(water) = p.robs1(4,water);
        cal(water) = p.rcalc(4,water);

        water = find(p.ctype2 == 101);
        watertop(water) = p.cprtop2(water);  profswater(:,water) = p.clwc(:,water);
        i2 = length(water);  lons(water) = p.rlon(water); lats(water) = p.rlat(water);
        wateramt(water) = p.cngwat2(water);
        obs(water) = p.robs1(4,water);
        cal(water) = p.rcalc(4,water);

        plot(watertop,'ko','linewidth',3); hold off
        title([num2str(mm) '/' num2str(dd) ' : ' num2str(ii) ' : ' num2str(i1) ' and ' num2str(i2)]); 

        all_waterprof = [all_waterprof profswater];
        all_watertop  = [all_watertop  watertop];
        all_wateramt  = [all_wateramt  wateramt];
        all_lat       = [all_lat lats]     ;
        all_lon       = [all_lon lons]     ;
        all_obs       = [all_obs obs]     ;
        all_cal       = [all_cal cal]     ;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%
      fname = [dir0 xthedir(ii).name];
      [h,ha,p0,pa] = rtpread(fname);
      oo = find(p0.landfrac <= 0.01);
      if length(oo) > 0
        [h,p] = subset_rtp(h,p0,[],[],oo);
        p.clwc = p0.clwc(:,oo);
        p.ciwc = p0.ciwc(:,oo);
        p.cc = p0.cc(:,oo);

        figure(2);
        pcolor(1:length(p.stemp),p.plevs,double(log10(p.clwc))); 
        shading flat; colorbar; set(gca,'ydir','reverse'); hold all

        watertop = nan(size(p.stemp));  profswater = nan(size(p.clwc)); wateramt = watertop;
        obs = watertop; cal = watertop;
        lons = nan(size(p.stemp)); lats = nan(size(p.stemp)); 

        water = find(p.ctype == 101);  
        watertop(water) = p.cprtop(water);  profswater(:,water) = p.clwc(:,water);
        i1 = length(water);  lons(water) = p.rlon(water); lats(water) = p.rlat(water);
        wateramt(water) = p.cngwat(water);
        obs(water) = p.robs1(4,water);
        cal(water) = p.rcalc(4,water);

        water = find(p.ctype2 == 101);
        watertop(water) = p.cprtop2(water);  profswater(:,water) = p.clwc(:,water);
        i2 = length(water);  lons(water) = p.rlon(water); lats(water) = p.rlat(water);
        wateramt(water) = p.cngwat2(water);
        obs(water) = p.robs1(4,water);
        cal(water) = p.rcalc(4,water);

        plot(watertop,'ko','linewidth',3); hold off
        title([num2str(mm) '/' num2str(dd) ' : ' num2str(ii) ' : ' num2str(i1) ' and ' num2str(i2)]); 

        xall_waterprof = [xall_waterprof profswater];
        xall_watertop  = [xall_watertop  watertop];
        xall_wateramt  = [xall_wateramt  wateramt];
        xall_lat       = [xall_lat lats]     ;
        xall_lon       = [xall_lon lons]     ;
        xall_obs       = [xall_obs obs]     ;
        xall_cal       = [xall_cal cal]     ;
      end

      ret(1)
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(1); clf

%  pcolor(1:length(all_lat),p.plevs(:,1),double(log10(all_waterprof))); 
%  shading flat; colorbar; set(gca,'ydir','reverse'); hold all
%  plot(all_watertop,'k','linewidth',3); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dn = -100 : 1 : 100; 
tall_obs = rad2bt(1231,all_obs);    tall_cal = rad2bt(1231,all_cal); 
tall_obs(tall_obs < 150) = NaN;     tall_cal(tall_cal < 150) = NaN; 
xtall_obs = rad2bt(1231,xall_obs);  xtall_cal = rad2bt(1231,xall_cal); 
xtall_obs(xtall_obs < 150) = NaN;   xtall_cal(xtall_cal < 150) = NaN; 

dn = -150 : 1: 150; 
  nnBIAS = nanhist(tall_obs - tall_cal,dn);
  nnxBIAS = nanhist(xtall_obs - xtall_cal,dn);
semilogy(dn,nnBIAS,'b',dn,nnxBIAS,'r')
  title('(b) bias OLD  (r) bias NEW')

plot_stratus_compare
