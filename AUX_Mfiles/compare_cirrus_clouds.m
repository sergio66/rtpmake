all_iceprof = [];
all_icetop  = [];
all_iceamt  = [];
all_lat = [];
all_lon = [];

xall_iceprof = [];
xall_icetop  = [];
xall_iceamt  = [];
xall_lat = [];
xall_lon = [];

for mm = 3 : 3
  for dd = 10 : 10
    dir0 = ['//strowdataN/s3/robinso2/PUYEHUE/CLOUDS_6chans/TYPEIS_rCfracType=-1_x10/'];
    xthedir = dir([dir0 'eracld_type0_allfov_0.375cut_2011_03_10*.rtp']);
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
        pcolor(1:length(p.stemp),p.plevs,double(log10(p.ciwc))); 
        shading flat; colorbar; set(gca,'ydir','reverse'); hold all

        icetop = nan(size(p.stemp));  profsice = nan(size(p.ciwc));  iceamt = icetop;
        lons = nan(size(p.stemp)); lats = nan(size(p.stemp)); 

        ice = find(p.ctype == 201);  
        icetop(ice) = p.cprtop(ice);  profsice(:,ice) = p.ciwc(:,ice);
        i1 = length(ice);  lons(ice) = p.rlon(ice); lats(ice) = p.rlat(ice);
        iceamt(ice) = p.cngwat(ice);

        ice = find(p.ctype2 == 201);
        icetop(ice) = p.cprtop2(ice);  profsice(:,ice) = p.ciwc(:,ice);
        i2 = length(ice);  lons(ice) = p.rlon(ice); lats(ice) = p.rlat(ice);
        iceamt(ice) = p.cngwat2(ice);

        plot(icetop,'ko','linewidth',3); hold off
        title([num2str(mm) '/' num2str(dd) ' : ' num2str(ii) ' : ' num2str(i1) ' and ' num2str(i2)]); 

        all_iceprof = [all_iceprof profsice];
        all_icetop  = [all_icetop  icetop];
        all_iceamt  = [all_iceamt  iceamt];
        all_lat       = [all_lat lats]     ;
        all_lon       = [all_lon lons]     ;
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
        pcolor(1:length(p.stemp),p.plevs,double(log10(p.ciwc))); 
        shading flat; colorbar; set(gca,'ydir','reverse'); hold all

        icetop = nan(size(p.stemp));  profsice = nan(size(p.ciwc)); iceamt = icetop;
        lons = nan(size(p.stemp)); lats = nan(size(p.stemp)); 

        ice = find(p.ctype == 201);  
        icetop(ice) = p.cprtop(ice);  profsice(:,ice) = p.ciwc(:,ice);
        i1 = length(ice);  lons(ice) = p.rlon(ice); lats(ice) = p.rlat(ice);
        iceamt(ice) = p.cngwat(ice);

        ice = find(p.ctype2 == 201);
        icetop(ice) = p.cprtop2(ice);  profsice(:,ice) = p.ciwc(:,ice);
        i2 = length(ice);  lons(ice) = p.rlon(ice); lats(ice) = p.rlat(ice);
        iceamt(ice) = p.cngwat2(ice);

        plot(icetop,'ko','linewidth',3); hold off
        title([num2str(mm) '/' num2str(dd) ' : ' num2str(ii) ' : ' num2str(i1) ' and ' num2str(i2)]); 

        xall_iceprof = [xall_iceprof profsice];
        xall_icetop  = [xall_icetop  icetop];
        xall_iceamt  = [xall_iceamt  iceamt];
        xall_lat       = [xall_lat lats]     ;
        xall_lon       = [xall_lon lons]     ;
      end
      ret(1)
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(1); clf

%  pcolor(1:length(all_lat),p.plevs(:,1),double(log10(all_iceprof))); 
%  shading flat; colorbar; set(gca,'ydir','reverse'); hold all
%  plot(all_icetop,'k','linewidth',3); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_cirrus_compare
