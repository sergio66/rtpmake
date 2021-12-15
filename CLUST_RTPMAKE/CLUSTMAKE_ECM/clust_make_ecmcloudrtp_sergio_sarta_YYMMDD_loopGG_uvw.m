%% have to hardcode  YYMMDD before starting, loops over granules
%% gets in few extra vars, such as vertical velocity, divergence, potential vorticity

%% see http://www.ecmwf.int/products/catalogue/I.html

% local running to test
% clustcmd -L clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_uvw.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 -p 4 clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_uvw.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 -p 4 clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_uvw.m 001:240

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2012 09 20];   %% for SNO with CrIS

yymmdd0  = [2011 01 11];   %% for JGR paper
yymmdd0  = [2011 06 12];   %% for JGR paper <<< 2011/06/11 is BAD
yymmdd0  = [2011 07 11];   %% for JGR paper

yymmdd0  = [2011 03 11];   %% for JGR paper

iaGlist  = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE

theinds = (1 : 2378)';
theinds = 1291;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = 9999;

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

if codeX == 0
  icestr = '_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end

icestruvw = ['uvw_pv_cloudy_airs_l1b_ecm' icestr '.'];
icestr    = [       'cloudy_airs_l1b_ecm' icestr '.'];

%%%%%%%%%%%%%%%%%%%%%%%%%

for ixx = 1 : length(iaGlist)
  ix = iaGlist(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  fnameOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fnameOUT = [fnameOUT icestruvw ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  xfnameOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  xfnameOUT = [xfnameOUT icestr ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  eeP = exist(fnameOUT);
  eeX = exist(xfnameOUT);

  if eeP == 0 & eeX > 0
    fprintf(1,' updating %s \n',xfnameOUT);

    [hhx,hhax,ppx, ppax] = rtpread(xfnameOUT);
  
    year  = yymmddgg(1);
    month = yymmddgg(2);
    day   = yymmddgg(3);
    gran  = yymmddgg(4);
    if mod(year,4) == 0
      mos = [31 29 31 30 31 30 31 31 30 31 30 31];  %% leap year
    else
      mos = [31 28 31 30 31 30 31 31 30 31 30 31];  %% normal year
      end
    days_so_far = 0;
    if month > 1
      days_so_far = sum(mos(1:month-1));
    end
    days_so_far = days_so_far + day;
    filename = ['/strowdataN/data/airs/Aqua_AIRS_Level1/AIRIBRAD.005/' ystr '/'];
    filename = ['/asl/data/airs/AIRIBRAD/' ystr '/'];
    filename = [filename num2str(days_so_far,'%03d') '/'];
    dir0 = filename;
    filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
    filename = [filename '.L1B.AIRS_Rad.v5*.hdf'];

    thedir = dir(filename);
    if length(thedir) == 1
      fname = [dir0 thedir.name];
    else
      fprintf(1,'%s \n',filename);
      disp('file does not exist');
      return
    end

    [a,b,c] = sdload_quiet(fname);
    p.rlat = [a.Latitude(:)'];
    p.rlon = [a.Longitude(:)'];
    p.rtime = [a.Time(:)'];
  
    [meantime, f, prof] = readl1b_all(fname);  
    p = prof;

    p.pobs = zeros(size(p.solazi));
    p.upwell = ones(size(p.solazi));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','hdf file',filename}};

    h.pfields=5; % (1=prof + 4=IRobs);

    h.nchan = length(theinds);
    h.ichan = theinds;;
    h.vchan = f(h.ichan);;

pXX = p;

    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};
    uvwfields = {'SKT','W'};
    uvwfields = {'U','V','W','PV','D','SP'};
    [h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,uvwfields); %%% add on ecm

    u_all = p.grib_U;
    v_all = p.grib_V;
    w_all = p.grib_W;

    p = rmfield(p,'grib_U');
    p = rmfield(p,'grib_V');
    p = rmfield(p,'grib_W');

    %% wz([1 2 3],:) = velocities at 250,500,850 mb

    for ii = 1 : length(p.rlon)
      pp = p.plevs(:,ii);

      if max(pp) < 300
        oo(ii) = length(pp);
        uz(1,ii) = NaN;
        vz(1,ii) = NaN;
        wz(1,ii) = NaN;
        fprintf(1,'hmmm 300 mb looks like index %5i could be high altitude! \n',ii);
      else
        dadiff = abs(pp-300);
        oo(ii) = find(dadiff == min(dadiff),1);
        uz(1,ii) = u_all(oo(ii),ii);
        vz(1,ii) = v_all(oo(ii),ii);
        wz(1,ii) = w_all(oo(ii),ii);
      end   

      if max(pp) < 500
        oo(ii) = length(pp);
        uz(2,ii) = NaN;
        vz(2,ii) = NaN;
        wz(2,ii) = NaN;
        fprintf(1,'hmmm 500 mb looks like index %5i could be high altitude! \n',ii);
      else
        dadiff = abs(pp-500);
        oo(ii) = find(dadiff == min(dadiff),1);
        uz(2,ii) = u_all(oo(ii),ii);
        vz(2,ii) = v_all(oo(ii),ii);
        wz(2,ii) = w_all(oo(ii),ii);
      end   

      if max(pp) < 850
        oo(ii) = length(pp);
        uz(3,ii) = NaN;
        vz(3,ii) = NaN;
        wz(3,ii) = NaN;
        fprintf(1,'hmmm 850 mb looks like index %5i could be high altitude! \n',ii);
      else
        dadiff = abs(pp-850);
        oo(ii) = find(dadiff == min(dadiff),1);
        uz(3,ii) = u_all(oo(ii),ii);
        vz(3,ii) = v_all(oo(ii),ii);
        wz(3,ii) = w_all(oo(ii),ii);
      end   

    end
    p.uz = uz;
    p.vz = vz;
    p.wz = wz;

    p0 = p;

    figure(1); scatter_coast(p.rlon,p.rlat,25,p.uz(2,:)); title('u at 500 mb')
    figure(2); scatter_coast(p.rlon,p.rlat,25,p.vz(2,:)); title('v at 500 mb')
    figure(3); scatter_coast(p.rlon,p.rlat,25,p.wz(2,:)); title('w at 500 mb')
    speed = sqrt(p.uz(2,:).^2 + p.vz(2,:).^2);
    figure(4); scatter_coast(p.rlon,p.rlat,25,speed); title('u/v speed at 500 mb')
    ind = 1 : 100 : 12150;
    hold on; 
      quiver(p.rlon(ind),p.rlat(ind),p.uz(2,ind),p.vz(2,ind),'color','k','linewidth',3); 
    hold off

    figure(1); scatter_coast(p.rlon,p.rlat,25,p.uz(1,:)); title('u at 300 mb')
    figure(2); scatter_coast(p.rlon,p.rlat,25,p.vz(1,:)); title('v at 300 mb')
    figure(3); scatter_coast(p.rlon,p.rlat,25,p.wz(1,:)); title('w at 300 mb')
    speed = sqrt(p.uz(2,:).^2 + p.vz(2,:).^2);
    figure(4); scatter_coast(p.rlon,p.rlat,25,speed); title('u/v speed at 300 mb')
    ind = 1 : 100 : 12150;
    hold on; 
      quiver(p.rlon(ind),p.rlat(ind),p.uz(1,ind),p.vz(1,ind),'color','k','linewidth',3); 
    hold off
   
    if sum(ppx.rtime - p.rtime) < eps
      ppx.uz = p.uz;
      ppx.vz = p.vz;
      ppx.wz = p.wz;
    else
      error('rtime DOES NOT MATCH')
    end
%    tobs = rad2bt(1231,p.robs1(1291,:));
%    tcld = rad2bt(1231,ppx.rcalc(1291,:));
%    plot(tobs,p.wz,'.',tcld,p.wz,'r.')

    %rtpwrite(fnamex,h,ha,p2x,pa)
    rtpwrite(xfnameOUT,hhx,hhax,ppx,ppax);

  end
end
