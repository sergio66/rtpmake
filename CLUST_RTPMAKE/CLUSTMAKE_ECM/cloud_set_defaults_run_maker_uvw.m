%{
https://www.mathworks.com/help/matlab/matlab_prog/suppress-warnings.html?requestedDomain=www.mathworks.com

Warning: HDFSD will be removed in a future release. Use MATLAB.IO.HDF4.SD instead.
> In hdfsd (line 261)
In sdload (line 116)
In cloud_set_defaults_run_maker (line 121)
In clustbatch_make_ecmcloudrtp_sergio_sarta_filelist (line 32)

>> w = warning('query','last')
w = identifier: 'MATLAB:imagesci:hdf:removalWarningHDFSD'
         state: 'on'
>> id = w.identifier;
>> warning('off',id)
>> lastwarn
ans = 
HDFSW will be removed in a future release. Use MATLAB.IO.HDFEOS.SW instead.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iv5or6 = 5;   %% AIRS L1B
iv5or6 = 6;   %% AIRS L1C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

add_the_paths_and_klayers_sarta_execs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iv5or6 == 5
  theinds = (1 : 2378)';
  theinds = 1291;  
else
  theinds = (1 : 2645)';
  theinds = 1520;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_run_sarta_options

% this is set in calling routine cloud_set_defaults_run_makerBLAH.m

%% icestr = ['NEWLANDFRAC/cloudy_airs_l1b_ecm' icestr '.'];
if iv5or6 == 5
  icestruvw = ['uvw_cloudy_airs_l1b_ecm' icestr '.'];
  icestr    = [    'cloudy_airs_l1b_ecm' icestr '.'];  
elseif iv5or6 == 6
  icestruvw = ['uvw_cloudy_airs_l1c_ecm' icestr '.'];
  icestr    = [    'cloudy_airs_l1c_ecm' icestr '.'];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ixx = 1 : length(iaGlist)
  ix = iaGlist(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  dout_set

  if ~exist(fdirOUT)
    mker = ['!mkdir -p ' fdirOUT];
    eval(mker);
    fprintf(1,'made %s \n',fdirOUT)
  end
  
  if iSlabCld_CumSumStrowORGeorge == 1
    fnameOUT    = [fdirOUT icestr    ystr '.' mstr '.' dstr '.' gstr '.rtp'];
    uvwfnameOUT = [fdirOUT icestruvw ystr '.' mstr '.' dstr '.' gstr '.mat'];    
  else
    fnameOUT    = [fdirOUT icestr    ystr '.' mstr '.' dstr '.' gstr '_cumsum_-1.rtp'];
    uvwfnameOUT = [fdirOUT icestruvw ystr '.' mstr '.' dstr '.' gstr '_cumsum_-1.mat'];    
  end

  eeP = exist(fnameOUT);
  eeX = exist(uvwfnameOUT);

  if eeX == 0 & eeP >= 0
    fprintf(1,' updating %s using %s \n',uvwfnameOUT,fnameOUT);
    if eeP > 0
      [hhx,hhax,ppx, ppax] = rtpread(fnameOUT);
    end
    
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

    read_in_L1B_or_L1C

    if iv5or6 == 5
      p.rlat = [a.Latitude(:)'];
      p.rlon = [a.Longitude(:)'];
      p.rtime = [a.Time(:)'];

      %[meantime, f, prof] = readl1b_all(fname);  %% has the old rtime 1993
      [meantime, f, prof] = readl1b_all(fname);  %% has the new rtime 1958    
      p = prof;
    elseif iv5or6 == 6
      p = gdata;
    end
    
    p.pobs = zeros(size(p.solazi));
    p.upwell = ones(size(p.solazi));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    plot(p.rlon,p.rlat,'.')

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
    
    %     [h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,uvwfields); %%% add on ecm
    %     u_all = p.grib_U;
    %     v_all = p.grib_V;
    %     w_all = p.grib_W;
    %     p = rmfield(p,'grib_U');
    %     p = rmfield(p,'grib_V');
    %     p = rmfield(p,'grib_W');

    %% wz([1 2 3],:) = velocities at 250,500,850 mb

    [p,h] = fill_ecmwf(p,h,[],1);
    p0 = p;

%%%%%%%%%%%%%%%%%%%%%%%%%

    [xyy,xmm,xdd,xhh] = tai2utcSergio(p.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
    time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
    co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
    p.co2ppm = co2ppm;
    run_sarta.co2ppm = p.co2ppm;
    fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),p.co2ppm(1));
    fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),p.co2ppm(end));

    p0 = p;

    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);

    p.rlon = wrapTo180(p.rlon);    
    add_the_DanZhou_emis

    %figure(1)
    %scatter_coast(p.rlon,p.rlat,10,p.nemis); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

    fnamex = fnameOUT;
    fnamex = uvwfnameOUT;

    %[h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    if ~exist(fnamex)
      %rtpwrite(fnamex,h,ha,p2x,pa)
      [xhd0,xpdmat] = get_richardson_number_levels(h,ha,p2,pa);      
      saver = ['save ' fnamex ' xhd0 xpdmat '];
      eval(saver)
      fprintf(1,'saved %s \n',fnamex)
      %% plot_richardson_PBLH      
    else
      fprintf(1,'%s already exists, not saving \n',fnamex)
    end

%    tobs = rad2bt(1231,p.robs1(1291,:));
%    tcld = rad2bt(1231,ppx.rcalc(1291,:));
%    plot(tobs,p.wz,'.',tcld,p.wz,'r.')

    %rtpwrite(fnamex,h,ha,p2x,pa)
    %rtpwrite(xfnameOUT,hhx,hhax,ppx,ppax);

  end
end
