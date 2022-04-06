%{


Hi Sergio,

Thanks for your patience on the AMIP runs. You can find the
model-means on Taki at /asl/s1/rkramer/CMIP6_modelmean_Sergio. Like
last time you’ll see a separate file for air temperature “ta”, surface
air temperature “tas” atmospheric specific humidity “hus” and for
surface specific humidity “huss”. You are looking for the versions
with “amip” in the filename.  Note the AMIP runs are only for Jan 1979
to December 2014, shorter than the fully coupled historical runs. So
keep that in mind in case your codes need some modifying.

I used same 11 models we used for the historical, hist-nat and
hist-GHG model-mean results so you can easily compare with AMIP. At
some point, especially if you want to publish these results, it may be
worth adding more models in for AMIP and historical runs. We stopped
at these 11 models because they are the only ones also available for
hist-nat and hist-GHG but I’m not sure you are planning to use those
extra two experiments anymore (used to decompose the historical
results into contributions from natural vs GHG forcing only). I don’t
think adding models will change your model-mean results, but it can’t
hurt to add more if they are available. At least for the sake of
completeness.

Let me know if you run into any issues! I saved the original amip
files for each model in that directory tree we created at
/asl/models/cmip6 in case you want to use them down the road.

 

Ryan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/asl/models/cmip6/[experiment name]/[variable type]/[temporal frequency]/[variable name]/[model name]/*files.nc
For example for the file tas_Amon_NorESM2-LM_historical_r1i1p1f1_gn_201001-201412.nc,  the path would be: 
/asl/models/cmip6/historical/atmos/mon/tas/NorESM2-LM

Included are netcdfs for surface-air temperature (tas), atmospheric
temperature (ta), surface specific humidity (huss), and atmospheric
specific humidity (hus) for the full historical run with all forcings
(historical), the historical run with only natural forcings
(hist-nat), and the historical run with only GHG forcing (hist-GHG).
The filenames should explain which is which but let me know if
anything is confusing. Essentially there is just one variable in each
file, but it includes the entire timeseries from 1850-2015. These
model-means were calculated from the dozen or so models listed in the
readme located under /asl/models/cmip6.

I set the lat and lon as the midpoint between the edges you gave me,
so the files are 64 x 72.  I haven’t really looked into the files too
closely so if you see anything funky when you do your analysis, let
me know. Now that I have all the code set up, it will be really quick
to make changes.

dir0 = '/asl/models/cmip6/';
dir0 = '/asl/s1/rkramer/CMIP6_modelmean_Sergio/';

expt_name = 'hist-GHG';   %% only GHG forcings
expt_name = 'hist-nat';   %% only natural forcings
expt_name = 'historical'; %% all forcings

var_type = 'atmos';       %% atmos only, no ocean

temporal = 'monthly';     %% monthly output

var_name{1} = 'huss';        %% surface     specific humidity
var_name{2} = 'hus';         %% atmospheric specific humidity
var_name{3} = 'tas';         %% surface     air temp
var_name{4} = 'ta';          %% atmospheric air temp

model_name = 'HadGEM3-GC31-LL';
model_name = 'GISS-E2-1-G';

for ii = 1 : 4
  thedir0 = [dir0 '/' expt_name '/' var_type '/' temporal '/' var_name{ii} '/' model_name '/'];
  x = dir([thedir0 '*200101-201412.nc']);
  fname = [thedir0 '/' x.name];
  s = read_netcdf_lls(fname);
  error('opo')
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/h4tools
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity
addpath /home/sergio/MATLABCODE/TROPOPAUSE

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 1

fip = mktempS('fx.ip.rtp');
fop = mktempS('fx.op.rtp');
frp = mktempS('fx.rp.rtp');

hus = read_netcdf_lls('/asl/s1/rkramer/CMIP6_modelmean_Sergio/hus_CMIP6_amip_modelmean.nc');
huss = read_netcdf_lls('/asl/s1/rkramer/CMIP6_modelmean_Sergio/huss_CMIP6_amip_modelmean.nc');
ta   = read_netcdf_lls('/asl/s1/rkramer/CMIP6_modelmean_Sergio/ta_CMIP6_amip_modelmean.nc');
tas  = read_netcdf_lls('/asl/s1/rkramer/CMIP6_modelmean_Sergio/tas_CMIP6_amip_modelmean.nc');

%% "time" is hours since noon at the midpoint of the first month

figure(1); pcolor(nanmean(tas.tas,3)'); colormap jet; shading flat; colorbar; title('mean ST')
figure(2); pcolor(nanmean(huss.huss,3)'); colormap jet; shading flat; colorbar; title('mean SH at surface')

allyears = [];
allmonths = [];

amip6_start = 1850;
amip6_start = 1979;

numtime = (2015-amip6_start)*12;   %% yes, the 1980 timesteps Ryan has saved

for ii = amip6_start:2014
  junkyear = ii*ones(1,12);
  junkmonth = [1 2 3 4 5 6 7 8 9 10 11 12];
  allyears = [allyears junkyear];
  allmonths = [allmonths junkmonth];
end

i2002 = (2002-amip6_start)*12 + 9;  %% 2002/09
i2014 = (2014-amip6_start)*12 + 8;  %% 2014/09
[allyears(i2002) allmonths(i2002) allyears(i2014) allmonths(i2014)]

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

if length(pa) == 0
  pa = {{'profiles','rtime','seconds since 1993'}};
  pa = {{'profiles','rtime','seconds since 1958'}};
end
if length(ha) == 0
  ha = {{'header','hdf file','amip6 stuff'}};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyuse = [2002 2002 2002 2002];
mmuse = [09   10   11   12  ];
for yy = 2003:2013
  yjunk = yy*ones(1,12);
  mjunk = [1 2 3 4 5 6 7 8 9 10 11 12];
  yyuse = [yyuse yjunk];
  mmuse = [mmuse mjunk];
end
yy = 2014; yjunk = yy*ones(1,8);
           mjunk = [1 2 3 4 5 6 7 8];
  yyuse = [yyuse yjunk];
  mmuse = [mmuse mjunk];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yySE = [2002 2014];  %% should do this
yyuse = ones(1,4)*yySE(1); 
mmuse = [09   10   11   12  ];

%% now go from 2003 to 2013
for yy = yySE(1)+1 : yySE(2)-1
  yjunk = yy*ones(1,12);
  mjunk = [1 2 3 4 5 6 7 8 9 10 11 12];
  yyuse = [yyuse yjunk];
  mmuse = [mmuse mjunk];
end
yy = yySE(2); yjunk = yy*ones(1,8);
           mjunk = [1 2 3 4 5 6 7 8];
  yyuse = [yyuse yjunk];
  mmuse = [mmuse mjunk];
[1:length(yyuse); yyuse; mmuse]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for ii = i2002 : i2014
for ii = JOB + (i2002-1)
  [allyears(ii) allmonths(ii) yyuse(ii-i2002+1) mmuse(ii-i2002+1)]  
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/AMIP6/Tile_Center/amip6_tile_center_monthly_timestep_' num2str(ii-i2002+1,'%03d') '.mat']; 
  
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/time
  p.rlon = wrapTo180(p.rlon);
%  [p,pa] = rtp_add_emis(p,pa);

  pnew_ip = p;
  hnew_ip = h;
  
  hnew_ip = rmfield(hnew_ip,'ngas');    hnew_ip.ngas = 1;
  hnew_ip = rmfield(hnew_ip,'glist');   hnew_ip.glist = 1;
  hnew_ip = rmfield(hnew_ip,'gunit');   hnew_ip.gunit = 21;
                                        hnew_ip.ptype = 0;
                                        hnew_ip.pfields = 3;
 
  %pnew_ip.solzen = solzen;
  %pnew_ip.satzen = satzen;
  %pnew_ip.scanang = saconv(p.satzen,p.zobs);  
  %pnew_ip.rtime   = rtime;

  pnew_ip = rmfield(pnew_ip,'stemp');
  pnew_ip = rmfield(pnew_ip,'ptemp');
  pnew_ip = rmfield(pnew_ip,'plevs');
  pnew_ip = rmfield(pnew_ip,'palts');
  pnew_ip = rmfield(pnew_ip,'mmw');
  pnew_ip = rmfield(pnew_ip,'gas_1');
  pnew_ip = rmfield(pnew_ip,'gas_2');
  pnew_ip = rmfield(pnew_ip,'gas_3');
  pnew_ip = rmfield(pnew_ip,'gas_4');
  pnew_ip = rmfield(pnew_ip,'gas_5');
  pnew_ip = rmfield(pnew_ip,'gas_6');
  pnew_ip = rmfield(pnew_ip,'gas_9');
  pnew_ip = rmfield(pnew_ip,'gas_12');

  pnew_ip = rmfield(pnew_ip,'cprtop');
  pnew_ip = rmfield(pnew_ip,'cprbot');
  pnew_ip = rmfield(pnew_ip,'cfrac');
  pnew_ip = rmfield(pnew_ip,'cngwat');
  pnew_ip = rmfield(pnew_ip,'cpsize');
  pnew_ip = rmfield(pnew_ip,'ctype');
  pnew_ip = rmfield(pnew_ip,'cprtop2');
  pnew_ip = rmfield(pnew_ip,'cprbot2');
  pnew_ip = rmfield(pnew_ip,'cfrac2');
  pnew_ip = rmfield(pnew_ip,'cngwat2');
  pnew_ip = rmfield(pnew_ip,'cpsize2');
  pnew_ip = rmfield(pnew_ip,'ctype2');
  pnew_ip = rmfield(pnew_ip,'cfrac12');

  pnew_ip = rmfield(pnew_ip,'rcalc');
  pnew_ip = rmfield(pnew_ip,'sarta_rclearcalc');
  
  junk = tas.tas(:,:,ii); junk = junk(:);                 pnew_ip.stemp = junk';
  junk = ta.ta(:,:,:,ii); junk = permute(junk,[3 1 2]);   junk = reshape(junk,19,72*64); pnew_ip.ptemp = junk;
  junk = hus.hus(:,:,:,ii); junk = permute(junk,[3 1 2]); junk = reshape(junk,19,72*64); pnew_ip.gas_1 = junk;
  junk = hus.plev * ones(1,4608);                         pnew_ip.plevs = junk/100;  %% hPa to mb
  junk = 19 * ones(1,4608);                               pnew_ip.nlevs = junk;
  pnew_ip.rtime = utc2taiSergio(yyuse(ii-i2002+1),mmuse(ii-i2002+1),15,12.0) * ones(1,4608);

  [xyy,xmm,xdd,xhh] = tai2utcSergio(pnew_ip.rtime);        % <<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
  time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
  co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
  pnew_ip.co2ppm = co2ppm;
  fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),pnew_ip.co2ppm(1));
  fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),pnew_ip.co2ppm(end));
  figure(1); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.co2ppm);
  figure(2); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp); title('stemp')
  pause(0.1)

  run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  %run_sarta.klayers_code = '/home/sergio/KLAYERS/BinV201/klayers_airs_wetwater_140levs';
  run_sarta.clear = +1;
  run_sarta.cloud = +1;
  run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC
  run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
  code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  code1 = '/home/sergio/SARTA_CLOUDY/BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  code1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
  run_sarta.sartaclear_code = code1;
  run_sarta.sartacloud_code = code1;
  run_sarta.co2ppm = co2ppm;

  rtpwrite(fip,hnew_ip,ha,pnew_ip,pa);  
  klayerser = ['!' run_sarta.klayers_code '      fin=' fip ' fout=' fop ' >& ugh']; eval(klayerser);
  sartaer   = ['!' run_sarta.sartacloud_code '   fin=' fop ' fout=' frp ' >& ugh']; eval(sartaer);
  [hnew_op,ha2,pnew_op,pa2] = rtpread(frp);
  %pnew_op.rcalc            = p2.rcalc;
  %pnew_op.sarta_rclearcalc = p2.sarta_rclearcalc;
  pnew_op.mmw = mmwater_rtp(hnew_op,pnew_op);

  %[hnew_ip,ha,pnew_ip,pa] = rtptrim_sartacloud(hnew_ip,ha,p2,pa);
  rmer = ['!/bin/rm '  fip ' ' fop ' ' frp]; eval(rmer);

  [Tw,Tw1km,Tdew,WBGT,RH,RH1km,colwater,TwSurf,RHSurf,TdewSurf] = layeramt2RH_wet_bulb_dew_point_temperature(hnew_op,pnew_op);
  %if isfield(pnew_ip,'d2m') & isfield(pnew_ip,'t2m')
  %  pnew_op.d2m = pnew_ip.d2m;
  %  pnew_op.t2m = pnew_ip.t2m;
  %  pnew_op.rh2m = airtemp_dewpointtemp_2_RH(pnew_op.t2m,pnew_op.d2m);
  %end
  pnew_op.Tw   = Tw;
  pnew_op.Tdew = Tdew;
  pnew_op.RH   = RH;
  pnew_op.mmw  = colwater;
  pnew_op.TwSurf   = TwSurf;
  pnew_op.TdewSurf = TdewSurf;
  pnew_op.RHSurf   = RHSurf;

  pnew_ip.rh = convert_humidity (pnew_ip.plevs*100,pnew_ip.ptemp,pnew_ip.gas_1,'specific humidity','relative humidity');

  yyuseII = yyuse(ii-i2002+1);
  mmuseII = mmuse(ii-i2002+1);
  comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compute_amip6_profile_rtpfiles.m';
  comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_CMIP6/clust_compute_amip6_profile_rtpfiles.m';
  saver = ['save ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2    yyuseII mmuseII'];
  if ~exist(fout)
    [yyuseII mmuseII]
    fprintf(1,'saving %s \n',fout)     
    eval(saver)
  else
    fprintf(1,'%s already exists ... not saving \n',fout) 
  end

  figure(1); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,pnew_op.stemp); colormap jet; title('stemp (K)');
  figure(2); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,pnew_op.mmw); colormap jet; title('col water (mm)');
  figure(3); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,pnew_op.spres); colormap jet; title('spres (mb)');
  figure(4); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,rad2bt(1231,pnew_op.rcalc(1520,:))); colormap jet; title('BT1231 (K)');

  i1419 = find(hnew_ip.vchan >= 1419,1); 
  figure(4); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,rad2bt(1419,pnew_op.rcalc(i1419,:))); colormap jet; title('BT1419 (K)');

  figure(3); plot(hnew_ip.vchan,std(rad2bt(hnew_ip.vchan,pnew_op.rcalc)'))
  figure(3); plot(hnew_ip.vchan,mean(rad2bt(hnew_ip.vchan,pnew_op.rcalc)'))
  pause(0.1)

end

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp];
eval(rmer);
