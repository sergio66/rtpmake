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

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlab2012/airs/readers
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

%addpath /home/strow/cress/Work/Rtp
%addpath /home/strow/Matlab/Grib

%addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
%addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm 
%addpath /asl/rtp_prod2/grib/                           %% for fill_ecmwf

%addpath  /asl/packages/rtp_prod2/grib
%addpath  /home/sbuczko1/git/rtp_prod2/grib

addpath /home/strow/git/rtp_prod2/grib
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB
addpath /home/sergio/MATLABCODE/TIME

theinds = (1 : 2378)';

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

%icestr = ['NEWLANDFRAC/cloudy_airs_l1b_ecm' icestr '.'];
icestr = ['cloudy_airs_l1b_ecm' icestr '.'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%this is input file name
rtpdir = '/asl/data/rtprod_airs/';
rtpdir = '/asl/rtp/rtprod_airs/';
rtpdir = '/home/sergio/MATLABCODE/SONDES/MAGIC/Sergio2016/';

file0   = [rtpdir 'all_era_match_sonde.rp.rtp'];
fileOUT = [rtpdir 'all_ecm_match_sonde.rp.rtp'];

fin0  = file0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% typically iaGList == JOB so this loop runs ONCE .. so you control muktiple files via the cluster

[h,ha,p,pa] = rtpread(fin0);
hORIG = h;
pORIG = p;

p = rmfield(p,'stemp');
p = rmfield(p,'ptemp');
p = rmfield(p,'plevs');
p = rmfield(p,'palts');
p = rmfield(p,'nlevs');
p = rmfield(p,'tcc');
p = rmfield(p,'cc');
p = rmfield(p,'clwc');
p = rmfield(p,'ciwc');

for ii = 1 : h.ngas
  str = ['gas_'  num2str(h.glist(ii))];
  str = ['p = rmfield(p,''' str ''');'];
  eval(str)
end

h = rmfield(h,'ngas');
h = rmfield(h,'glist');
h = rmfield(h,'gunit');
h.ngas = 0;

[yy,mm,dd] = tai2utcSergio(p.rtime);

eeP = exist(fileOUT);

if eeP == 0
  fprintf(1,' making %s \n',fileOUT);
  toucher = ['!touch ' fileOUT];
  eval(toucher)

  which fill_ecmwf
  disp('calling fill_ecmwf')

  %%% EASY PEASY
  %%% [p,h] = fill_ecmwf(p,h);

  %% COMPLICATED
  iCnt = 0;
  for yyx = 2012 : 2013
    for mmx = 1 : 12
      if mmx ~= 12
        woo = find(yy == yyx & mm == mmx);
      else
        woo = find(yy == yyx & mm == mmx & dd < 31);
      end
      if length(woo) > 0
        iCnt = iCnt + 1;
        [hx,px] = subset_rtp(h,p,[],[],woo);
        [px,hx] = fill_ecmwf(px,hx);
        fprintf(1,'%4i/%2i : length(woo), length(stemp) = %4i %4i \n',yyx,mmx,length(woo),length(px.stemp))
	saveYY(iCnt) = yyx;
	saveMM(iCnt) = mmx;	
	saveWOO(iCnt) = length(woo);
	saveLEN(iCnt) = length(px.stemp);
	if iCnt == 1
	  hall = hx;
	  pall = px;
	else
	  [hall,pall] = cat_rtp(hall,pall,hx,px);
	end
      end
    end
  end

  save replace_fields_with_ecmwf.mat hall pall
  
  p = pall;
  h = hall;
  
  p0 = p;

  %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
  %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
  %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
  %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
  addpath /asl/rtp_prod2/emis/
  addpath /asl/rtp_prod2/util/
  [p,pa] = rtp_add_emis(p,pa);
  %figure(1)
  %scatter_coast(p.rlon,p.rlat,10,p.nemis); 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  run_sarta.ForceNewSlabs = +1;
  [p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

  fnamex = fileOUT;
  [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
  rtpwrite(fnamex,h,ha,p2x,pa)
else
  fprintf(1,' %s already exists \n',fileOUT)
end
