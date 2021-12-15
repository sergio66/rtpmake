function [p1,p2] = convert_era2ecm(file1,file2)

fprintf(1,'in = %s out = %s \n',file1,file2)

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ECM
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
%% sarta = '/home/sergio/SARTA_CLOUDY//BinV201/sarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_newdirs';
sarta   = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_pge_v6_tunmlt';
sarta   = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

if exist(file2)
  iY = input('your output file exists!!!! proceed (+1/-1) Y/N : ');
  if iY < 0
    error('stoppping')
  end
end

[h,ha,p1,pa] = rtpread(file1);
if h.ptype == 1
  error('need h.ptype == 0')
end

iF = 0;
nwp = {'gas_1','gas_2','gas_3','gas_4','gas_5','gas_6','gas_9','gas_12','ptemp','stemp','cc','tcc','clwc','ciwc','nlevs','plevs'};
slab = {'cngwat','cngwat2','cprtop','cprtop2','cprbot','cprbot2','cpsize','cpsize','cfrac','cfrac2','cfrac12','ctype','ctype2'};
junk = {'orig_ctop','orig_ctop2','sarta_lvlODice','sarta_lvlODwater','sarta_rclearcalc'};

list = nwp;
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p1,xfield)
    iF = iF + 1;
    p1 = rmfield(p1,xfield);
 end
end
list = slab;
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p1,xfield)
    iF = iF + 1;
    p1 = rmfield(p1,xfield);
 end
end
list = junk;
for ii = 1 : length(list)
  xfield = list{ii};
  if isfield(p1,xfield)
    iF = iF + 1;
    p1 = rmfield(p1,xfield);
 end
end
fprintf(1,'removed %3i fields \n',iF);

[h,ha,p2,pa] = quickmake_oneECM_file(h,ha,p1,pa);
rtpwrite(file2,h,ha,p2,pa);