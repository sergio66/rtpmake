function [] = esrl_mipas_ch4_to_klayers(rdate,opts)

% esrl_mipas_ch4_to_klayers.m
%
% Merge ESRL GML CH4 and MIPAS CH4 
%       MIPAS CH4 avalable for period 2002-2014 only
% Retain original pressure/height before passing to klayers
%       unlike original fill_ch4.m matching to ECMWF model fields
%       which sacrifices upper levels.

addpath /home/chepplew/myLib/matlib/readers       % read_netcdf

% ===================================================
% home of the MIPAS CH4 - get year/month of filenames
% ===================================================
mip = struct;
mip.home = '/home/chepplew/data/MIPAS/CH4/';
mip.list = dir([mip.home 'MIPAS-E_IMK.20*_CH4_61.nc']);
for ifn=1:length(mip.list)
  fname = mip.list(ifn).name;
  junk = strsplit(fname,{'-','_','.'});
  mip.fndnum(ifn) = datenum(junk{4},'yyyymm');
end

% ==================================================
% The ESRL Reference CH4 distribution at 2010.12.31
% ==================================================
ch4 = struct;
ch4.ref_fn = '/home/chepplew/data/GML_CH4/ch4_molefractions/20101231.nc';

[sac aac] = read_netcdf(ch4.ref_fn);

% Total CH4 ppbv
ch4.tot = sac.bgrnd + sac.fossil + sac.agwaste + sac.natural + sac.bioburn + sac.ocean;

% Get global mean (near) surface ppbv at this reference point:
ch4.gbl_mn0 = mean(ch4.tot(:,:,1:3,:),[1 2 3 4]);

% Average the 8 time & longitude samples and get the number of levels:
ch4.tot_tmn = squeeze(nanmean(ch4.tot,4));
ch4.plevs   = squeeze(nanmean(sac.pressure,4));
ch4.nlevs   = size(sac.pressure,3);

% ============================================
% Load CH4 global annual mean trend from file:
% ============================================
ch4.gblmn_fn = '/home/chepplew/data/GML_CH4/trends/ch4_annmean_gl.txt';
D = importdata(ch4.gblmn_fn,' ',62);
ch4.annmean  = D.data(:,2);
ch4.annyrs   = datenum(num2str(D.data(:,1)),'yyyy');   % ! January or June?
clear D;

% =========================================
% Check requested date and get derivatives
% =========================================
% %rdate = '2003/11/21';
try
  dnum = datenum(rdate,'yyyy/mm/dd');
catch ME
  error(ME.identifier)
  return
end
dtime = datetime(dnum,'convertfrom','datenum');
uyr   = year(dtime);
umn   = month(dtime);
udy   = day(dtime);
ujdy  = day(dtime,'dayofyear');

% =========================================
% scale GML reference profile as needed
% =========================================
xdim = ch4.annyrs - ch4.annyrs(1);
stim = dnum - ch4.annyrs(1);
junk = interp1(xdim, ch4.annmean, stim, 'linear','extrap');

% =========================================
% scale MIPAS profiles if outside of period
% =========================================
if(dnum < mip.fndnum(1))
   

elseif( dnum > mip.fndnum(end))
  stim = mip.fndnum(end) - ch4.annyrs(1);
  junk = interp1(xdim, ch4.annmean, stim, 'linear','extrap');


% ===========================================
%  match the MIPAS and ESRL profiles
% ===========================================



