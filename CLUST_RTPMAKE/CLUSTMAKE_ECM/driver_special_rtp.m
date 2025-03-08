addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta/

%% for rtp ECM or ERA, see cris_l1c_to_rtp_sergio.m

addpath ~motteler/shome/cris/ccast/source/          % fixmyQC
addpath /asl/matlib/rtptools                        % set_attr
addpath /asl/matlib/time
addpath /asl/matlib/aslutil                         % int2bits
addpath /home/sbuczko1/git/rtp_prod2/chirp/util/uniform_clear/
addpath /home/sbuczko1/git/rtp_prod2/grib/          % fill_era
addpath /home/sbuczko1/git/rtp_prod2/emis
addpath /home/sbuczko1/git/rtp_prod2/util           % genscratchpath
addpath /home/sbuczko1/git/rtp_prod2/cris/util      % guard_ind
addpath /home/sbuczko1/git/rtp_prod2/cris/util/uniform_clear
addpath /home/chepplew/projects/chirp               % cat_rtp_clh
addpath /home/sergio/MATLABCODE

system_slurm_stats

if ~exist('iInterp')
  iInterp = -1;
  iInterp = +1;
end

if ~exist('iERAorECM')
  iERAorECM = +1; %% till June 2019
  iERAorECM = -1; %% after June 2019
end

if ~exist('iSNPPorJ1orJ2')
  iSNPPorJ1orJ2 = +0; %% SuomiNPP
  iSNPPorJ1orJ2 = +1; %% J1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fin  = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr.ip.rtp';
fout = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr.ip_new.rtp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

[h,ha,p,pa] = rtpread(fin);
[yy,mm,dd,hh] = tai2utcSergio(p.rtime);

h0 = h;
p0 = p;
pattr = pa;
  
if h.ptype == 1
  error('only handling input levels file')
end
if h.ptype == 0
  for ii = 1 : h.ngas
    junkstr = ['gas_' num2str(h.glist(ii))];
    p = rmfield(p,junkstr);
  end
  p = rmfield(p,'plevs');
  %p = rmfield(p,'palts');
  p = rmfield(p,'ptemp');
  p = rmfield(p,'stemp');
  p = rmfield(p,'wspeed');
  p = rmfield(p,'spres');

  h.ngas = 0;  
  h = rmfield(h,'glist');
  h = rmfield(h,'gunit');
  h.pfields = 4; % 4 = IR obs

end
hIN = h;
pIN = p;

ee = dir([fout]);
if length(ee) == 0

  fn_rtp1 = mktempS('fx.ip.rtp');
  fn_rtp2 = mktempS('fx.op.rtp');
  fn_rtp3 = mktempS('fx.rp.rtp');
  fn_rtp4 = mktempS('fx.xp.rtp');

  %%%%%%%%%%%%%%%%%%%%%%%%%
  cfg.model = 'ecmwf';  
  switch cfg.model
    case 'ecmwf'
      if iInterp <= 0
        [p,h,pattr]  = fill_ecmwf(p,h,pattr);
      else
        new_8_ECMfiles_interp_analysis
      end

    case 'era'
      if iInterp <= 0
        [p,h,pattr]  = fill_era(p,h,pattr);
        %[p,h,pattr]  = fill_era_interp(p,h,pattr);
      else
        new_4_ERAfiles_interp_analysis
      end

    case 'merra'
      [p,h,pattr]  = fill_merra(p,h,pattr);

  end
  %%%%%%%%%%%%%%%%%%%%%%%%%

%   i900 = find(h.vchan >= 900,1);
%   tobs = rad2bt(900,p.robs1(i900,:));
%   tclr = rad2bt(900,p.sarta_rclearcalc(i900,:));
%   tcld = rad2bt(900,p.rcalc(i900,:));
%   addpath /home/sergio/MATLABCODE/PLOTTER
%   figure(1); clf; scatter_coast(p.rlon,p.rlat,25,tobs); title('BT 900 obs FSR'); cx1 = caxis; colormap jet
%   figure(2); clf; scatter_coast(p.rlon,p.rlat,25,tclr); title('BT 900 clr');     cx2 = caxis; colormap jet
%   figure(3); clf; scatter_coast(p.rlon,p.rlat,25,tcld); title('BT 900 cld');     cx3 = caxis; colormap jet
% 
%   cx(1) = min([cx1(1)  cx2(1) cx3(1)]);
%   cx(2) = max([cx1(2)  cx2(2) cx3(2)]);
%   figure(1); caxis(cx); 
%   figure(2); caxis(cx); 
%   figure(3); caxis(cx); 

%   %rtpwrite(fout,h,ha,p,pa)
  semilogy(nanmean(p0.ptemp,2),nanmean(p0.plevs,2),'b.-',nanmean(p.ptemp,2),nanmean(p.plevs,2),'r'); set(gca,'ydir','reverse')
  semilogy(nanmean(p0.ptemp,2)-nanmean(p.ptemp,2),nanmean(p.plevs,2),'r'); set(gca,'ydir','reverse')
  loglog(nanmean(p0.gas_1,2),nanmean(p0.plevs,2),'b.-',nanmean(p.gas_1,2),nanmean(p.plevs,2),'r'); set(gca,'ydir','reverse')
  fprintf(1,'DONE : %s written out  \n',[fout])  
else
  fprintf(1,'%s already exists \n',[fout])
end
