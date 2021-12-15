addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/PLOTTER

figure(1); clf;
figure(2); clf;

ii0 = input('Enter gran to check : ');

while ii0 > 0 & ii0 < 241
  f0 = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/10/31/cloudy_airs_l1b_ecm_sarta_baum_ice.2018.10.31.' num2str(ii0,'%03d') '.rtp'];
  fx = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/10/31/cloudy_airs_l1b_ecm_sarta_baum_ice.2018.10.31.' num2str(ii0,'%03d') '_pert_tcc_model_1.rtp'];

  i0 = -1;
  if exist(f0)
    thedir = dir(f0);
    if thedir.bytes > 4e8
      i0 = +1;
    end
  end

  ix = -1;
  if exist(fx)
    thedir = dir(fx);
    if thedir.bytes > 4e8
      ix = +1;
    end
  end

  if i0 > 0 & ix > 0
    [h,ha,p0,pa] = rtpread(f0);
    [h,ha,px,pa] = rtpread(fx);
    tobs = rad2bt(h.vchan,p0.robs1);
    tcal0 = rad2bt(h.vchan,p0.rcalc);
    tcalx = rad2bt(h.vchan,px.rcalc);

    figure(1);
    plot(h.vchan,nanmean(tobs'-tcal0'),'b',h.vchan,nanmean(tobs'-tcalx'),'r',h.vchan,nanmean(tcal0'-tcalx'),'k',h.vchan,nanstd(tobs'-tcal0'),'c',h.vchan,nanstd(tobs'-tcalx'),'m',h.vchan,nanstd(tcal0'-tcalx'),'g')
    hl = legend('tobs-tcal0','tobs-tcalx','tcal0-tcalx','location','best'); set(hl,'fontsize',10); grid

    figure(2); scatter_coast(p0.rlon,p0.rlat,10,p0.stemp);

    ii0 = input('Enter gran to check (1:240) : ');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue'); pause
tobs1231 = [];
tclr1231 = [];
tcal1231_0 = [];
tcal1231_x = [];
stemp0 = [];
stempx = [];
tcc0 = [];
tccx = [];
rlat = [];
rlon = [];
gran = [];

for ii = 1 : 240
  fprintf(1,'.');

  f0 = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/10/31/cloudy_airs_l1b_ecm_sarta_baum_ice.2018.10.31.' num2str(ii,'%03d') '.rtp'];
  fx = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/10/31/cloudy_airs_l1b_ecm_sarta_baum_ice.2018.10.31.' num2str(ii,'%03d') '_pert_tcc_model_1.rtp'];

  i0 = -1;
  if exist(f0)
    thedir = dir(f0);
    if thedir.bytes > 4e8
      i0 = +1;
    end
  end

  ix = -1;
  if exist(fx)
    thedir = dir(fx);
    if thedir.bytes > 4e8
      ix = +1;
    end
  end

  if i0 > 0 & ix > 0
    [h,ha,p0,pa] = rtpread(f0);
    [h,ha,px,pa] = rtpread(fx);
    tobs = rad2bt(h.vchan,p0.robs1);
    tclr = rad2bt(h.vchan,p0.sarta_rclearcalc);
    tcal0 = rad2bt(h.vchan,p0.rcalc);
    tcalx = rad2bt(h.vchan,px.rcalc);

    tobs1231   = [tobs1231   tobs(1520,:)];
    tclr1231   = [tclr1231   tclr(1520,:)];
    tcal1231_0 = [tcal1231_0 tcal0(1520,:)];
    tcal1231_x = [tcal1231_x tcalx(1520,:)];
    stemp0 = [stemp0 p0.stemp];
    stempx = [stempx px.stemp];
    tcc0 = [tcc0 p0.tcc];
    tccx = [tccx px.tcc];
    gran = [gran ones(size(p0.stemp))*ii];
    rlon = [rlon p0.rlon];    
    rlat = [rlat p0.rlat];    

    bias0(ii)  = nanmean(tobs1231-tcal1231_0);
    biasx(ii)  = nanmean(tobs1231-tcal1231_x);
    bias0x(ii) = nanmean(tcal1231_0-tcal1231_x);
    std0(ii)  = nanstd(tobs1231-tcal1231_0);
    stdx(ii)  = nanstd(tobs1231-tcal1231_x);
    std0x(ii) = nanstd(tcal1231_0-tcal1231_x);
    meanobs1231(ii) = nanmean(tobs1231);
    meancld1231(ii) = nanmean(tobs1231-tclr1231);
    meanrlat(ii) = mean(p0.rlat);

%    whos tobs1231 tcal1231_* rlon rlat stemp* tcc*

    if mod(ii,10) == 0
      figure(2); plot(meanrlat,bias0,'b.',meanrlat,biasx,'r.',meanrlat,bias0x,'g.')
      plotaxis2;
      pause(0.1)
    end
  end
end
fprintf(1,'\n');
save compare_NWPtcc_vs_homebrewtcc_rads_cumsum_9999.mat t*1231* stemp* tcc* gran rlon rlat bias* std* mean*
a = load('compare_NWPtcc_vs_homebrewtcc_rads_cumsum_9999.mat')

dbt = -50 : 0.1 : 50;
nx = hist(tobs1231-tcal1231_x,dbt);
n0 = hist(tobs1231-tcal1231_0,dbt);
n0x = hist(tcal1231_0-tcal1231_x,dbt);
figure(2); plot(dbt,n0,dbt,nx,dbt,n0x); grid
figure(2); semilogy(dbt,n0,'b',dbt,nx,'r',dbt,n0x,'g','linewidth',2); grid

figure(1); 
  subplot(211); plot(meanrlat,bias0,'bx',meanrlat,biasx,'ro',meanrlat,bias0x,'gd'); 
  plotaxis2; ylabel('bias (K)');
  title('1231 cm-1 average granule stats')
  subplot(212); plot(meanrlat,std0,'bx',meanrlat,stdx,'ro',meanrlat,std0x,'dg'); plotaxis2; ylabel('std (K)')
  xlabel('Latitude');
