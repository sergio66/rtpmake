N = 29; N1 = 1;
firstORend = 0;
iPrint = -1;

yy0 = 2002; mm0 = 09; dd0 = 15;
thedateS(1,:) = [yy0 mm0 dd0];

[yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,N,firstORend,iPrint);
dd1 = 15;
thedateE(1,:) = [yy1 mm1 dd1];

for ii = 2 : JOB
  yy0 = yy1; mm0 = mm1; dd0 = dd1;
  thedateS(ii,:) = [yy0 mm0 dd0];  

  [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,N,firstORend,iPrint);
  dd1 = 15;
  thedateE(ii,:) = [yy1 mm1 dd1];

  fprintf(1,'ii=%4i   Start %4i/%2i/%2i  End %4i/%2i/%2i \n',ii,thedateS(ii,:),thedateE(ii,:))
end  

[thedateS thedateE]
[yyM,mmM,ddM] = addNdays(thedateS(JOB,1),thedateS(JOB,2),thedateS(JOB,3),8,firstORend,iPrint);

fprintf(1,'JOB = %3i spans %4i/%2i/%2i to %4i/%2i/%2i both ends inclusive \n',JOB,[thedateS(JOB,:) thedateE(JOB,:)]);
fprintf(1,'          midpoint %4i/%2i/%2i \n',yyM,mmM,ddM);

rtime0 = utc2taiSergio(thedateS(JOB,1),thedateS(JOB,2),thedateS(JOB,3),0.00);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now we need to get overpass times and solzen angles, just set scanang to 22 deg
load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72.mat');
load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_457_64x72.mat');
monitor_memory_whos;

thedata.avgrtime_desc = nanmean(squeeze(nanmean(thedata.rtime_desc,1)),1);
thedata.avgrtime_asc  = nanmean(squeeze(nanmean(thedata.rtime_asc,1)),1);
if iDorA > 0
  JOB_23timesteps = find(thedata.avgrtime_desc >= rtime0,1);
else
  JOB_23timesteps = find(thedata.avgrtime_asc >= rtime0,1);
end
fprintf(1,'JOB = %3i in terms of months corresponds to 16 day timestep number %3i \n',JOB,JOB_23timesteps);
fprintf(1,'really should set JOB --> JOB_23timesteps below!!! \n')
JOBx = JOB_23timesteps; %% this should really be uncommented and     used; after  Feb 2024
JOBx = JOB;             %% this should reallt be   commented and not used; before Feb 2024 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% from comment, see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/driver_loop_get_asc_desc_solzen_time.m
if iDorA > 0
  rlon   = thedata.rlon_desc(:,:,JOBx);     rlon = rlon(:)';
  rlat   = thedata.rlat_desc(:,:,JOBx);     rlat = rlat(:)';
  solzen = thedata.solzen_desc(:,:,JOBx);   solzen = solzen(:)';
  satzen = thedata.satzen_desc(:,:,JOBx);   satzen = satzen(:)';
  hour   = thedata.hour_desc(:,:,JOBx);     hour = hour(:)';
  bt1231 = thedata.bt1231_desc(:,:,JOBx,2); bt1231 = bt1231(:)';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  rtime  = thedata.rtime_desc(:,:,JOBx);    rtime = rtime(:)';
    [mooY,mooM,mooD,mooH] = tai2utcSergio(rtime);
    junk = [JOB JOBx iDorA round(mean(mooY)) round(mean(mooM)) round(mean(mooD)) mean(mooH) yyM mmM ddM];
    fprintf(1,'JOB = %3i ... using JOBx = %3i in iDorA = %2i we get mean thedata.rtime_Xsc yy/mm/dd = %4i/%2i/%2i at HH = %8.4f instead of %4i/%2i/%2i \n',junk);
 
    rtime0 = utc2taiSergio(mooY(1),mooM(1),mooD(1),0.00);  drtime =  (rtime-rtime0);
  
    rtimex = utc2taiSergio(thedateS(JOB,1),thedateS(JOB,2),thedateS(JOB,3),0.00) + drtime;
    [xmooY,xmooM,xmooD,xmooH] = tai2utcSergio(rtimex);
    junk = [round(mean(xmooY)) round(mean(xmooM)) round(mean(xmooD)) mean(xmooH) yyM mmM ddM];
    fprintf(1,'  then we reset to get new rtime yy/mm/dd = %4i/%2i/%2i at HH = %8.4f hopefully close to %4i/%2i/%2i \n',junk);  
  rtime = rtimex;
  %%%%%%%%%%%%%%%%%%%%%%%%%

elseif iDorA < 0
  rlon   = thedata.rlon_asc(:,:,JOBx);     rlon = rlon(:)';
  rlat   = thedata.rlat_asc(:,:,JOBx);     rlat = rlat(:)';
  solzen = thedata.solzen_asc(:,:,JOBx);   solzen = solzen(:)';
  satzen = thedata.satzen_asc(:,:,JOBx);   satzen = satzen(:)';
  hour   = thedata.hour_asc(:,:,JOBx);     hour = hour(:)';
  bt1231 = thedata.bt1231_asc(:,:,JOBx,2); bt1231 = bt1231(:)';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  rtime  = thedata.rtime_asc(:,:,JOBx);    rtime = rtime(:)';
    [mooY,mooM,mooD,mooH] = tai2utcSergio(rtime);
    junk = [JOB JOBx iDorA round(mean(mooY)) round(mean(mooM)) round(mean(mooD)) mean(mooH) yyM mmM ddM];
    fprintf(1,'JOB = %3i ... using JOBx = %3i in iDorA = %2i we get mean thedata.rtime_Xsc yy/mm/dd = %4i/%2i/%2i at HH = %8.4f instead of %4i/%2i/%2i \n',junk);

    rtime0 = utc2taiSergio(mooY(1),mooM(1),mooD(1),0.00);  drtime =  (rtime-rtime0);
  
    rtimex = utc2taiSergio(thedateS(JOB,1),thedateS(JOB,2),thedateS(JOB,3),0.00) + drtime;
    [xmooY,xmooM,xmooD,xmooH] = tai2utcSergio(rtimex);
    junk = [round(mean(xmooY)) round(mean(xmooM)) round(mean(xmooD)) mean(xmooH) yyM mmM ddM];
    fprintf(1,'  then we reset to get new rtime yy/mm/dd = %4i/%2i/%2i at HH = %8.4f hopefully close to %4i/%2i/%2i \n',junk);  
  
  rtime = rtimex;
  %%%%%%%%%%%%%%%%%%%%%%%%%
end

figure(1); scatter_coast(rlon,rlat,50,solzen); colormap jet
figure(2); scatter_coast(rlon,rlat,50,satzen); colormap jet
figure(3); scatter_coast(rlon,rlat,50,hour); colormap jet
figure(4); scatter_coast(rlon,rlat,50,bt1231); colormap jet
figure(5); scatter_coast(rlon,rlat,50,rlon-p.rlon); colormap jet
figure(6); scatter_coast(rlon,rlat,50,rlat-p.rlat); colormap jet
