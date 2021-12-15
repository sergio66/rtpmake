function [p,FOV_ECM,fxx,alltOBS,alltECM,closestBT] = loop_get_bestL2orECM_FAST_5deglat(head0,prof0,g,iSwapAllFields,iStemp_ColWV,tempx,iBadDCCOnly);  %% swaps all cloud fields

%% iSwapAllFields = 0 ==> swap all fields
%% iSwapAllFields = 1 ==> swap all fields EXCEPT for cprtop for water

%% accounts for ST, WV and LF

%{
%% testing as a script
head0 = h;
prof0 = p;
g = g;
iSwapAllFields = -1;
iStemp_ColWV = iStemp_ColWV;
tempx = tempx;     %% bunch of temporary filenames
%}


addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
RH0 = layeramt2RH(head0,prof0);
prof1 = prof0; prof1.ptemp = prof1.ptemp + 1;
RH1 = layeramt2RH(head0,prof1);

mmw = mmwater_rtp(head0,prof0);

ind_trop = tropopause_rtp(head0,prof0);
for ii = 1 : length(prof0.stemp)
  junk = prof0.ptemp(:,ii);
  tropopauseT(ii) = junk(ind_trop(ii));
  junk = prof0.plevs(:,ii);
  tropopauseP(ii) = junk(ind_trop(ii));
end
  
h = head0;
p = prof0;

pcalc = p.rcalc; [mmjunk,nnjunk] = size(pcalc); fprintf(1,'here 0 size(p.rcalc) = %5i x %5i \n',mmjunk,nnjunk)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch_ind_xyz = get_retrieval_chans(h,g,iStemp_ColWV);

ch_ind_use = g;
ch_ind_use = ch_ind_xyz;

%ch_ind_use = ch_ind_use(ch_ind_use < 1620);        %% keep to LW and MW
zooLW = find(ch_ind_use >=  246 & ch_ind_use <= 900);   %% keep to 720  cm-1 < ch < 960 cm-1
zooMW = find(ch_ind_use >= 1431 & ch_ind_use <= 1864);  %% keep to 1320 cm-1 < ch < 1620 cm-1
oo = union(zooLW,zooMW(1:3:end));

%% use only window chans for this cloud stuff
oo = find((h.vchan(ch_ind_xyz) > 743 & h.vchan(ch_ind_xyz) < 970) | (h.vchan(ch_ind_xyz) > 1180 & h.vchan(ch_ind_xyz) < 1280));

addpath /home/sergio/MATLABCODE/DUSTFLAG
freqsNchans;
ch_ind_use = intersect(g,cind1(1:5));
ch_ind_use = cind1(1:5);;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'  for mapping cloudfields, will be looking at %4i chans \n',length(ch_ind_use));
plot(head0.vchan(ch_ind_use),nanmean(rad2bt(head0.vchan(ch_ind_use),prof0.robs1(ch_ind_use,:))'),'o-'); title('mean over chans chosen for cloud mapping');
pause(0.1);

alltOBS = real(rad2bt(h.vchan(ch_ind_use),p.robs1(ch_ind_use,:)));
alltECM = real(rad2bt(h.vchan(ch_ind_use),p.rcalc(ch_ind_use,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%
[hxx,pxx] = subset_rtp_allcloudfields(h,p,[],ch_ind_use,[]);  %% go from 2378 chans to about 5 chans
%hxx.ichan
%hxx.vchan
%%%%%%%%%%%%%%%%%%%%%%%%%

i899 = find(hxx.vchan >= 899,1);
tobs899 = rad2bt(900,pxx.robs1(i899,:));
cldforce = pxx.stemp - tobs899;

tcal899 = rad2bt(900,pxx.rcalc(i899,:));
figure(1); clf; scatter(cldforce,max(pxx.cfrac,pxx.cfrac2),10,pxx.rlat,'filled'); colorbar
  xlabel('cldforce = stemp-BTobs900'); ylabel('cfrac'); title('colorbar = rlat')
xx = -10:100; plot(xx,max(0,tanh((xx-60)/5))); grid
cfracXFORC = max(0,tanh((cldforce-55)/5));  %% force cldfrac = 1 if cldforc > 55 K
check_these_force = find(cldforce > 55);
check_these_force = find(cfracXFORC > 0.9);
check_these_force = find(cfracXFORC > 0.5);

[tcld1,tcld2] = do_cld1_cld2(pxx,tobs899,tcal899,iBadDCCOnly);
%disp('ret 0'); pause

iNotOK = +1;
iFix = 0;
%% before used to give up at iFix == 10
while iFix < 12 & iNotOK > 0
  iFix = iFix + 1;
  fprintf(1,' doing n = %2i try at checking clouds \n',iFix)
  [pxx,iNotOK] = check_for_errors(pxx,-1,iFix);  %% see possible pitfalls in clouds
end
if iFix >= 12 & iNotOK > 0
  %disp('oops, could not fix cprtop vs cprbot vs spres'); %keyboard
  error('oops, could not fix cprtop vs cprbot vs spres')
end	  
if hxx.ptype == 0
  rtpwrite(tempx.tempfile1,hxx,[],pxx,[]);
  klayerser = ['!' tempx.klayers ' fin=' tempx.tempfile1 ' fout=' tempx.tempfile2 ' >& ' tempx.tempfile4]; eval(klayerser)
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
else
  rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
end
eval(['!tail -100 ' tempx.tempfile4])
[hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);

iDoSarta = -1;
newalltECM = real(rad2bt(hxx.vchan,pxx.rcalc));
if abs(sum(sum(newalltECM-alltECM))) > 1e-5
  plot(hxx.vchan,nanmean(newalltECM'-alltECM'),hxx.vchan,nanstd(newalltECM'-alltECM'))
  disp('warning : the SARTA cloud calcs in prof0 are (slightly) different than newest/right now calcs ...')
  [mmm,nnn] = size(alltECM);
  fprintf(1,'mean diff over %4i chans, %5i profs = %8.6f K \n',mmm,nnn,abs(sum(sum(newalltECM-alltECM))) / (mmm * nnn));
  alltECM = newalltECM;
end

% also see /home/sergio/MATLABCODE/matrix_nearest_neighbour.m
% also see /home/sergio/MATLABCODE/matrix_nearest_neighbour.m
% also see /home/sergio/MATLABCODE/matrix_nearest_neighbour.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% closest point in B for each point in A
%% uses up memory
% distances = squeeze(sum(bsxfun(@minus,alltOBS(:,1:1000)',permute(alltECM(:,1:1000)',[3 2 1])).^2,2));
% distances = squeeze(sum(bsxfun(@minus,alltOBS',permute(alltECM',[3 2 1])).^2,2));
% [~,ik] = min(distances,[],2)
% Make an array the size of A containing these closest points in B:
% closestBT = alltECM(ik,:)

% OR

lff = ones(size(prof0.landfrac));
oo = find(prof0.landfrac == 0); lff(oo) = -1;
oo = find(prof0.landfrac >  0); lff(oo) = +1;   
xalltOBS = [alltOBS; prof0.stemp; 250+mmw; 250+lff*50];

%% loop this iN times
iN = 5;
iN = 2;
iN = 3;
for loop = 1 : iN
  if iDoSarta > 0
    %% avoid running SARTA the first time, since you've already done it above
    pxx = p;
    pxx = sanity_check(pxx);    
    rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
    [hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);
    alltECM = real(rad2bt(hxx.vchan,pxx.rcalc));

    tcal899 = rad2bt(900,pxx.rcalc(i899,:));
    figure(1); clf; scatter(cldforce,max(pxx.cfrac,pxx.cfrac2),10,pxx.rlat,'filled'); colorbar
    xlabel('cldforce = stemp-BTobs900'); ylabel('cfrac'); title('colorbar = rlat')
    xx = -10:100; plot(xx,max(0,tanh((xx-60)/5))); grid
    [tcld1,tcld2] = do_cld1_cld2(pxx,tobs899,tcal899,iBadDCCOnly);
    %disp('ret N'); pause

    %% now mark which of the cldforc > 55 have tcld1 > tobs900
    tBTD_DCC_bad = 2;
    oops0 = find(cfracXFORC > 0.5);
    oops1 = find(cfracXFORC > 0.5 & ((tcld1 - tobs899) > tBTD_DCC_bad));  %% bad dcc
    oops2 = find(cfracXFORC < 0.5 & ((tcld1 - tobs899) > tBTD_DCC_bad));  %% bad mid-lower trop clouds
    if length(oops1) > 0 | length(oops2) > 0    
      rcalc = p.rcalc;
      robs1 = p.robs1;
      calflag = p.calflag;      
      fprintf(1,' found %5i of %5i potential "incorrect cloudset" DCC   ie where cfracXFORC > 0.5 & (tcld1 - tobs899 > 2) \n',length(oops1),length(oops0));
      fprintf(1,' found        %5i potential "incorrect cloudset" OTHER ie where cfracXFORC < 0.5 & (tcld1 - tobs899 > 2) \n',length(oops2));
      %save junk.mat hxx pxx cfracXFORC tcld1 tobs899 tBTD_DCC_bad
      pxx = move_cld1_cld2(hxx,pxx,oops1,oops2,tobs899,tempx,tBTD_DCC_bad,RH0,RH1);
      alltECM = real(rad2bt(hxx.vchan,pxx.rcalc));
      p = pxx;   %% HAVE TO DO THIS since we changed some of the cprtop and ptemp
      p.robs1 = robs1;
      p.rcalc = rcalc;
      p.calflag = calflag;      
      
      tcal899 = rad2bt(900,pxx.rcalc(i899,:));
      figure(1); clf; scatter(cldforce,max(pxx.cfrac,pxx.cfrac2),10,pxx.rlat,'filled'); colorbar
      xlabel('cldforce = stemp-BTobs900'); ylabel('cfrac'); title('colorbar = rlat')
      xx = -10:100; plot(xx,max(0,tanh((xx-60)/5))); grid
      [tcld1,tcld2] = do_cld1_cld2(pxx,tobs899,tcal899,iBadDCCOnly);
      %disp('ret N (move cld)'); pause
    end
    
  end
  iDoSarta = +1;
  xalltECM = [alltECM; prof0.stemp; 250+mmw; 250+lff*50];

  %[damin0, ind0] = min(pdist2(alltECM',alltOBS','euclidean'));  toc2a = toc;    %%% ORIG
  [damin, ind] = min(pdist2(xalltECM',xalltOBS','euclidean')); toc2b = toc;    %%% NEW
  fprintf(1,'   loop %2i of 5 --> min = %8.6f \n',loop,sum(damin))
  
  %latspread = max(p.rlat)-min(p.rlat);
  %if latspread > 15  %% llooks like completely ascending or descending granule
  %  disp('A/D granule ...)  
  %else
  %  disp('polar granule ...)
  %  [~, ind] = min(pdist2(xalltECM',xalltOBS','euclidean')); toc2c = toc;    %%% NEW  
  %end
  %error('kjslkjgslkjgs 22')

  FOV_ECM = ind;
  figure(1); clf; plot(FOV_ECM); title('mapping indices to BEST cloud calcs')
  figure(2); clf; plot(1:length(p.stemp),tropopauseT,'b',1:length(p.stemp),tropopauseT(ind),'r')
  figure(3); clf; plot(1:length(p.stemp),tropopauseP,'b',1:length(p.stemp),tropopauseP(ind),'r')

  closestBT = alltECM(:,ind);
  fxx       = h.vchan(ch_ind_use);

  figure(2); clf
  plot(h.vchan(ch_ind_use),nanmean(alltOBS'-alltECM'),'b',h.vchan(ch_ind_use),nanstd(alltOBS'-alltECM'),'c--',...
     h.vchan(ch_ind_use),nanmean(alltOBS'-closestBT'),'r',h.vchan(ch_ind_use),nanstd(alltOBS'-closestBT'),'m--')
  xlabel('wavenumber cm-1'); ylabel('BT(K)'); title('from ECM (b/c) and remapped (r/m)'); grid

  disp('  remapping cloud fields .....')
  p.ctype = p.ctype(FOV_ECM);
  p.cfrac = p.cfrac(FOV_ECM);
  p.cprtop = p.cprtop(FOV_ECM);
  p.cprbot = p.cprbot(FOV_ECM);
  p.cngwat = p.cngwat(FOV_ECM);
  p.cpsize = p.cpsize(FOV_ECM);

  p.ctype2 = p.ctype2(FOV_ECM);
  p.cfrac2 = p.cfrac2(FOV_ECM);
  p.cprtop2 = p.cprtop2(FOV_ECM);
  p.cprbot2 = p.cprbot2(FOV_ECM);
  p.cngwat2 = p.cngwat2(FOV_ECM);
  p.cpsize2 = p.cpsize2(FOV_ECM);

  p.cfrac12 = p.cfrac12(FOV_ECM);

  figure(1); clf; scatter(cldforce,max(p.cfrac,p.cfrac2),10,p.rlat,'filled'); colorbar
    xlabel('cldforce = stemp-BTobs900'); ylabel('cfrac'); title('BEFORE map (colorbar = rlat)')
  %% find the highest cloud
  boo = find(p.cprtop < p.cprtop2 & p.cprtop > 0);
    p.cfrac(boo) = max(cfracXFORC(boo),p.cfrac(boo));
  boo = find(p.cprtop2 < p.cprtop & p.cprtop2 > 0);
    p.cfrac2(boo) = max(cfracXFORC(boo),p.cfrac2(boo));
  figure(2); clf; scatter(cldforce,max(p.cfrac,p.cfrac2),10,p.rlat,'filled'); colorbar
    xlabel('cldforce = stemp-BTobs900'); ylabel('cfrac'); title('AFTER (colorbar = rlat)')
  %disp('ret'); pause
  pause(0.1)
  
  if isfield(p,'clwc'); p.clwc = p.clwc(:,FOV_ECM); end
  if isfield(p,'ciwc'); p.ciwc = p.ciwc(:,FOV_ECM); end
  if isfield(p,'cc');   p.cc   = p.cc(:,FOV_ECM);   end

  if isfield(p,'orig_ctop');  p.orig_ctop  = p.orig_ctop2(FOV_ECM); end
  if isfield(p,'orig_ctop2'); p.orig_ctop2 = p.orig_ctop2(FOV_ECM); end
  if isfield(p,'sarta_lvlODice');   p.sarta_lvlODice =   p.sarta_lvlODice(:,FOV_ECM); end
  if isfield(p,'sarta_lvlODwater'); p.sarta_lvlODwater = p.sarta_lvlODwater(:,FOV_ECM); end

  if isfield(p,'cemis');   p.cemis  = p.cemis(:,FOV_ECM); end
  if isfield(p,'crho');    p.crho   = p.crho(:,FOV_ECM); end
  if isfield(p,'cstemp');  p.cstemp   = p.cstemp(:,FOV_ECM); end
  if isfield(p,'cemis2');  p.cemis2  = p.cemis2(:,FOV_ECM); end
  if isfield(p,'crho2');   p.crho2   = p.crho2(:,FOV_ECM); end
  if isfield(p,'cstemp3'); p.cstemp2 = p.cstemp2(:,FOV_ECM); end

  p.rcalcBestCld = p.rcalc(:,FOV_ECM); %%% <<<< this is assumed best calc after swap, who knows what will happen after
                                     %%% <<<< re-running sarta-cloudy

  pcalc=p.rcalc; [mmjunk,nnjunk] = size(pcalc); fprintf(1,'here loop %2i size(p.rcalc) = %5i x %5i \n',loop,mmjunk,nnjunk)
  pxx = p;
end

pnew0 = p;

if iSwapAllFields == +1
  %% keep cloud top for water the same
  disp('You Asked To Keep cloud top for water the same')
  Bonk = find(p.ctype == 101 & p.cngwat > 0 & prof0.ctype == 101 & prof0.cngwat > 0)
    p.cprtop(bonk) = prof0.cprtop(bonk);
    p.cprbot(bonk) = prof0.cprbot(bonk);
    p.cngwat(bonk) = prof0.cngwat(bonk);
  bonk = find(p.ctype2 == 101 & p.cngwat2 > 0 & prof0.ctype2 == 101 & prof0.cngwat2 > 0)    
    p.cprtop2(bonk) = prof0.cprtop2(bonk);
    p.cprbot2(bonk) = prof0.cprbot2(bonk);
    p.cngwat2(bonk) = prof0.cngwat2(bonk);
end

%p.rcalcBestCld = p.rcalc(:,FOV_ECM); %%% <<<< this is assumed best calc after swap, who knows what will happen after
%                                     %%% <<<< re-running sarta-cloudy

pnew1 = p;

iSwapPlot = -1;
if iSwapPlot > 0

  donk = find(p.ctype == 101 & p.cngwat > 0); 
  donk0 = find(prof0.ctype == 101 & prof0.cngwat > 0); 
  donknew = find(pnew0.ctype == 101 & pnew0.cngwat > 0); 
  figure(1); clf; pp = 0:25:1000; 
  plot(pp,hist(p.cprtop(donk),pp)/sum(hist(p.cprtop(donk),pp)),'r', ...
       pp,hist(prof0.cprtop(donk0),pp)/sum(hist(prof0.cprtop(donk0),pp)),'b',...
       pp,hist(pnew0.cprtop(donknew),pp)/sum(hist(pnew0.cprtop(donknew),pp)),'k')
  donk = find(p.ctype2 == 101 & p.cngwat2 > 0); 
  donk0 = find(prof0.ctype2 == 101 & prof0.cngwat2 > 0); 
  donknew = find(pnew0.ctype2 == 101 & pnew0.cngwat2 > 0); 
  figure(2); clf; pp = 0:25:1000; 
  plot(pp,hist(p.cprtop2(donk),pp)/sum(hist(p.cprtop2(donk),pp)),'r',...
       pp,hist(prof0.cprtop2(donk0),pp)/sum(hist(prof0.cprtop2(donk0),pp)),'b',...
       pp,hist(pnew0.cprtop2(donknew),pp)/sum(hist(pnew0.cprtop2(donknew),pp)),'k')

% debug_smoother
% keyboard

  figure(1); clf; scatter(p.rlon,p.rlat,30,rad2bt(1231,p.robs1(1291,:))); title('BT1231 obs')
  figure(2); clf; scatter(p.rlon,p.rlat,30,p.cprtop); title('new cprtop'); caxis([0 1000]); colorbar
  figure(3); clf; scatter(p.rlon,p.rlat,30,p.cprtop2); title('new cprtop2'); caxis([0 1000]); colorbar
  figure(4); clf; scatter(p.rlon,p.rlat,30,log10(p.cngwat)); title('new log10cngwat)')
  figure(5); clf; scatter(p.rlon,p.rlat,30,log10(p.cngwat2)); title('new log10(cngwat2)')
  figure(6); clf; scatter(p.rlon,p.rlat,30,p.cfrac); title('new cfrac')
  figure(7); clf; scatter(p.rlon,p.rlat,30,p.cfrac2); title('new cfrac2')

end

pBefore = p;
figure(1); colormap jet; bwah1 = find(p.ctype == 101); bwah2 = find(p.ctype2 == 101); plot(bwah1,p.cpsize(bwah1),bwah2,p.cpsize2(bwah2))
  scatter(p.rlon(bwah1),p.rlat(bwah1),30,p.cngwat(bwah1))
  scatter(p.rlon(bwah2),p.rlat(bwah2),30,p.cngwat(bwah2))  

%[pafter1,BTx] = smooth_cloud_fields(h,p,tempx);   %% smooth the cloud fields
%figure(2); colormap jet; bwah1 = find(pafter1.ctype == 101); bwah2 = find(pafter1.ctype2 == 101); plot(bwah1,pafter1.cpsize(bwah1),bwah2,pafter1.cpsize2(bwah2))
%  scatter(pafter1.rlon(bwah1),pafter1.rlat(bwah1),30,pafter1.cngwat(bwah1))
%  scatter(pafter1.rlon(bwah2),pafter1.rlat(bwah2),30,pafter1.cngwat(bwah2))  

[pafter2,BTx] = smooth_cloud_fields_ctype_careful(h,p,tempx);   %% smooth the cloud fields
figure(3); colormap jet; bwah1 = find(pafter2.ctype == 101); bwah2 = find(pafter2.ctype2 == 101); plot(bwah1,pafter2.cpsize(bwah1),bwah2,pafter2.cpsize2(bwah2))
  scatter(pafter2.rlon(bwah1),pafter2.rlat(bwah1),30,pafter2.cngwat(bwah1))
  scatter(pafter2.rlon(bwah2),pafter2.rlat(bwah2),30,pafter2.cngwat(bwah2))  

p = pafter2;  %% should be good

pxx = p;
for loop = 1 : iN
  if iDoSarta > 0
    %% avoid running SARTA the first time, since you've already done it above  
    pxx = p;    
    pxx = sanity_check(pxx);        
    rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
    [hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);
    alltECM = real(rad2bt(hxx.vchan,pxx.rcalc));

    tcal899 = rad2bt(900,pxx.rcalc(i899,:));
    figure(1); clf; scatter(cldforce,max(pxx.cfrac,pxx.cfrac2),10,pxx.rlat,'filled'); colorbar
    xlabel('cldforce = stemp-BTobs900'); ylabel('cfrac'); title('colorbar = rlat')
    xx = -10:100; plot(xx,max(0,tanh((xx-60)/5))); grid
    [tcld1,tcld2] = do_cld1_cld2(pxx,tobs899,tcal899,iBadDCCOnly);    
    %disp('ret N'); pause

    %% now mark which of the cldforc > 55 have tcld1 > tobs900
    tBTD_DCC_bad = 2;
    oops0 = find(cfracXFORC > 0.5);
    oops1 = find(cfracXFORC > 0.5 & ((tcld1 - tobs899) > tBTD_DCC_bad));  %% bad dcc
    oops2 = find(cfracXFORC < 0.5 & ((tcld1 - tobs899) > tBTD_DCC_bad));  %% bad mid-lower trop clouds
    if length(oops1) > 0 | length(oops2) > 0
      rcalc = p.rcalc;
      robs1 = p.robs1;
      calflag = p.calflag;
      fprintf(1,' found %5i of %5i potential "incorrect cloudset" DCC ie where cfracXFORC > 0.5 & (tcld1 - tobs899 > 2) \n',length(oops1),length(oops0));
      fprintf(1,' found        %5i potential "incorrect cloudset" OTHER ie where cfracXFORC < 0.5 & (tcld1 - tobs899 > 2) \n',length(oops2));      
      pxx = move_cld1_cld2(hxx,pxx,oops1,oops2,tobs899,tempx,tBTD_DCC_bad,RH0,RH1);      
      alltECM = real(rad2bt(hxx.vchan,pxx.rcalc));
      p = pxx;   %% HAVE TO DO THIS since we changed some of the cprtop and ptemp
      p.robs1 = robs1;
      p.rcalc = rcalc;
      p.calflag = calflag;
      
      tcal899 = rad2bt(900,pxx.rcalc(i899,:));
      figure(1); clf; scatter(cldforce,max(pxx.cfrac,pxx.cfrac2),10,pxx.rlat,'filled'); colorbar
      xlabel('cldforce = stemp-BTobs900'); ylabel('cfrac'); title('colorbar = rlat')
      xx = -10:100; plot(xx,max(0,tanh((xx-60)/5))); grid
      [tcld1,tcld2] = do_cld1_cld2(pxx,tobs899,tcal899,iBadDCCOnly);
      %disp('ret N (move cld)'); pause
    end

  end
  iDoSarta = +1;
  xalltECM = [alltECM; prof0.stemp; 250+mmw; 250+lff*50];

  %[damin0, ind0] = min(pdist2(alltECM',alltOBS','euclidean'));  toc2a = toc;    %%% ORIG
  [damin, ind] = min(pdist2(xalltECM',xalltOBS','euclidean')); toc2b = toc;    %%% NEW
  fprintf(1,'   after smoothing cloud fields, loop %2i of 5 --> min = %8.6f \n',loop,sum(damin))
  
  %latspread = max(p.rlat)-min(p.rlat);
  %if latspread > 15  %% llooks like completely ascending or descending granule
  %  disp('A/D granule ...)  
  %else
  %  disp('polar granule ...)
  %  [~, ind] = min(pdist2(xalltECM',xalltOBS','euclidean')); toc2c = toc;    %%% NEW  
  %end
  %error('kjslkjgslkjgs 22')

  FOV_ECM = ind;
  figure(1); clf; plot(FOV_ECM); title('mapping indices to BEST cloud calcs')
  figure(2); clf; plot(1:length(p.stemp),tropopauseT,'b',1:length(p.stemp),tropopauseT(ind),'r')
  figure(3); clf; plot(1:length(p.stemp),tropopauseP,'b',1:length(p.stemp),tropopauseP(ind),'r')

  closestBT = alltECM(:,ind);
  fxx       = h.vchan(ch_ind_use);

  figure(2); clf
  plot(h.vchan(ch_ind_use),nanmean(alltOBS'-alltECM'),'b',h.vchan(ch_ind_use),nanstd(alltOBS'-alltECM'),'c--',...
     h.vchan(ch_ind_use),nanmean(alltOBS'-closestBT'),'r',h.vchan(ch_ind_use),nanstd(alltOBS'-closestBT'),'m--')
  xlabel('wavenumber cm-1'); ylabel('BT(K)'); title('from ECM (b/c) and remapped (r/m)'); grid

  disp('  after smoothing cloud fields remapping cloud fields .....')
  p.ctype = p.ctype(FOV_ECM);
  p.cfrac = p.cfrac(FOV_ECM);
  p.cprtop = p.cprtop(FOV_ECM);
  p.cprbot = p.cprbot(FOV_ECM);
  p.cngwat = p.cngwat(FOV_ECM);
  p.cpsize = p.cpsize(FOV_ECM);

  p.ctype2 = p.ctype2(FOV_ECM);
  p.cfrac2 = p.cfrac2(FOV_ECM);
  p.cprtop2 = p.cprtop2(FOV_ECM);
  p.cprbot2 = p.cprbot2(FOV_ECM);
  p.cngwat2 = p.cngwat2(FOV_ECM);
  p.cpsize2 = p.cpsize2(FOV_ECM);

  p.cfrac12 = p.cfrac12(FOV_ECM);
  
  if isfield(p,'clwc'); p.clwc = p.clwc(:,FOV_ECM); end
  if isfield(p,'ciwc'); p.ciwc = p.ciwc(:,FOV_ECM); end
  if isfield(p,'cc');   p.cc   = p.cc(:,FOV_ECM);   end

  if isfield(p,'orig_ctop');  p.orig_ctop  = p.orig_ctop2(FOV_ECM); end
  if isfield(p,'orig_ctop2'); p.orig_ctop2 = p.orig_ctop2(FOV_ECM); end
  if isfield(p,'sarta_lvlODice');   p.sarta_lvlODice =   p.sarta_lvlODice(:,FOV_ECM); end
  if isfield(p,'sarta_lvlODwater'); p.sarta_lvlODwater = p.sarta_lvlODwater(:,FOV_ECM); end

  if isfield(p,'cemis');   p.cemis  = p.cemis(:,FOV_ECM); end
  if isfield(p,'crho');    p.crho   = p.crho(:,FOV_ECM); end
  if isfield(p,'cstemp');  p.cstemp   = p.cstemp(:,FOV_ECM); end
  if isfield(p,'cemis2');  p.cemis2  = p.cemis2(:,FOV_ECM); end
  if isfield(p,'crho2');   p.crho2   = p.crho2(:,FOV_ECM); end
  if isfield(p,'cstemp3'); p.cstemp2 = p.cstemp2(:,FOV_ECM); end

  pxx = p;
end

p.rcalcBestCld = p.rcalc(:,FOV_ECM); %%% <<<< this is assumed best calc after swap, who knows what will happen after
                                     %%% <<<< re-running sarta-cloud

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% all this is fancy plotting and debugging

if iSwapPlot > 0
  figure(1); clf; scatter(p.rlon,p.rlat,30,rad2bt(1231,p.robs1(1291,:))); title('BT1231 obs')
  figure(2); clf; scatter(p.rlon,p.rlat,30,p.cprtop); title('new cprtop'); caxis([0 1000]); colorbar
  figure(3); clf; scatter(p.rlon,p.rlat,30,p.cprtop2); title('new cprtop2'); caxis([0 1000]); colorbar
  figure(4); clf; scatter(p.rlon,p.rlat,30,log10(p.cngwat)); title('new log10cngwat)')
  figure(5); clf; scatter(p.rlon,p.rlat,30,log10(p.cngwat2)); title('new log10(cngwat2)')
  figure(6); clf; scatter(p.rlon,p.rlat,30,p.cfrac); title('new cfrac')
  figure(7); clf; scatter(p.rlon,p.rlat,30,p.cfrac2); title('new cfrac2')

  pp = 0:25:1000; 
  figure(1); clf; 
  donk = find(p.ctype == 101 & p.cngwat > 0); 
  donk0 = find(prof0.ctype == 101 & prof0.cngwat > 0); 
  donknew = find(pnew0.ctype == 101 & pnew0.cngwat > 0); 
  plot(pp,hist(p.cprtop(donk),pp)/sum(hist(p.cprtop(donk),pp)), 'r',...
       pp,hist(prof0.cprtop(donk0),pp)/sum(hist(prof0.cprtop(donk0),pp)),'b',...
       pp,hist(pnew0.cprtop(donknew),pp)/sum(hist(pnew0.cprtop(donknew),pp)),'k')
  figure(2); clf; 
  donk = find(p.ctype2 == 101 & p.cngwat2 > 0); 
  donk0 = find(prof0.ctype2 == 101 & prof0.cngwat2 > 0); 
  donknew = find(pnew0.ctype2 == 101 & pnew0.cngwat2 > 0); 
  plot(pp,hist(p.cprtop2(donk),pp)/sum(hist(p.cprtop2(donk),pp)), 'r',...
       pp,hist(prof0.cprtop2(donk0),pp)/sum(hist(prof0.cprtop2(donk0),pp)),'b',...
       pp,hist(pnew0.cprtop2(donknew),pp)/sum(hist(pnew0.cprtop2(donknew),pp)),'k')
end

RHF = layeramt2RH(head0,p);

%{
%% when testing hurricane irma, 2017/09/06 g177
addpath /home/sergio/MATLABCODE/COLORMAP
stc = barnet_redgreenblue;
zah = find(p.atrack == 41);
figure(1); pcolor(p.rlon(zah),p.plevs(1:100,1),RH0(:,zah)); shading flat; colorbar; set(gca,'yscale','log'); set(gca,'ydir','reverse');
ax = axis; axis([ax(1) ax(2) 10 1000]); colormap(stc); caxis([0 120])

figure(2); pcolor(p.rlon(zah),p.plevs(1:100,1),RHF(:,zah)); shading flat; colorbar; set(gca,'yscale','log'); set(gca,'ydir','reverse');
ax = axis; axis([ax(1) ax(2) 10 1000]); colormap(stc); caxis([0 120])
error('lkjdg')
%}

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDebug = -1;
if iDebug > 0
  figure(2); clf
  plot(h.vchan(ch_ind_use),nanmean(alltOBS'-alltECM'),'b',h.vchan(ch_ind_use),nanstd(alltOBS'-alltECM'),'c--',...
     h.vchan(ch_ind_use),nanmean(alltOBS'-closestBT'),'k',h.vchan(ch_ind_use),nanstd(alltOBS'-closestBT'),'k--',...
     h.vchan(ch_ind_use),nanmean(alltOBS'-BTx.allchanscprtopfix(ch_ind_use,:)'),'r',...
     h.vchan(ch_ind_use),nanstd(alltOBS'-BTx.allchanscprtopfix(ch_ind_use,:)'),'m--','linewidth',2)
  xlabel('wavenumber cm-1'); ylabel('BT(K)'); title('from ECM (b/c) and remapped (r/m)'); grid

  sigh = find(abs(BTx.tobs-BTx.cprtopfix) == max(abs(BTx.tobs-BTx.cprtopfix)));
  [BTx.tobs(sigh)-BTx.cprtopfix(sigh) BTx.tobs(sigh)]
  [p.cprtop(sigh) p.cngwat(sigh) p.ctype(sigh) prof0.cprtop(sigh) prof0.cngwat(sigh) prof0.ctype(sigh)]
  [p.cprtop2(sigh) p.cngwat2(sigh) p.ctype2(sigh) prof0.cprtop2(sigh) prof0.cngwat2(sigh) prof0.ctype2(sigh)]
  [p.cfrac(sigh) p.cfrac2(sigh) prof0.cfrac(sigh) prof0.cfrac2(sigh) single(p.nlevs(sigh)) single(p.nlevs(ind(sigh)))]
  oo = find(h.vchan(ch_ind_use) >= 1231,1);
  [alltOBS(oo,sigh) alltECM(oo,ind(sigh)) BTx.cprtopfix(sigh)]
  plot(hxx.vchan,alltOBS(:,sigh),'k',hxx.vchan,alltECM(:,ind(sigh)),'b',hxx.vchan,BTx.allchanscprtopfix(ch_ind_use,sigh),'r')
  plot(p.ptemp(:,sigh),1:101,p.ptemp(:,ind(sigh)),1:101,'r'); set(gca,'ydir','reverse'); axis([180 300 1 101]); grid
  semilogx(p.gas_1(:,sigh),1:101,p.gas_1(:,ind(sigh)),1:101,'r'); set(gca,'ydir','reverse'); axis([1e14 1e23 1 101]); grid

  [hxx,pxx] = subset_rtp_allcloudfields(head0,prof0,[],[],sigh);
    rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
    [hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);  
    testprof0.robs = pxx.robs1;
    testprof0.rcalc0 = pxx.rcalc;
  [hxx,pxx] = subset_rtp_allcloudfields(head0,prof0,[],[],ind(sigh));
    rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
    [hxx,haxx,pxxS,paxx] = rtpread(tempx.tempfile3);  
    testprof0.rcalcS = pxxS.rcalc;
  [hxx,pxx] = subset_rtp_allcloudfields(h,p,[],[],sigh);
    rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
    [hxx,haxx,pxxF,paxx] = rtpread(tempx.tempfile3);  
    testprof0.rcalcF = pxxF.rcalc;

  %% compare profile with which we swapped cloud fields (pxxS) and what we want (pxxF)
  [pxxS.stemp pxxF.stemp pxxS.cfrac pxxF.cfrac pxxS.cfrac2 pxxF.cfrac2 pxxS.cfrac12 pxxF.cfrac12]
  [pxxS.cprtop  pxxF.cprtop  pxxS.cprbot pxxF.cprbot   pxxS.cngwat pxxF.cngwat   pxxS.cpsize pxxF.cpsize]
  [pxxS.cprtop2 pxxF.cprtop2 pxxS.cprbot2 pxxF.cprbot2 pxxS.cngwat2 pxxF.cngwat2 pxxS.cpsize2 pxxF.cpsize2]
  plot(hxx.vchan(g),rad2bt(hxx.vchan(g),testprof0.robs(g)),'k',...
       hxx.vchan(g),rad2bt(hxx.vchan(g),testprof0.rcalc0(g)),'b',...
       hxx.vchan(g),rad2bt(hxx.vchan(g),testprof0.rcalcS(g)),'g',...     
       hxx.vchan(g),rad2bt(hxx.vchan(g),testprof0.rcalcF(g)),'r')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
