function [pjunk,g,FOV_ECM,fxx,alltOBS,alltECM,closestBT] = subset_best_cloud(h,ha,xp,pa,tempx)

%see /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/loop_get_bestL2orECM_FAST_5deglat.m

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD

g = dogoodchan;

%{
%% technically I should do this since I want the list of non-popped channels for this granule
%% bit look below, I end up using cind(1:5) so am using a default set of window channels ....
g = findgoodchan(h,ha,xp,pa);;
%}


%{
iStemp_ColWV = 5;
ch_ind_xyz = get_retrieval_chans(h,g,iStemp_ColWV);

ch_ind_use = g;
ch_ind_use = ch_ind_xyz;

%ch_ind_use = ch_ind_use(ch_ind_use < 1620);        %% keep to LW and MW
zooLW = find(ch_ind_use >=  246 & ch_ind_use <= 900);   %% keep to 720  cm-1 < ch < 960 cm-1
zooMW = find(ch_ind_use >= 1431 & ch_ind_use <= 1864);  %% keep to 1320 cm-1 < ch < 1620 cm-1
oo = union(zooLW,zooMW(1:3:end));

%% use only window chans for this cloud stuff
oo = find((h.vchan(ch_ind_xyz) > 743 & h.vchan(ch_ind_xyz) < 970) | (h.vchan(ch_ind_xyz) > 1180 & h.vchan(ch_ind_xyz) < 1280));
%}

addpath /home/sergio/MATLABCODE/DUSTFLAG
freqsNchans;
ch_ind_use = intersect(g,cind1(1:5));         
ch_ind_use = cind1(1:5);;  %%% <<<<<<<<<<<<<<< we assume these are the chans being used in loop_get_bestL2orECM_FAST_5deglat.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now run through klayers and sarta for select chans
%% then run through the cloud swapper
xp.robsqual = xp.robs1;
[hxx,pxx] = subset_rtp_allcloudfields(h,xp,[],ch_ind_use,[]);

iprtp = mktemp('junk.ip.rtp');
oprtp = mktemp('junk.op.rtp');
rprtp = mktemp('junk.rp.rtp');
uprtp = mktemp('junkugh');
rtpwrite(iprtp,hxx,ha,pxx,pa);
klayerser = ['!time ' tempx.klayers   ' fin=' iprtp ' fout=' oprtp ' >& ' uprtp]; eval(klayerser)
sartaer   = ['!time ' tempx.sarta     ' fin=' oprtp ' fout=' rprtp];              eval(sartaer)

[hxx,haxx,pxx,paxx] = rtpread(rprtp);
rmer = ['!/bin/rm ' iprtp ' ' oprtp ' ' rprtp ' ' uprtp]; eval(rmer);
tobsx = rad2bt(hxx.vchan(5),pxx.robs1(5,:));
tcalx = rad2bt(hxx.vchan(5),pxx.rcalc(5,:));
plot(tobsx-tcalx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempx.tempfile1 = iprtp;
tempx.tempfile2 = oprtp;
tempx.tempfile3 = rprtp;
tempx.tempfile4 = uprtp;

xp.rcalc(hxx.ichan,:) = pxx.rcalc;
xp.robsqual(hxx.ichan,:) = pxx.rcalc;

pjunk = pxx;
pjunk.rcalc = zeros(2378,length(pxx.stemp));
pjunk.robs1 = zeros(2378,length(pxx.stemp));
pjunk.robs1 = xp.robs1;
pjunk.rcalc(hxx.ichan,:) = pxx.rcalc;
pjunk.robsqual = pjunk.rcalc;
pjunk.calflag = xp.calflag;

hjunk = h;
hjunk.pfields = hjunk.pfields + 2;
  hjunk.ngas = hxx.ngas;
  hjunk.glist = hxx.glist;
  hjunk.gunit = hxx.gunit;
  hjunk.ptype = 1;
[pjunk,FOV_ECM,fxx,alltOBS,alltECM,closestBT] = loop_get_bestL2orECM_FAST_5deglat(hjunk,pjunk,g,0,4,tempx,0);