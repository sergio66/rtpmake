% local running to test
% clustcmd -L jpl_eof_sarta_co2_srf.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 jpl_eof_sarta_co2_srf.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 jpl_eof_sarta_co2_srf.m 001:240

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gg = JOB;    %% granule

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m130_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
code2 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
code3 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m150_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

co2 = [370 385 400];
srf = [130 140 150];

%% randomly subset and then set CO2
data = 1 : 12150;

iCnt = 0;
for ii = 1 : length(co2)
  for jj = 1 : length(srf)
    iCnt = iCnt + 1;
    y = datasample(data,1350,'Replace',false);   %% need to take 1350 unique samples from data
    rand_ind(iCnt,:) = y;
    data = setdiff(data,y);     %% remove y from the data
  end
end
if (sum(unique(rand_ind(:)) - (1:12150)') ~= 0)
  error('hmm, problems in random init')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirIN = ['/asl/data/rtprod_airs/2009/03/01/'];
fIN = ['cloudy_airs_l1b_ecm_sarta_baum_ice.2009.03.01.'];
fIN = [fIN num2str(gg,'%03d') '.rtp'];
fIN = [dirIN fIN];
[h0,ha,p0,pa] = rtpread(fIN);

fip = mktemp('fx.ip.rtp');
fop = mktemp('fx.op.rtp');
frp = mktemp('fx.rp.rtp');
fugh1 = mktemp('fugh1');
fugh2 = mktemp('fugh2');

iCnt = 0;
for ii = 1 : length(co2)
  for jj = 1 : length(srf)

    iCnt = iCnt + 1;
    dx = [num2str(ii) num2str(jj)];
    dirOUT = ['/asl/data/rtprod_airs/2009/03/01/JPL_EOF_' dx '/'];
    fOUT   = ['jpl_' dx '_'  num2str(gg,'%03d') '.rtp'];
    fOUT   = [dirOUT fOUT];

    fprintf(1,'gran %3i  co2loop %3i ppm     srf %3i \n',gg,co2(ii),srf(jj))
    fprintf(1,'  %s \n',fOUT)
    figure(1); plot(1:12150,1:12150,'b.',rand_ind(iCnt,:),rand_ind(iCnt,:),'ro');
    xyz = rand_ind(iCnt,:);
    figure(1); plot(p0.xtrack,p0.atrack,'b.',p0.xtrack(xyz),p0.atrack(xyz),'ro');
    pause(0.1)

    [h,p] = subset_rtp_allcloudfields(h0,p0,[],[],rand_ind(iCnt,:));
    p.co2ppm = ones(size(p.rlat)) * co2(ii);
    if ~isfield(p,'emis')
      [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
    end
    h = rmfield(h,'vchan');
    p = rmfield(p,'robs1');
    p = rmfield(p,'calflag');
    p = rmfield(p,'rcalc');
    h.nchan = 2834;
    h.ichan = (1:2834)';
    rtpwrite(fip,h,ha,p,pa);

    klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' fugh1]; 
    eval(klayerser);

    if jj == 1
      sarta = code1;
    elseif jj == 2
      sarta = code2;
    elseif jj == 3
      sarta = code3;
    end

    sartaer = ['!time ' sarta ' fin=' fop ' fout=' frp ' >& ' fugh2]; 
    eval(sartaer);

    [hx,hax,px,pax] = rtpread(frp);
    rtpwrite(fOUT,hx,hax,px,pax);

%{
    figure(2)
      tcal0 = rad2bt(1231,p.rcalc(1231,:));
      tcalF = rad2bt(1231,px.rcalc(1231,:));
      plot(tcalF-tcal0)

    figure(3)
      tcal0 = rad2bt(h.vchan,p.rcalc);
      tcalF = rad2bt(hx.vchan,px.rcalc);
      plot(h.vchan,tcalF-tcal0)
%}

  end
end

rmer = ['!/bin/rm ' fip ' ' fop  ' '  frp ' ' fugh1 ' ' fugh2];
eval(rmer);