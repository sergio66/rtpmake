function p2 = driver_sarta_clear_rtp(h,ha,p,pa,run_sarta);

addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/

klayers = run_sarta.klayers;
sarta   = run_sarta.sarta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fip = mktemp('fip.rtp');
fop = mktemp('fop.rtp');
frp = mktemp('frp.rtp');
ugh1 = mktemp('ugh1');
ugh2 = mktemp('ugh2');

if h.ptype == 0
  rtpwrite(fip,h,ha,p,pa);
  klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' > ' ugh1];
  eval(klayerser);
elseif h.ptype == 1
  rtpwrite(fop,h,ha,p,pa);
else
  error('need h.ptype == 0 or 1');
end

sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' > ' ugh2];
eval(sartaer);

[h2,ha2,p2,pa2] = rtpread(frp);

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh1 ' ' ugh2];
eval(rmer);


