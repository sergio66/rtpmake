sarta = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
sartaer = ['!' sarta ' fin=' fop ' fout=' frp];

if isfield(globalavg,'cloudxtra')
  cloudaxtra = globalavg.cloudxtra;
  globalavg = rmfield(globalavg,'cloudxtra');
end

globalclr = globalavg;
globalclr.ctype = -9999 * ones(size(globalclr.stemp));
globalclr.cfrac = 0 * ones(size(globalclr.stemp));
globalclr.cngwat = 0 * ones(size(globalclr.stemp));
globalclr.ctype2 = -9999 * ones(size(globalclr.stemp));
globalclr.cfrac2 = 0 * ones(size(globalclr.stemp));
globalclr.cngwat2 = 0 * ones(size(globalclr.stemp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pjunk = globalavg;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
bt0 = rad2bt(hjunk.vchan,pjunk.rcalc);

pjunk = globalavg;
pjunk.stemp = pjunk.stemp + 1;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
btST = rad2bt(hjunk.vchan,pjunk.rcalc);

pjunk = globalavg;
pjunk.ptemp = pjunk.ptemp + 1;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
btT = rad2bt(hjunk.vchan,pjunk.rcalc);

pjunk = globalavg;
pjunk.gas_1 = pjunk.gas_1 * 1.1;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
btWV = rad2bt(hjunk.vchan,pjunk.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pjunk = globalclr;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
clr_bt0 = rad2bt(hjunk.vchan,pjunk.rcalc);

pjunk = globalclr;
pjunk.stemp = pjunk.stemp + 1;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
clr_btST = rad2bt(hjunk.vchan,pjunk.rcalc);

pjunk = globalclr;
pjunk.ptemp = pjunk.ptemp + 1;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
clr_btT = rad2bt(hjunk.vchan,pjunk.rcalc);

pjunk = globalclr;
pjunk.gas_1 = pjunk.gas_1 * 1.1;
rtpwrite(fop,hnew_op,ha2,pjunk,pa2);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
clr_btWV = rad2bt(hjunk.vchan,pjunk.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9); plot(hjunk.vchan,clr_btST-clr_bt0,'b',hjunk.vchan,btST-bt0,'r'); xlim([650 1650])
  title('Compare (b)clr (r) cld jacs for STEMP')
figure(10); plot(hjunk.vchan,clr_btT-clr_bt0,'b',hjunk.vchan,btT-bt0,'r'); xlim([650 1650])
  title('Compare (b)clr (r) cld jacs for COLTEMP')
figure(11); plot(hjunk.vchan,(clr_btWV-clr_bt0)/log(1.1),'b',hjunk.vchan,(btWV-bt0)/log(1.1),'r'); xlim([650 1650])
  title('Compare (b)clr (r) cld jacs for COLWV')
figure(12); plot(hjunk.vchan,1/15*(clr_btWV-clr_bt0)/log(1.1) + clr_btT-clr_bt0,'b',hjunk.vchan,1/15*(btWV-bt0)/log(1.1) + btT-bt0,'r'); xlim([650 1650])
  title('Compare (b)clr (r) cld jacs for COLT + 1/15*COLWV')

jacWV = (btWV-bt0)/log(1.1);
jacT  = (btT-bt0);
jacST = (btST-bt0);

clr_jacWV = (clr_btWV-clr_bt0)/log(1.1);
clr_jacT  = (clr_btT-clr_bt0);
clr_jacST = (clr_btST-clr_bt0);

ix = 10; %% 95 pctile
ix = 08; %% 80 pctile
figure(13); plot(hjunk.vchan,0.0025*jacWV(:,ix) + 0.025*jacT(:,ix) + 0.02*jacST(:,ix),'b',hjunk.vchan,0.0025*clr_jacWV(:,ix) + 0.025*clr_jacT(:,ix) + 0.02*clr_jacST(:,ix),'r'); xlim([650 1650])
  title('Compare (b)clr (r) cld jacs for COLT + 0.0025*COLWV')
