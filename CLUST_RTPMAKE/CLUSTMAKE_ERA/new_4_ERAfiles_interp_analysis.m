[xyyh,xmmh,xddh,xhhh] = tai2utcSergio(mean(p.rtime));
fprintf(1,'mean rtime in AIRS L1B/C data corresponds to following : %4i %2i %2i %6.3f \n',xyyh,xmmh,xddh,xhhh);

[xyyh,xmmh,xddh,xhhh] = tai2utcSergio(p.rtime);
daysSince2002 = change2days(xyyh,xmmh,xddh,2002);

numdays = max(daysSince2002)-min(daysSince2002);
if numdays > 1
  disp('warning this code can only handle a 1 day spread in rtime eg 2012/02/24 h08 to 2012/03/24 h07')
end

numdays = max(p.rtime) - min(p.rtime);
numdays = numdays/60/60/24;
if numdays > 1
  disp('warning this code can only handle a 1 day spread in rtime eg 2012/02/24 h08 to 2012/03/24 h07')
end

booX = nan(size(p.rtime));
boo = find(xhhh >= 0 & xhhh <= 6);
  tB1(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),0.0);
  tB2(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),6.0);
  booX(boo) = 1;
boo = find(xhhh > 6 & xhhh <= 12);
  tB1(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),6.0);
  tB2(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),12.0);      
  booX(boo) = 2;
boo = find(xhhh > 12 & xhhh <= 18);
  tB1(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),12.0);
  tB2(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),18.0);      
  booX(boo) = 3;
boo = find(xhhh > 18 & xhhh <= 24);
  tB1(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),18.0);
  tB2(boo) = utc2taiSergio(xyyh(boo),xmmh(boo),xddh(boo),24.0);
  booX(boo) = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFileSpanHour = 6;  %% 4 files so each spans 6 hours
frac1 = 1 - (p.rtime-tB1)/(iFileSpanHour*60*60);
frac1(frac1 < 0) = 0;
frac1(frac1 > 1) = 1;    
frac2 = 1-frac1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pClosest,hClosest] = fill_era(p,h);
pB1 = p; pB1.rtime = ones(size(pB1.rtime)) .* tB1;
  [pB1,hB1] = fill_era(pB1,h);
pB2 = p; pB2.rtime = ones(size(pB1.rtime)) .* tB2;      
[pB2,hB2] = fill_era(pB2,h);

scatter_coast(pClosest.rlon,pClosest.rlat,30,pClosest.stemp-pB1.stemp);
p = pClosest;

%p.sst   = frac1 .* pB1.sst + frac2 .* pB2.sst;
p.spres = frac1 .* pB1.spres + frac2 .* pB2.spres;    
p.stemp = frac1 .* pB1.stemp + frac2 .* pB2.stemp;
p.wspeed = frac1 .* pB1.wspeed + frac2 .* pB2.wspeed;
p.wsource = frac1 .* pB1.wsource + frac2 .* pB2.wsource;
p.tcc = frac1 .* pB1.tcc + frac2 .* pB2.tcc;
p.plat = frac1 .* pB1.plat + frac2 .* pB2.plat;
p.plon = frac1 .* pB1.plon + frac2 .* pB2.plon;

frac1matr = ones(mean(pB1.nlevs),1) * frac1;
frac2matr = ones(mean(pB1.nlevs),1) * frac2;        
p.ptemp = frac1matr .* pB1.ptemp + frac2matr .* pB2.ptemp;        
p.gas_1 = frac1matr .* pB1.gas_1 + frac2matr .* pB2.gas_1;
p.gas_3 = frac1matr .* pB1.gas_3 + frac2matr .* pB2.gas_3;
p.cc = frac1matr .* pB1.cc + frac2matr .* pB2.cc;
p.clwc = frac1matr .* pB1.clwc + frac2matr .* pB2.clwc;
p.ciwc = frac1matr .* pB1.ciwc + frac2matr .* pB2.ciwc;            
p.plevs = frac1matr .* pB1.plevs + frac2matr .* pB2.plevs;

%p.nlevs = frac1 .* pB1.nlevs + frac2 .* pB2.nlevs;
p.nlevs = pB1.nlevs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.ngas = 2;
h.gunit = [21 21]';
h.glist = [01 03]';
h.ptype = 0;
