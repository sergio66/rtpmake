fprintf(1,'mean p.rtime = %8.6f \n',mean(p.rtime));
[xyyh,xmmh,xddh,xhhh] = tai2utcSergio(mean(p.rtime));
fprintf(1,'mean rtime in AIRS L1B/C data corresponds to following : %4i %2i %2i %6.3f \n',xyyh,xmmh,xddh,xhhh);
%fprintf(1,'when hhh = %6.3f is same as %2i : %2i \n',hhh,floor(xhhh),floor((xhhh-floor(xhhh))*60));

if xhhh >= 0 & xhhh <= 3
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,0.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,3.0);
elseif xhhh > 3 & xhhh <= 6
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,3.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,6.0);      
elseif xhhh > 6 & xhhh <= 9
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,6.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,9.0);      
elseif xhhh > 9 & xhhh <= 12
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,9.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,12.0);
elseif xhhh > 12 & xhhh <= 15
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,12.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,15.0);
elseif xhhh > 15 & xhhh <= 18
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,15.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,18.0);
elseif xhhh > 18 & xhhh <= 21
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,18.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,21.0);
elseif xhhh > 21 & xhhh <= 24
  tB1 = utc2taiSergio(xyyh,xmmh,xddh,21.0);
  tB2 = utc2taiSergio(xyyh,xmmh,xddh,24.0);
end

iFileSpanHour = 3;  %% 8 files so each spans 3 hours
frac1 = 1 - (p.rtime-tB1)/(iFileSpanHour*60*60);
frac1(frac1 < 0) = 0;
frac1(frac1 > 1) = 1;    
frac2 = 1-frac1;

[pClosest,hClosest] = fill_ecmwf(p,h);
pB1 = p; pB1.rtime = ones(size(pB1.rtime)) * tB1;
  [pB1,hB1] = fill_ecmwf(pB1,h);
pB2 = p; pB2.rtime = ones(size(pB1.rtime)) * tB2;      
[pB2,hB2] = fill_ecmwf(pB2,h);

scatter_coast(pClosest.rlon,pClosest.rlat,30,pClosest.stemp-pB1.stemp);
p = pClosest;
p.sst   = frac1 .* pB1.sst + frac2 .* pB2.sst;
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
