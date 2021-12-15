function p = fracB1B2(h,p0,frac1,frac2,oo_0_6,tB1_0_6,tB2_0_6,iERAorECM);

if nargin ==  7
  iERAorECM = +1; %% default
end

if abs(iERAorECM) ~= 1
  iERAorECM
  error('need ERA or ECM');
end

%% this is to add profiles from two ERA files together, linearly weighted
%{    
%% from AIRS L1 files
       rlon: [1x384 double]
        rlat: [1x384 double]
       rtime: [1x384 double]
    landfrac: [1x384 double]
       salti: [1x384 double]
        zobs: [1x384 double]
     scanang: [1x384 double]
      satzen: [1x384 double]
      solzen: [1x384 double]
      upwell: [1x384 double]
        udef: [1x384 single]

%% from ERA/ECM
      wspeed: [1x384 double]
       spres: [1x384 single]
       stemp: [1x384 single]
     wsource: [1x384 single]
         tcc: [1x384 single]
        plat: [1x384 single]
        plon: [1x384 single]
       ptemp: [60x384 single]
       gas_1: [60x384 single]
       gas_3: [60x384 single]
          cc: [60x384 single]
        clwc: [60x384 single]
        ciwc: [60x384 single]
       plevs: [60x384 single]
       nlevs: [1x384 int32]
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = p0;
if length(oo_0_6) > 0
  miaow = oo_0_6;
  [hjunk,pjunk] = subset_rtp(h,p,[],[],miaow);
  pB1 = pjunk; pB1.rtime = ones(size(pB1.rtime)) .* tB1_0_6(miaow);
  if iERAorECM == +1
    [pB1,hB1] = fill_era(pB1,h);
  else
    [pB1,hB1] = fill_ecmwf(pB1,h);
  end

  pB2 = pjunk; pB2.rtime = ones(size(pB1.rtime)) .* tB2_0_6(miaow);      
  if iERAorECM == +1
    [pB2,hB2] = fill_era(pB2,h);
  else
    [pB2,hB2] = fill_ecmwf(pB2,h);
  end
  
  %p.sst   = frac1 .* pB1.sst + frac2 .* pB2.sst;
  p.spres(miaow)   = frac1(miaow) .* pB1.spres   + frac2(miaow) .* pB2.spres;    
  p.stemp(miaow)   = frac1(miaow) .* pB1.stemp   + frac2(miaow) .* pB2.stemp;    
  p.wspeed(miaow)  = frac1(miaow) .* pB1.wspeed  + frac2(miaow) .* pB2.wspeed;    
  p.wsource(miaow) = frac1(miaow) .* pB1.wsource + frac2(miaow) .* pB2.wsource;    
  p.tcc(miaow)     = frac1(miaow) .* pB1.tcc     + frac2(miaow) .* pB2.tcc;    
  p.plat(miaow)    = frac1(miaow) .* pB1.plat    + frac2(miaow) .* pB2.plat;    
  p.plon(miaow)    = frac1(miaow) .* pB1.plon    + frac2(miaow) .* pB2.plon;    

  frac1matr = ones(mean(pB1.nlevs),1) * frac1(miaow);
  frac2matr = ones(mean(pB1.nlevs),1) * frac2(miaow);        
  p.ptemp(:,miaow) = frac1matr .* pB1.ptemp + frac2matr .* pB2.ptemp;        
  p.gas_1(:,miaow) = frac1matr .* pB1.gas_1 + frac2matr .* pB2.gas_1;
  p.gas_3(:,miaow) = frac1matr .* pB1.gas_3 + frac2matr .* pB2.gas_3;
  p.cc(:,miaow)    = frac1matr .* pB1.cc    + frac2matr .* pB2.cc;
  p.clwc(:,miaow)  = frac1matr .* pB1.clwc  + frac2matr .* pB2.clwc;
  p.ciwc(:,miaow)  = frac1matr .* pB1.ciwc  + frac2matr .* pB2.ciwc;            
  p.plevs(:,miaow) = frac1matr .* pB1.plevs + frac2matr .* pB2.plevs;

  %p.nlevs = frac1 .* pB1.nlevs + frac2 .* pB2.nlevs;
  %p.nlevs(miaow) = pB1.nlevs(miaow);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
