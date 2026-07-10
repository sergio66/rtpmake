function [p,bady0] = fix_nan_emis(p0);

p = p0;
bady0 = [];

[mm,nn] = size(p.emis);
[badx,bady] = find(isnan(p.emis) | isinf(p.emis));
bady = unique(bady);

if length(bady) > 0
  bady0 = bady;
  fprintf(1,'round 1 : found %5i bad emiss \n',length(bady));
  %printarray(bady);
  for pppps = 1 : length(bady)
    pppp = bady(pppps);
    nemiss = p.nemis(pppp);
    efreqs = p.efreq(:,pppp);    
    emisss = p.emis(:,pppp);
    rhosss = p.rho(:,pppp);
    efreqs(nemiss+1:mm) = 0;
    emisss(nemiss+1:mm) = 0;
    rhosss(nemiss+1:mm) = 0;
    p.efreq(:,pppp) = efreqs;
    p.emis(:,pppp)  = emisss;
    p.rho(:,pppp)   = rhosss;
  end
end

[badx,bady] = find(isnan(p.emis) | isinf(p.emis));
bady = unique(bady);
if length(bady) > 0
  fprintf(1,'round 2 : still found %5i bad emiss \n',length(bady));
  for pppps = 1 : length(bady)
    pppp = bady(pppps);
    nemiss = p.nemis(pppp);
    efreqs = p.efreq(1:nemiss,pppp);    
    emisss = p.emis(1:nemiss,pppp);
    rhosss = p.rho(1:nemiss,pppp);
    good = isfinite(emisss) & isfinite(rhosss);
    bad  = ~isfinite(emisss) |  ~isfinite(rhosss);    
    emisss(bad) = interp1(efreqs(good),emisss(good),efreqs(bad),[],'extrap');
    rhosss(bad) = (1-emisss(bad))/pi;
    p.efreq(1:nemiss,pppp) = efreqs;
    p.emis(1:nemiss,pppp)  = emisss;
    p.rho(1:nemiss,pppp)   = rhosss;
  end
end

[badx,bady] = find(isnan(p.emis) | isinf(p.emis));
bady = unique(bady);
if length(bady) > 0
  fprintf(1,'round 3 : oh oh still found %5i bad emiss \n',length(bady));
  disp('now what????')
end

bad = find(p.emis > 1 | p.emis < 0 | p.rho < 0);
if length(bad) > 0
  fprintf(1,'round F0 : checking 0 <= emis <= 1 found %5i bad emiss which is %8.4f percent --- clipping \n',length(bad),length(bad)*100/(mm*nn));  
  p.emis(p.emis > 1) = 1.0;
  p.emis(p.emis < 0) = 0.0;
  p.rho = (1-p.emis)/pi;
end

bad = find(p.emis > 1 | p.emis < 0 | p.rho < 0);
if length(bad) > 0
  fprintf(1,'round F1 : checking 0 <= emis <= 1 found %5i bad emiss which is %8.4f percent --- clipping \n',length(bad),length(bad)*100/(mm*nn));  
  p.emis(p.emis > 1) = 1.0;
  p.emis(p.emis < 0) = 0.0;
  p.rho = (1-p.emis)/pi;
end
