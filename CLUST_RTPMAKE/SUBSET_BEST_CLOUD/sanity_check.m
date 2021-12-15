function prof = sanity_check(pIN0);

prof = pIN0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove cloud fractions less than some minimum

%disp('checking cfrac after smoothing ....')

cmin = 0.00000001;
cmax = 1.00000000;
hcmin = 0.5*cmin;

ix=find(prof.cfrac < cmin);
prof.cfrac(ix)  = 0;
prof.cfrac12(ix)= 0;
prof.cngwat(ix) = 0;
prof.ctype(ix)  = -1;
prof.cprtop(ix) = -9999;
prof.cprbot(ix) = -9999;

ix=find(prof.cfrac2 < cmin);
prof.cfrac2(ix)  = 0;
prof.cfrac12(ix) = 0;
prof.cngwat2(ix) = 0;
prof.ctype2(ix)  = -1;
prof.cprtop2(ix) = -9999;
prof.cprbot2(ix) = -9999;

%{
%%% WOW did not have this in the 
ix=find(prof.ctype == 101 & prof.cprtop < 400);
prof.cfrac(ix)  = 0;
prof.cfrac12(ix) = 0;
prof.cngwat(ix) = 0;
prof.ctype(ix)  = -1;
prof.cprtop(ix) = -9999;
prof.cprbot(ix) = -9999;

ix=find(prof.ctype2 == 101 & prof.cprtop2 < 400);
prof.cfrac2(ix)  = 0;
prof.cfrac12(ix) = 0;
prof.cngwat2(ix) = 0;
prof.ctype2(ix)  = -1;
prof.cprtop2(ix) = -9999;
prof.cprbot2(ix) = -9999;
%}

ix=find(prof.cngwat < 0.00001);
prof.cfrac(ix)  = 0;
prof.cfrac12(ix)= 0;
prof.cngwat(ix) = 0;
prof.ctype(ix)  = -1;
prof.cprtop(ix) = -9999;
prof.cprbot(ix) = -9999;

ix=find(prof.cngwat2 < 0.00001);
prof.cfrac2(ix)  = 0;
prof.cfrac12(ix) = 0;
prof.cngwat2(ix) = 0;
prof.ctype2(ix)  = -1;
prof.cprtop2(ix) = -9999;
prof.cprbot2(ix) = -9999;

ix = find(prof.cfrac12 >= hcmin & prof.cfrac12 < cmin);
prof.cfrac12(ix) = cmin;
ix = find(prof.cfrac12 < hcmin);
prof.cfrac12(ix) = 0;
junk = prof.cfrac(ix) + prof.cfrac2(ix);
ii = ix( find(junk > 1) );
ii1 = ii( find(prof.cfrac(ii) > prof.cfrac2(ii)) );
ii2 = setdiff(ii,ii1);
prof.cfrac(ii1) = prof.cfrac(ii1)-hcmin;
prof.cfrac2(ii2) = prof.cfrac2(ii2)-hcmin;

cfracx  = prof.cfrac - prof.cfrac12;
cfrac2x = prof.cfrac2 - prof.cfrac12;
ctot = cfracx + cfrac2x + prof.cfrac12;
ix = find(ctot > 1);
%% prof.cfrac12(ix) =  max(prof.cfrac(ix) + prof.cfrac2(ix) - 1 + eps,0);
%% prof.cfrac12(ix) =  max(1-cfracx(ix)-cfrac2x(ix) + eps,0);   %% new Jan 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% black cloud check
oo = find(prof.ctype < 100);
prof.cfrac12(oo) = 0.0;
prof.cngwat(oo)  = 0.0;
prof.cfrac(oo)  = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now check cprtop
%disp('---> checking cprtop vs cprbot vs spres')
iNotOK = +1;
iFix = 0;
%% before used to give up at iFix == 10
while iFix < 15 & iNotOK > 0
  iFix = iFix + 1;
  [prof,iNotOK] = check_for_errors(prof,-1,iFix);  %% see possible pitfalls in clouds
   % fprintf(1,' did n=%2i try at checking clouds \n',iFix)
end

if iFix >= 12 & iNotOK > 0
  %disp('oops, could not fix cprtop vs cprbot vs spres')
  %keyboard
  error('oops, could not fix cprtop vs cprbot vs spres')
end
