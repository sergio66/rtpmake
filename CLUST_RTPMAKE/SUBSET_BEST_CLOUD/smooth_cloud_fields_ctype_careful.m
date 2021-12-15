function [p,BTx] = smooth_cloud_fields(h,pIN0,tempx)

%% note that pIN and pBestCld have EXACTLY the same fields, except "rcalc" where pBestCld has the swapped "best" rcalc estimate,
%% which is not necessarily true (since no geo profile included)

pIN = sanity_check(pIN0);

pIN = sanity_check(pIN);

pOrigSwapped = pIN;
p = pIN;

freqsNchans;
woo = find(h.vchan >= 1231,1);
woo = find(h.vchan >= 899,1);

woo = 758; %% 900 cm-1
[Y,wooA,wooB] = intersect(h.ichan,758); woo = wooA;

if length(woo) == 0
  error('oh oh no channel 758 (900 cm-1)')
end

BTx.tobs    = rad2bt(h.vchan(woo),pIN.robs1(woo,:));         %% BT1231 obs
BTx.tcalc0  = rad2bt(h.vchan(woo),pIN.rcalc(woo,:));         %% BT1231 calc, straight from  ~sergio/MATLABCODE/matlib/clouds/sarta/driver_sarta_cloud
BTx.pseudobest = rad2bt(h.vchan(woo),p.rcalcBestCld(woo,:)); %% BT1231 calc, after swapping in best cloud calc, not necessarily true (since no geo profile included)

[hxx,pxx] = subset_rtp(h,pIN,[],woo,[]);
if hxx.ptype == 0
  rtpwrite(tempx.tempfile1,hxx,[],pxx,[]);
  klayerser = ['!' tempx.klayers ' fin=' tempx.tempfile1 ' fout=' tempx.tempfile2 ' >& ' tempx.tempfile4]; eval(klayerser)
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
else
  rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
end
[hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);
BTx.truebest = rad2bt(hxx.vchan,pxx.rcalc); %% BT1231 calc, after swapping in best cloud calc, true (geo profile included)

if length(p.cfrac) == 12150 | length(p.cfrac) == 1350
  disp('smoothing new ice cloud fields')
  
  thetype1 = p.ctype;
  thetype1 = reshape(thetype1,90,135);
  thetype2 = p.ctype2;
  thetype2 = reshape(thetype2,90,135);

  %% (y)es this is EQUAL to
  ywah1_101 = find(thetype1 == 101);   ywah2_101 = find(thetype2 == 101);
  ywah1_201 = find(thetype1 == 201);   ywah2_201 = find(thetype2 == 201);
  
  %% (n)o this is NOT EQUAL to
  nwah1_101 = find(thetype1 ~= 101);   nwah2_101 = find(thetype2 ~= 101);
  nwah1_201 = find(thetype1 ~= 201);   nwah2_201 = find(thetype2 ~= 201);

  clear x x0
  x0 = double(p.cfrac); 
    figure(7); clf; scatter(p.rlon,p.rlat,30,x0);
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
    x = single(x);
    x = x(:); x(x<0) = 0; x(x>1) = 1; x = x';
    zz = find(p.cfrac == 0); x(zz) = 0; 
    figure(8); clf; scatter(p.rlon,p.rlat,30,x);
    p.cfrac = x; %disp('ret to continue'); pause

  clear x x0
  x0 = double(p.cfrac2); 
    figure(7); clf; scatter(p.rlon,p.rlat,30,x0);
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);
    x = single(x);    
    x = x(:); x(x<0) = 0; x(x>1) = 1; x = x';
    zz = find(p.cfrac2 == 0); x(zz) = 0; 
    figure(8); clf; scatter(p.rlon,p.rlat,30,x); p.cfrac2 = x; %disp('ret to continue'); pause

  clear x x0
  x0 = p.cfrac12; figure(7); clf; scatter(p.rlon,p.rlat,30,x0);
    x0 = reshape(x0,90,135);
    x0 = smoothn(x0);
    x = x0;
    x = x(:);  x(x<0) = 0; x(x>1) = 1; x = x';
    zz = find(p.cfrac12 == 0); x(zz) = 0; 
    figure(8); clf; scatter(p.rlon,p.rlat,30,x); p.cfrac12 = x; %disp('ret to continue'); pause

  clear x x0
  x0 = double(p.cprtop);
    figure(7); clf; scatter(p.rlon,p.rlat,30,x0);
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
    x = single(x);        
    x = x(:); x(x < 0) = p.cprtop(x<0); 
    bad = find(x > p.spres'-50); x(bad) = p.spres(bad)-60; x = x';
    figure(8); clf; scatter(p.rlon,p.rlat,30,x); p.cprtop = x; p.cprbot = p.cprtop+50;  
    figure(7); caxis([0 1000]); colorbar; figure(8); caxis([0 1000]); colorbar; %disp('ret to continue'); pause

  clear x x0
  x0 = double(p.cngwat);
    figure(7); clf; scatter(p.rlon,p.rlat,30,log10(x0));
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
    x = single(x);        
    x = x(:); x(x < 0) = 0; x = x';
    zz = find(p.cngwat == 0); x(zz) = 0; 
    figure(8); clf; scatter(p.rlon,p.rlat,30,log10(x)); p.cngwat = x; %disp('ret to continue'); pause

  clear x x0
  x0 = double(p.cpsize);
    figure(7); clf; scatter(p.rlon,p.rlat,30,x0);
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
    x = single(x);        
    x = x(:); x(x < 0) = 0; x = x';
    figure(8); clf; scatter(p.rlon,p.rlat,30,x); p.cpsize = x; %disp('ret to continue'); pause

  clear x x0
  x0 = double(p.cprtop2);
    figure(7); clf; scatter(p.rlon,p.rlat,30,x0);
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);    
    x = x(:); x(x < 0) = p.cprtop2(x<0);
    x = single(x);        
    bad = find(x > p.spres'-50); x(bad) = p.spres(bad)-60; x = x'; 
    figure(8); clf; scatter(p.rlon,p.rlat,30,x); p.cprtop2 = x; p.cprbot2 = p.cprtop2+50;
    figure(7); caxis([0 1000]); colorbar; figure(8); caxis([0 1000]); colorbar; %disp('ret to continue'); pause

  clear x x0
  x0 = double(p.cngwat2);
    figure(7); clf; scatter(p.rlon,p.rlat,30,log10(x0));
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);
    x = single(x);        
    x = x(:); x(x < 0) = 0; x = x';
    zz = find(p.cngwat2 == 0); x(zz) = 0; 
    figure(8); clf; scatter(p.rlon,p.rlat,30,log10(x)); p.cngwat2 = x; %disp('ret to continue'); pause

  clear x x0
  x0 = double(p.cpsize2);
    figure(7); clf; scatter(p.rlon,p.rlat,30,x0);
    x0 = reshape(x0,90,135);
    xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
    xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
    x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);
    x = single(x);        
    x = x(:); x(x < 0) = 0; x = x';
    figure(8); clf; scatter(p.rlon,p.rlat,30,x); p.cpsize2 = x; %disp('ret to continue'); pause

elseif length(p.atrack) > 10000 & mod(length(p.atrack),12150) == 0
  %% ie we are sure N grans have been stuck together
  boo = unique(p.findex);
  for bboo = 1 : length(boo)
    aha = find(p.findex == boo(bboo));
    
    inds = (p.atrack(aha)-1)*90 + p.xtrack(aha);
    ninds = setdiff(1:12150,inds);     %% make sure no problems if less than 12150 profiles to map to 90x135 and then smooth!!

    clear x x0
    x0 = zeros(1,12150);
    x0(inds) = p.cfrac(aha); x0(ninds) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,x0(inds));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
      x = x(:); x(x<0) = 0; x(x>1) = 1; x = x';
      zz = find(p.cfrac(aha) == 0); x(inds(zz)) = 0; 
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,x(inds)); p.cfrac(aha) = x(inds); %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);
    x0(inds) = p.cfrac2(aha); x0(ninds) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,x0(inds));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);    
      x = x(:); x(x<0) = 0; x(x>1) = 1; x = x';
      zz = find(p.cfrac2(aha) == 0); x(inds(zz)) = 0; 
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,x(inds)); p.cfrac2(aha) = x(inds); %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);
    x0(inds) = p.cfrac12(aha); x0(ninds) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,x0(inds));
      x0 = reshape(x0,90,135);
      x = x0;
      x = smoothn(x);  
      x = x(:); x(x<0) = 0; x(x>1) = 1; x = x';
      zz = find(p.cfrac12(aha) == 0); x(inds(zz)) = 0; 
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,x(inds)); p.cfrac12(aha) = x(inds); %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);
    x0(inds) = p.cprtop(aha); x0(ninds) = NaN; x0(x0<0) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,x0(inds));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
      x = x(:); x(x < 0) = p.cprtop(aha(x<0)); 
      bad = find(x(inds) > p.spres(aha)'-50); x(bad) = p.spres(aha(bad))-60; x = x';
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,x(inds)); p.cprtop(aha) = x(inds); p.cprbot(aha) = p.cprtop(aha)+50;  
      figure(7); caxis([0 1000]); colorbar; figure(8); caxis([0 1000]); colorbar; %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);
    x0(inds) = p.cprtop2(aha); x0(ninds) = NaN; x0(x0<0) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,x0(inds));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);    
      x = x(:); x(x < 0) = p.cprtop2(aha(x<0)); 
      bad = find(x(inds) > p.spres(aha)'-50); x(bad) = p.spres(aha(bad))-60; x = x'; 
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,x(inds)); p.cprtop2(aha) = x(inds); p.cprbot2(aha) = p.cprtop2(aha)+50;
      figure(7); caxis([0 1000]); colorbar; figure(8); caxis([0 1000]); colorbar; %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);    
    x0(inds) = p.cngwat(aha); x0(ninds) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,log10(x0(inds)));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
      x = x(:); x(x < 0) = 0; x = x';
      zz = find(p.cngwat(aha) == 0); x(inds(zz)) = 0; 
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,log10(x(inds))); p.cngwat(aha) = x(inds); %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);    
    x0(inds) = p.cngwat2(aha); x0(ninds) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,log10(x0(inds)));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);    
      x = x(:); x(x < 0) = 0; x = x';
      zz = find(p.cngwat2(aha) == 0); x(inds(zz)) = 0; 
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,log10(x(inds))); p.cngwat2(aha) = x(inds); %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);    
    x0(inds) = p.cpsize(aha); x0(ninds) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,x0(inds));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah1_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah1_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah1_101) = xjunk_101(ywah1_101); x(ywah1_201) = xjunk_201(ywah1_201);
      x = x(:); x(x < 0) = 0; x = x';
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,x(inds)); p.cpsize(aha) = x(inds); %disp('ret to continue'); pause

    clear x x0
    x0 = zeros(1,12150);    
    x0(inds) = p.cpsize2(aha); x0(ninds) = NaN;
      figure(7); clf; scatter(p.rlon(aha),p.rlat(aha),30,x0(inds));
      x0 = reshape(x0,90,135);
      xjunk_101 = x0; xjunk_101(nwah2_101) = NaN; xjunk_101 = smooth2a(xjunk_101,1);
      xjunk_201 = x0; xjunk_201(nwah2_201) = NaN; xjunk_201 = smooth2a(xjunk_201,1);        
      x = x0; x(ywah2_101) = xjunk_101(ywah2_101); x(ywah2_201) = xjunk_201(ywah2_201);    
      x = x(:); x(x < 0) = 0; x = x';
      figure(8); clf; scatter(p.rlon(aha),p.rlat(aha),30,x(inds)); p.cpsize2(aha) = x(inds); %disp('ret to continue'); pause
  end
else
  disp('oops not smoothing since significantly less than 12150 FOVS OR not a multiple of 12150')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(p.cfrac) == 12150 | length(p.cfrac) == 1350
  disp('averaging SMOOTHED and ORIG cloud fields');
  p.cfrac = min(max((p.cfrac + pOrigSwapped.cfrac)/2,0),1);
  p.cfrac2 = min(max((p.cfrac2 + pOrigSwapped.cfrac2)/2,0),1);
  p.cfrac12 = min(max((p.cfrac12 + pOrigSwapped.cfrac12)/2,0),1);
  p.cngwat = min(max((p.cngwat + pOrigSwapped.cngwat)/2,0),300);
  p.cngwat2 = min(max((p.cngwat2 + pOrigSwapped.cngwat2)/2,0),300);
  p.cprtop = min(max((p.cprtop + pOrigSwapped.cprtop)/2,0),1000);
  p.cprtop2 = min(max((p.cprtop2 + pOrigSwapped.cprtop2)/2,0),1000);
end

%% now need to do sanity checks >>>>>>>>>>>>>>>>>>>>>>>>
p = sanity_check(p);

[hxx,pxx] = subset_rtp(h,p,[],woo,[]);
if hxx.ptype == 0
  rtpwrite(tempx.tempfile1,hxx,[],pxx,[]);
  klayerser = ['!' tempx.klayers ' fin=' tempx.tempfile1 ' fout=' tempx.tempfile2 ' >& ' tempx.tempfile4]; eval(klayerser)
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
else
  rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
end
[hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);
BTx.after_smooth = rad2bt(hxx.vchan,pxx.rcalc); %% BT1231 calc, after swapping in best cloud calc, true (geo profile included, after smoothing)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now improve cprtop >>>>>>>>>>>>>>>>>>>>>>
iAdjCPTOP = -1;
if iAdjCPTOP > 0
  p = adjust_cptop(h,p);                   %% if stemp-BT1231obs > 30, adjust the ice cloud field
  p = sanity_check(p);

  [hxx,pxx] = subset_rtp(h,p,[],woo,[]);
  if hxx.ptype == 0
    rtpwrite(tempx.tempfile1,hxx,[],pxx,[]);
    klayerser = ['!' tempx.klayers ' fin=' tempx.tempfile1 ' fout=' tempx.tempfile2 ' >& ' tempx.tempfile4]; eval(klayerser)
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
  else
    rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
  end
end
[hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);  %% whtether or not you ran sarta again, just read in tempfile3
BTx.cprtopfix = rad2bt(hxx.vchan,pxx.rcalc); %% BT1231 calc, after swapping in best cloud calc, true (geo profile included, after smoothing)

figure(8); clf
  plot(1:length(p.stemp),[BTx.tobs; BTx.tcalc0; BTx.pseudobest; BTx.truebest; BTx.after_smooth; BTx.cprtopfix])
  plot(1:length(p.stemp),BTx.tobs,'k',1:length(p.stemp),BTx.tcalc0,'b',1:length(p.stemp), BTx.cprtopfix,'r')
  plot(BTx.tobs,BTx.tcalc0,'b.',BTx.tobs, BTx.cprtopfix,'r.',BTx.tobs,BTx.tobs,'k')
  xlabel('BT 1231 obs'); ylabel('BT1231 calc (b) orig (r) after swap')
figure(9); clf
  plot(1:length(p.stemp),[BTx.tobs-BTx.tcalc0; BTx.tobs-BTx.pseudobest; BTx.tobs-BTx.truebest; BTx.tobs-BTx.after_smooth; BTx.tobs-BTx.cprtopfix],'linewidth',2)
  hl = legend('tcalc0','pseudobest','truebest','after smooth','cprtopfix'); set(hl,'fontsize',10)
figure(10); clf
  dBT = -100 : 5 : +100;
  n0 = hist(BTx.tobs-BTx.tcalc0,dBT);
  np = hist(BTx.tobs-BTx.pseudobest,dBT);  %%% >>> this is always gonna look great, as this comes from swapping in best calc w/o worrying!!!
  nt = hist(BTx.tobs-BTx.truebest,dBT);
  ns = hist(BTx.tobs-BTx.after_smooth,dBT);  
  nc = hist(BTx.tobs-BTx.cprtopfix,dBT);
  plot(dBT,[n0; np; nt; ns; nc],'linewidth',2); grid
  hl = legend('tcalc0','pseudobest','truebest','after smooth','cprtopfix'); set(hl,'fontsize',10)  
  maxx = max(max([n0; nt; ns; nc]));
  axis([-75 +75 0 maxx]); grid on;
  pause(0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iAllChan = -1;
if iAllChan > 0
  hxx = h;
  pxx = p;
  if hxx.ptype == 0
    rtpwrite(tempx.tempfile1,hxx,[],pxx,[]);
    klayerser = ['!' tempx.klayers ' fin=' tempx.tempfile1 ' fout=' tempx.tempfile2 ' >& ' tempx.tempfile4]; eval(klayerser)
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
  else
    rtpwrite(tempx.tempfile2,hxx,[],pxx,[]);
    sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)
  end
  [hxx,haxx,pxx,paxx] = rtpread(tempx.tempfile3);  %% whtether or not you ran sarta again, just read in tempfile3
  BTx.allchanscprtopfix = rad2bt(hxx.vchan,pxx.rcalc); %% BT1231 calc, after swapping in best cloud calc, true (geo profile included, after smoothing)
  BTx.allchans_vchan = hxx.vchan;
end
