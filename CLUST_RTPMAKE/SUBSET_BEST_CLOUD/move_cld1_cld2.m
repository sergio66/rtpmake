function  pyy = move_cld1_cld2(hxx,pxx,oops1,oops2,tobs899,tempx,tBTD_DCC_bad,RH0,RH1);

%% hxx,pxx are rtp structures
%%    oops1 = find(cfracXFORC > 0.5 & ((tcld1 - tobs899) > tBTD_DCC_bad));  %% bad dcc
%%    oops2 = find(cfracXFORC < 0.5 & ((tcld1 - tobs899) > tBTD_DCC_bad));  %% bad mid-lower trop clouds

addpath /home/sergio/MATLABCODE/matlib/clouds/sarta/

pxx0 = pxx.cprtop;

pyy = pxx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(oops1) > 0
  tcld1    = ones(size(pxx.stemp))*NaN;
  ycld1    = ones(size(pxx.stemp))*NaN;
  tcld2    = ones(size(pxx.stemp))*NaN;
  newcldp1 = ones(size(pxx.stemp))*NaN;
  newclds1 = ones(size(pxx.stemp))*NaN;
  
  for iii = 1 : length(oops1)
    ii = oops1(iii);
    nlevs = pxx.nlevs(ii);
    xp = pxx.plevs(1:nlevs,ii);
    xt = pxx.ptemp(1:nlevs,ii);
    xw = pxx.gas_1(1:nlevs,ii);
    const_rh_fac = ones(size(xw));
    if pxx.cfrac(ii) > 0 & pxx.cngwat(ii) > 0
      tcld1(ii) = interp1(log(xp),xt,log(pxx.cprtop(ii)));
    end
    if pxx.cfrac2(ii) > 0 & pxx.cngwat2(ii) > 0
      tcld2(ii) = interp1(log(xp),xt,log(pxx.cprtop2(ii)));
    end

    minp = max(50,pxx.cprtop(ii)-50);
    kkk = find(xp <= pxx.cprtop(ii) & xp > minp);
    ttt = xt(kkk);
    aha = find(ttt <= (tobs899(ii) - tBTD_DCC_bad));

    if length(aha) == 0
      %% now we need to reduce p.ptemp at cprtop so it is 2 K cooler than tobs899
      %% first do linear pert from CPRTOP to 50 mb (or 50 mb above CPRTOP, whichever is larger)
      minp = max(50,pxx.cprtop(ii)-50);
      kkk = find(xp <= pxx.cprtop(ii) & xp >= minp);    
      slope = -(2*tBTD_DCC_bad) / abs(pxx.cprtop(ii) - minp);
      %[slope minp pxx.cprtop(ii)]
      deltaT = slope*(xp(kkk)-(pxx.cprtop(ii))) + -2*tBTD_DCC_bad;
      const_rh_fac(kkk) = keep_RH_const(xt(kkk),xt(kkk) + deltaT);
      xt(kkk) = xt(kkk) + deltaT;
      
      %% then do linear pert from CPRTOP to 500
      minp = 500;
      kkk = find(xp <= 500 & xp >= pxx.cprtop(ii));
      slope = (2*tBTD_DCC_bad) / abs(pxx.cprtop(ii) - 500);
      %[slope 500 pxx.cprtop(ii)]    
      deltaT = slope*(xp(kkk)-(pxx.cprtop(ii))) + -2*tBTD_DCC_bad;
      const_rh_fac(kkk) = keep_RH_const(xt(kkk),xt(kkk) + deltaT);      
      xt(kkk) = xt(kkk) + deltaT;
    
      %figure(5); semilogy(pxx.ptemp(1:nlevs,ii),pxx.plevs(1:nlevs,ii),'b',xt,pxx.plevs(1:nlevs,ii),'r');
      %set(gca,'ydir','reverse'); axis([180 300 10 1000]); grid
      %disp('redid ptemp');
      %pause(0.1)
    
      pyy.ptemp(1:nlevs,ii) = xt;
      pyy.gas_1(1:nlevs,ii) = xw .* const_rh_fac;      
      
      newcldp1(ii) = pxx.cprtop(ii);
      if pyy.ctype(ii) == 201
        newclds1(ii) = ice_size(tcld1(ii)-2*tBTD_DCC_bad,0);
      else
        newclds1(ii) = pxx.cpsize(ii);	
      end
    else
      aha = xp(kkk(aha));
      aha = max(aha);
      newcldp1(ii) = aha;
      if pyy.ctype(ii) == 201      
        newclds1(ii) = ice_size(tcld1(ii),0);
      else
        newclds1(ii) = pxx.cpsize(ii);	
      end
    end

    xp = pxx.plevs(1:nlevs,ii);
    yt = pyy.ptemp(1:nlevs,ii);
    if pxx.cfrac(ii) > 0 & pxx.cngwat(ii) > 0
      ycld1(ii) = interp1(log(xp),yt,log(pxx.cprtop(ii)));
    end
  end

  pyy.cpsize(oops1) = newclds1(oops1);
  pyy.cprtop(oops1) = newcldp1(oops1);
  pyy.cprbot(oops1) = newcldp1(oops1) + 50;

  figure(1); clf; plot(pxx.cprtop,pxx.cprtop - newcldp1,'o'); grid
  figure(1); clf; plot(pxx.cprtop,newcldp1,'o',newcldp1,newcldp1); grid; xlabel('old cprtop'); ylabel('new cprtop')
  figure(2); clf; plot(pxx.cprtop-pyy.cprtop,tcld1-ycld1,'.'); grid; xlabel('old-newcprtop'); ylabel('old-new tcld1')
  figure(3); clf; plot(pxx.cprtop,pxx.cpsize,'o'); grid  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(oops2) > 0
  tcld1    = ones(size(pxx.stemp))*NaN;
  ycld1    = ones(size(pxx.stemp))*NaN;
  tcld2    = ones(size(pxx.stemp))*NaN;
  newcldp2 = ones(size(pxx.stemp))*NaN;
  newclds2 = ones(size(pxx.stemp))*NaN;  
  for iii = 1 : length(oops2)
    ii = oops2(iii);
    nlevs = pxx.nlevs(ii);
    xp = pxx.plevs(1:nlevs,ii);
    xt = pxx.ptemp(1:nlevs,ii);
    xw = pxx.gas_1(1:nlevs,ii);
    const_rh_fac = ones(size(xw));
    yes = find(xt > 0 & xp > 0);
    if pxx.cfrac(ii) > 0 & pxx.cngwat(ii) > 0 & pxx.cprtop(ii) > 0
      %fprintf(1,'%5i %5i %8.6f \n',iii,ii,pxx.cprtop(ii));
      %xp
      %xt
      tcld1(ii) = interp1(log(xp(yes)),xt(yes),log(pxx.cprtop(ii)));
    end
    if pxx.cfrac2(ii) > 0 & pxx.cngwat2(ii) > 0 & pxx.cprtop2(ii) > 0
      tcld2(ii) = interp1(log(xp(yes)),xt(yes),log(pxx.cprtop2(ii)));
    end

    minp = max(50,pxx.cprtop(ii)-50);
    minp = 50;
    kkk = find(xp <= pxx.cprtop(ii) & xp > minp);
    ttt = xt(kkk);
    aha = find(ttt <= (tobs899(ii) - tBTD_DCC_bad));

    if length(aha) == 0
      fprintf(1,' OHOH oops2 \n');
      figure(5); semilogy(pxx.ptemp(1:nlevs,ii),pxx.plevs(1:nlevs,ii),'b');
      set(gca,'ydir','reverse'); axis([180 300 10 1000]); grid
      ax = axis; line([ax(1) ax(2)],        [pxx.cprtop(ii) pxx.cprtop(ii)],'color','k');
      ax = axis; line([tcld1(ii) tcld1(ii)],[ax(3) ax(4)],'color','k');
      if pyy.ctype(ii) == 201      
        newclds2(ii) = ice_size(tcld1(ii),0);
      else
        newclds2(ii) = pxx.cpsize(ii);	
      end            
    else
      aha = xp(kkk(aha));
      aha = max(aha);
      newcldp2(ii) = aha;
      if pyy.ctype(ii) == 201      
        newclds2(ii) = ice_size(tcld1(ii),0);
      else
        newclds2(ii) = pxx.cpsize(ii);
      end      
    end

    xp = pxx.plevs(1:nlevs,ii);
    yt = pyy.ptemp(1:nlevs,ii);
    yes = find(yt > 0 & xp > 0);    
    if pxx.cfrac(ii) > 0 & pxx.cngwat(ii) > 0
      ycld1(ii) = interp1(log(xp(yes)),yt(yes),log(pxx.cprtop(ii)));
    end
  end

  pyy.cpsize(oops2) = newclds2(oops2);
  pyy.cprtop(oops2) = newcldp2(oops2);
  pyy.cprbot(oops2) = newcldp2(oops2) + 50;

  figure(3); clf; plot(pxx.cprtop,pxx.cprtop - newcldp2,'o'); grid
  figure(3); clf; plot(pxx.cprtop,newcldp2,'o',newcldp2,newcldp2); grid; xlabel('old cprtop'); ylabel('new cprtop')
  figure(4); clf; plot(pxx.cprtop-pyy.cprtop,tcld1-ycld1,'.'); grid; xlabel('old-newcprtop'); ylabel('old-new tcld1')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(oops1) > 0
  [hjunk,pjunk] = subset_rtp(hxx,pyy,[],[],oops1);
  rtpwrite(tempx.tempfile2,hjunk,[],pjunk,[]);
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)

  [hjunk,hajunk,pjunk,pajunk] = rtpread(tempx.tempfile3);
  pyy.rcalc(:,oops1) = pjunk.rcalc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

if length(oops2) > 0
  [hjunk,pjunk] = subset_rtp(hxx,pyy,[],[],oops2);
  rtpwrite(tempx.tempfile2,hjunk,[],pjunk,[]);
  sartaer   = ['!' tempx.sarta   ' fin=' tempx.tempfile2 ' fout=' tempx.tempfile3 ' >& ' tempx.tempfile4]; eval(sartaer)

  [hjunk,hajunk,pjunk,pajunk] = rtpread(tempx.tempfile3);
  pyy.rcalc(:,oops2) = pjunk.rcalc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

