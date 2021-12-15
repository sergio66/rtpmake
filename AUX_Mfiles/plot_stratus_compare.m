figure(1); clf
  simplemap(all_lat,all_lon,all_watertop)

for ii = 1 : length(all_lon);
  lala = all_waterprof(:,ii); 
  if sum(isfinite(lala)) > 1
    maxindex(ii) = find(lala == max(lala),1);
  else
    maxindex(ii) = NaN;
  end
end

addpath /strowdata1/shared/sergio/MATLABCODE/SHOWSTATS/
plevs = p.plevs(:,1);

figure(1); clf; plot(maxindex,all_watertop,'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(maxindex,all_watertop,1:1:37,700:25:1000,1); 
pcolor(plevs(x),y,log10(n))
xlabel('CWC Max Index (mb)'); ylabel('Cprtop Water (mb)'); colorbar
hl = line([700 1000],[700 1000],'linewidth',2); set(hl,'color','k')
title('ALL FOVS')

stratus = find(all_lon >= -110 & all_lon <= -75 & all_lat >= -20 & all_lat <= 0);
figure(2); clf; plot(maxindex(stratus),all_watertop(stratus),'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(maxindex(stratus),all_watertop(stratus),1:1:37,700:25:1000,1); 
pcolor(plevs(x),y,log10(n))
xlabel('CWC Max Index (mb)'); ylabel('Cprtop Water (mb)'); colorbar
hl = line([700 1000],[700 1000],'linewidth',2); set(hl,'color','k')
title('S. AMERICA STRATUS')

figure(3); clf; plot(maxindex,all_watertop,'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(nansum(all_waterprof),all_wateramt,(0:1:50)*1e-4,0:10:1000,1); 
pcolor(x,y,log10(n))
xlabel('CWC SUM'); ylabel('Cngwat (g/m2)'); colorbar
hl = line([0 0.005],[0 2000],'linewidth',2); set(hl,'color','k')
title('ALL FOVS')

figure(4); clf; plot(maxindex(stratus),all_watertop(stratus),'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(nansum(all_waterprof(:,stratus)),all_wateramt(stratus),(0:1:50)*1e-4,0:10:1000,1); 
pcolor(x,y,log10(n))
xlabel('CWC SUM'); ylabel('Cngwat (g/m2)'); colorbar
hl = line([0 0.005],[0 2000],'linewidth',2); set(hl,'color','k')
title('S.AMERICA FOVS')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); clf
  simplemap(xall_lat,xall_lon,xall_watertop)

for ii = 1 : length(xall_lon);
  lala = xall_waterprof(:,ii); 
  if sum(isfinite(lala)) > 1
    xmaxindex(ii) = find(lala == max(lala),1);
  else
    xmaxindex(ii) = NaN;
  end
end

addpath /strowdata1/shared/sergio/MATLABCODE/SHOWSTATS/
plevs = p.plevs(:,1);

figure(5); clf; plot(xmaxindex,xall_watertop,'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(xmaxindex,xall_watertop,1:1:37,700:25:1000,1); 
pcolor(plevs(x),y,log10(n))
xlabel('XCWC Max Index (mb)'); ylabel('Cprtop Water (mb)'); colorbar
hl = line([700 1000],[700 1000],'linewidth',2); set(hl,'color','k')
title('XALL FOVS')

stratus = find(xall_lon >= -110 & xall_lon <= -75 & xall_lat >= -20 & xall_lat <= 0);
figure(6); clf; plot(xmaxindex(stratus),xall_watertop(stratus),'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(xmaxindex(stratus),xall_watertop(stratus),1:1:37,700:25:1000,1); 
pcolor(plevs(x),y,log10(n))
xlabel('XCWC Max Index (mb)'); ylabel('Cprtop Water (mb)'); colorbar
hl = line([700 1000],[700 1000],'linewidth',2); set(hl,'color','k')
title('XS. AMERICA STRATUS')

figure(7); clf; plot(xmaxindex,xall_watertop,'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(nansum(xall_waterprof),xall_wateramt,(0:1:50)*1e-4,0:10:1000,1); 
pcolor(x,y,log10(n))
xlabel('XCWC SUM'); ylabel('Cnwat (g/m2)'); colorbar
hl = line([0 0.005],[0 2000],'linewidth',2); set(hl,'color','k')
title('XALL FOVS')

figure(8); clf; plot(xmaxindex(stratus),xall_watertop(stratus),'o'); 
%  [n x y]=myhist2d(fx,fy,limitsx, limitsy, LinearOrLog)
[n x y]=myhist2d(nansum(xall_waterprof(:,stratus)),xall_wateramt(stratus),(0:1:50)*1e-4,0:10:1000,1); 
pcolor(x,y,log10(n))
xlabel('XCWC SUM'); ylabel('Cnwat (g/m2)'); colorbar
hl = line([0 0.005],[0 2000],'linewidth',2); set(hl,'color','k')
title('XS.AMERICA FOVS')

figure(9)
dn = 0:0.01:2; nnCPRTOP = nanhist(xall_watertop./all_watertop,dn); 
               nnCNGWAT = nanhist(xall_wateramt./all_wateramt,dn); 
plot(dn,nnCPRTOP,'b',dn,nnCNGWAT,'r')
  title('(b) cprtop NEW/OLD  (r) cngwat NEW/OLD')