ocean = find(xpdmat.landfrac == 0);
land  = find(xpdmat.landfrac > 0.90);
if length(ocean)/length(xpdmat.landfrac) < 0.01
  ocean = land;
elseif length(land)/length(xpdmat.landfrac) < 0.01
  land = ocean;
end

if max(xpdmat.zPBLH_Ri) > 100
  xscale = 1000;  %% km --> m
else
  xscale = 1;     %% already in m
end

figure(3); clf; scatter_coast(xpdmat.rlon,xpdmat.rlat,50,xpdmat.zPBLH_Ri/xscale-ypdmat.zPBLH_Ri);        title('PBLH Ri ECMWF-UMBC');  colormap(usa2); caxis([-1 +1]*1.500)
figure(4); clf; dz = (-2000 : 100 : +2000)/1000; plot(dz,histc(xpdmat.zPBLH_Ri/xscale-ypdmat.zPBLH_Ri,dz));     title('PBLH Ri ECMWF-UMBC');  grid;
  fprintf(1,'mean(NWP = UMBC) = %8.6f m, std(NWP - UMBC) = %8.6f \n',mean(xpdmat.zPBLH_Ri/xscale-ypdmat.zPBLH_Ri),std(xpdmat.zPBLH_Ri/xscale-ypdmat.zPBLH_Ri))

figure(5); clf; plot(-diff(xpdmat.zalts(:,6000))/1000,xpdmat.zalts(2:91,6000)/1000,'bo-',-diff(ypdmat.zalts(:,6000))/1000,ypdmat.zalts(2:101,6000)/1000,'rx-')
  axis([0 2 0 5]); ylabel('hgt [km]'); xlabel('laayer thickness or diff between levels [km]'); legend('NWP levels xpdmat','UMBC layers ypdmat');
figure(6); clf; plot(xpdmat.Ri(:,6000),xpdmat.zalts(:,6000)/1000,'bo-',ypdmat.Ri(:,6000),ypdmat.zalts(:,6000)/1000,'rx-')
  axis([-2 2 0 5]); ylabel('hgt [km]'); xlabel('Ri []'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');

figure(5); clf; plot(-diff(nanmean(xpdmat.zalts(:,ocean),2))/1000,nanmean(xpdmat.zalts(2:91,ocean),2)/1000,'bo-',-diff(nanmean(ypdmat.zalts(:,ocean),2))/1000,nanmean(ypdmat.zalts(2:101,ocean),2)/1000,'rx-')
  axis([0 2 0 5]); ylabel('hgt [km]'); xlabel('laayer thickness or diff between levels [km]'); legend('NWP levels xpdmat','UMBC layers ypdmat');
figure(6); clf; plot(nanmean(xpdmat.Ri(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(ypdmat.Ri(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([-2 2 0 5]); ylabel('hgt [km]'); xlabel('Ri []'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');

figure(7); clf; plot(nanmean(xpdmat.ptemp(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(ypdmat.ptemp(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([200 300 0 15]); ylabel('hgt [km]'); xlabel('T(z) [K]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');
  axis([290 300 0 2])
figure(8); clf; semilogx(nanmean(xpdmat.gas_1(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(ypdmat.gg(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([0.001 0.1 0 15]); ylabel('hgt [km]'); xlabel('WV MR [g/g]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');
  axis([0.01 0.02 0 2])

lev_speed = sqrt((xpdmat.u).^2 + (xpdmat.v).^2);
lay_speed = sqrt((ypdmat.u).^2 + (ypdmat.v).^2);
figure(9); clf; plot(nanmean(lev_speed(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(lay_speed(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([0 10 0 15]); ylabel('hgt [km]'); xlabel('windspeed [m/s]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');
  axis([0 10 0 2])

figure(1); clf;  colormap(jet); scatter_coast(xpdmat.rlon,xpdmat.rlat,50,xpdmat.zPBLH_Ri/xscale - ypdmat.salti);    title('PBLH Ri from ECMWF');  cx = caxis;
figure(2); clf;  colormap(jet); scatter_coast(xpdmat.rlon,xpdmat.rlat,50,ypdmat.zPBLH_Ri - ypdmat.salti);           title('PBLH Ri from UMBC');   caxis(cx);
%[min(ypdmat.salti) max(ypdmat.salti) min(xpdmat.zPBLH_Ri/xscale) max(xpdmat.zPBLH_Ri/xscale) min(ypdmat.zPBLH_Ri) max(ypdmat.zPBLH_Ri)]

figure(10); clf; colormap(jet); scatter_coast(xpdmat.rlon,xpdmat.rlat,50,ypdmat.zPBLH_gg - ypdmat.salti);           title('PBLH gg from UMBC');   caxis(cx);
figure(11); clf; colormap(jet); scatter_coast(xpdmat.rlon,xpdmat.rlat,50,ypdmat.zPBLH_rh - ypdmat.salti);           title('PBLH rh from UMBC');   caxis(cx);
figure(12); clf; colormap(jet); scatter_coast(xpdmat.rlon,xpdmat.rlat,50,ypdmat.zPBLH_Tvir - ypdmat.salti);         title('PBLH Tvir from UMBC'); caxis(cx);
figure(13); clf; colormap(jet); scatter_coast(xpdmat.rlon,xpdmat.rlat,50,ypdmat.zPBLH_Tpot - ypdmat.salti);         title('PBLH Tpot from UMBC'); caxis(cx);

figure(1); caxis([0 4]);
figure(2); caxis([0 4]);
figure(10); caxis([0 4]);
figure(11); caxis([0 4]);
figure(12); caxis([0 4]);
figure(13); caxis([0 4]);

%{
figure(7); clf; semilogy(nanmean(xpdmat.ptemp(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'bo-',nanmean(ypdmat.ptemp(:,ocean),2),nanmean(ypdmat.plays(:,ocean),2),'rx-')
  axis([250 300 500 1050]); ylabel('hgt [km]'); xlabel('T(z) [K]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat'); set(gca,'ydir','reverse')
figure(8); clf; semilogx(nanmean(xpdmat.gas_1(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'bo-',nanmean(ypdmat.gg(:,ocean),2),nanmean(ypdmat.plays(:,ocean),2),'rx-')
  axis([0.001 0.1 500 1050]); ylabel('hgt [km]'); xlabel('WV MR [g/g]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat'); set(gca,'ydir','reverse')
%}

