days = (1:240)*365/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pp = 1 : length(thepoints)
  if mod(pp,1000) == 0
    fprintf(1,'+')
  elseif mod(pp,100) == 0
    fprintf(1,'.')
  end

  for ii = 1 : 240

    %%%%%%%%%%%%%%%%%%%%%%%%%
    junkST0 = squeeze(stemp(ii,:,:,pp));
    [junkST,IST] = sort(junkST0(:));
    stemp_mean(ii) = mean(junkST);
    stemp_Q90(ii)  = mean(junkST(indQ90));

    junk = squeeze(stemp(ii,11,6,pp));
    stemp_cntr(ii) = junk;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    junkBT0 = squeeze(bt1231clr(ii,:,:,pp));
    [junkBT,IBT] = sort(junkBT0(:));
    bt1231clr_mean(ii) = mean(junkBT);
    bt1231clr_Q90(ii)  = mean(junkBT(indQ90));

    junk = squeeze(bt1231clr(ii,11,6,pp));
    bt1231clr_cntr(ii) = junk;

    %% based on BT1231clr Q90
    stemp_Q90_clr(ii) = mean(junkST0(IBT(indQ90)));

    %%%%%%%%%%%%%%%%%%%%%%%%%
    junkCT0 = squeeze(bt1231cld(ii,:,:,pp));
    [junkCT,ICT] = sort(junkCT0(:));
    bt1231cld_mean(ii) = mean(junkCT);
    bt1231cld_Q90(ii)  = mean(junkCT(indQ90));

    junk = squeeze(bt1231cld(ii,11,6,pp));
    bt1231cld_cntr(ii) = junk;

    %% based on BT1231cld Q90
    stemp_Q90_cld(ii) = mean(junkST0(ICT(indQ90)));

    %%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    figure(1); plot(1:252,junkST0(IST),1:252,junkBT0(IBT),1:252,junkCT0(ICT)); legend('ST','BT clr','BT cld','location','best','fontsize',10); title('sort key : themselves')
    figure(2); plot(1:252,junkST0(IST),1:252,junkBT0(IST),1:252,junkCT0(IST)); legend('ST','BT clr','BT cld','location','best','fontsize',10); title('sort key : ST')
    figure(3); plot(1:252,junkST0(IBT),1:252,junkBT0(IBT),1:252,junkCT0(IBT)); legend('ST','BT clr','BT cld','location','best','fontsize',10); title('sort key : BT Clr')
    figure(4); plot(1:252,junkST0(ICT),1:252,junkBT0(ICT),1:252,junkCT0(ICT)); legend('ST','BT clr','BT cld','location','best','fontsize',10); title('sort key : BT Cld')
    %}
    
  end

  %%%%%

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231clr_cntr,4);
  bt1231clr_trend_cntr(pp) = B(2);
  bt1231clr_trend_cntr_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231clr_mean,4);
  bt1231clr_trend_mean(pp) = B(2);
  bt1231clr_trend_mean_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231clr_Q90,4);
  bt1231clr_trend_Q90(pp) = B(2);
  bt1231clr_trend_Q90_err(pp) = err.se(2);

  %%%%%

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231cld_cntr,4);
  bt1231cld_trend_cntr(pp) = B(2);
  bt1231cld_trend_cntr_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231cld_mean,4);
  bt1231cld_trend_mean(pp) = B(2);
  bt1231cld_trend_mean_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231cld_Q90,4);
  bt1231cld_trend_Q90(pp) = B(2);
  bt1231cld_trend_Q90_err(pp) = err.se(2);

  %%%%%

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_mean,4);
  skt_trend_mean(pp) = B(2);
  skt_trend_mean_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_Q90,4);
  skt_trend_Q90(pp) = B(2);
  skt_trend_Q90_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_cntr,4);
  skt_trend_cntr(pp) = B(2);
  skt_trend_cntr_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_Q90_clr,4);
  skt_trend_Q90_clr(pp) = B(2);
  skt_trend_Q90_clr_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_Q90_cld,4);
  skt_trend_Q90_cld(pp) = B(2);
  skt_trend_Q90_cld_err(pp) = err.se(2);

  %%%%%

  %plot(1:240,[stemp_mean; stemp_Q90; stemp_cntr; stemp_Q90_clr; stemp_Q90_cld]); title('5 different ways of computing \newline STEMP for a tile');
  %[skt_trend_mean(pp) skt_trend_Q90(pp) skt_trend_cntr(pp) skt_trend_Q90_cldr(pp) skt_trend_Q90_cld(pp)]
  %legend('mean over tile','Q90 for the tile','title center','based on BT1231clr Q90','based on BT1231cld Q90','location','best','fontsize',10);

end

load llsmap5

figure(1); clf; pcolor(reshape(skt_trend_mean,   72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); title('SKT trend : mean SKT');
figure(2); clf; pcolor(reshape(skt_trend_Q90,    72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); title('SKT trend : Q90 SKT');
figure(3); clf; pcolor(reshape(skt_trend_cntr,   72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); title('SKT trend : center SKT');
figure(4); clf; pcolor(reshape(skt_trend_Q90_clr,72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); title('SKT trend : Q90 BT1231 Clr SKT');
figure(5); clf; pcolor(reshape(skt_trend_Q90_cld,72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); title('SKT trend : Q90 BT1231 Cld SKT');

figure(1); clf; pcolor(reshape(skt_trend_mean - skt_trend_mean,   72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); title('SKT trend : mean SKT - mean SKT');
figure(2); clf; pcolor(reshape(skt_trend_mean - skt_trend_Q90,    72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); title('SKT trend : mean SKT - Q90 SKT');
figure(3); clf; pcolor(reshape(skt_trend_mean - skt_trend_cntr,   72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); title('SKT trend : mean SKT - center SKT');
figure(4); clf; pcolor(reshape(skt_trend_mean - skt_trend_Q90_clr,72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); title('SKT trend : mean SKT - Q90 BT1231 Clr SKT');
figure(5); clf; pcolor(reshape(skt_trend_mean - skt_trend_Q90_cld,72,length(thepoints)/72)); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); title('SKT trend : mean SKT - Q90 BT1231 Cld SKT');

figure(6); clf; 
junk =        nanmean(reshape(skt_trend_mean,   72,length(thepoints)/72),1);
junk = [junk; nanmean(reshape(skt_trend_Q90,    72,length(thepoints)/72),1)];
junk = [junk; nanmean(reshape(skt_trend_cntr,   72,length(thepoints)/72),1)];
junk = [junk; nanmean(reshape(skt_trend_Q90_clr,72,length(thepoints)/72),1)];
junk = [junk; nanmean(reshape(skt_trend_Q90_cld,72,length(thepoints)/72),1)];
plot(Ylat,junk,'linewidth',2); xlabel('Latitude'); ylabel('dSKT/dt [K/yr]');
legend('SKT mean','SKT Q90','SKT tile center','SKT from Q90 BT1231 clr','SKT from Q90 BT1231 cld','location','best','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iTorWV > 0
  find_Tz_trends_tile_center_mean_quantile_etc
elseif iTorWV < 0
  find_WVz_trends_tile_center_mean_quantile_etc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


whos *trend*
comment0 = 'see Readme_tilecenter_vs_hottest10percent_vs_allpintsintile_for_trends : indices into the 4608 tiles are from /home/sergio/git/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/test_read_eraI.m';
comment1 = 'created using data from   clust_loop_make_monthly_tile_273points.m --> clust_find_hottest_10percent_from_ERA5clearcalcs.m';
comment2 = 'then                      driver_load_in_summary_10percent_from_ERA5clearcalc.m WITH PROFILE OPTION --> clust_find_hottest_10percent_from_ERA5clearcalcs --> driver_load_in_summary_10percent_from_ERA5clearcalc_profile.m';

if iTorWV > 0
  fsaveout = ['trends_summary_10percent_from_ERA5clearcalc_T_region' num2str(iSub,'%02d') '.mat'];
elseif iTorWV < 0
  fsaveout = ['trends_summary_10percent_from_ERA5clearcalc_WV_region' num2str(iSub,'%02d') '.mat'];
end

if ~exist(fsaveout)
  fprintf(1,'saving to %s \n',fsaveout)
  saver = ['save ' fsaveout ' *trend* comment* Ylat yindex iTorWV iSub thepoints'];
  eval(saver)
else
  fprintf(1,'outfile already exists, not saving %s \n',fsaveout)  
end

disp('now put together using driver_put_together_5_regions_ERA5_tile_trends.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
do_XX_YY_from_X_Y
aslmap(1,rlat65,rlon73,smoothn((reshape(skt_trend_mean',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 SKT trend mean')
aslmap(2,rlat65,rlon73,smoothn((reshape(skt_trend_Q90',72,64)')  ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 SKT trend Q90')
aslmap(3,rlat65,rlon73,smoothn((reshape(skt_trend_cntr',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 SKT trend cntr')
aslmap(4,rlat65,rlon73,smoothn((reshape(skt_trend_mean'-skt_trend_Q90',72,64)')  ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.25/10); colormap(llsmap5); title('ERA5 SKT trend : mean-Q90')
aslmap(5,rlat65,rlon73,smoothn((reshape(skt_trend_mean'-skt_trend_cntr',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.25/10); colormap(llsmap5); title('ERA5 SKT trend : mean-cntr')

figure(6); dst = (-1:0.05:+1)*0.025; plot(dst,histc(skt_trend_mean-skt_trend_Q90,dst),'b',dst,histc(skt_trend_mean-skt_trend_cntr,dst),'r'); plotaxis2;
title('SKT trend diff mean-X'); hl = legend('Q90','CNTR','location','best','fontsize',10);
x1 = skt_trend_mean-skt_trend_Q90;
x2 = skt_trend_mean-skt_trend_cntr;
printarray([mean(x1) std(x1) mean(x2) std(x2)],'mean/std for meanSKT-Q90SKT and  meanSKT-cntrSKT')
%}
