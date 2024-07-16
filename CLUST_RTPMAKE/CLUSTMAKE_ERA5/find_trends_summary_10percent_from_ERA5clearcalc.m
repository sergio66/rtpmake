days = (1:240)*365/12;

for pp = 1 : 4608
  if mod(pp,1000) == 0
    fprintf(1,'+')
  elseif mod(pp,100) == 0
    fprintf(1,'.')
  end

  for ii = 1 : 240
    junk = squeeze(bt1231clr(ii,:,:,pp));
    junk = sort(junk(:));
    bt1231clr_mean(ii) = mean(junk);
    bt1231clr_Q90(ii)  = mean(junk(indQ90));

    junk = squeeze(bt1231clr(ii,11,6,pp));
    bt1231clr_cntr(ii) = junk;

    %%%%%%%%%%%%%%%%%%%%%%%%%

    junk = squeeze(bt1231cld(ii,:,:,pp));
    junk = sort(junk(:));
    bt1231cld_mean(ii) = mean(junk);
    bt1231cld_Q90(ii)  = mean(junk(indQ90));

    junk = squeeze(bt1231cld(ii,11,6,pp));
    bt1231cld_cntr(ii) = junk;

    %%%%%%%%%%%%%%%%%%%%%%%%%

    junk = squeeze(stemp(ii,:,:,pp));
    junk = sort(junk(:));
    stemp_mean(ii) = mean(junk);
    stemp_Q90(ii)  = mean(junk(indQ90));

    junk = squeeze(stemp(ii,11,6,pp));
    stemp_cntr(ii) = junk;

  end

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231clr_cntr,4);
  bt1231clr_trend_cntr(pp) = B(2);
  bt1231clr_trend_cntr_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231clr_mean,4);
  bt1231clr_trend_mean(pp) = B(2);
  bt1231clr_trend_mean_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231clr_Q90,4);
  bt1231clr_trend_Q90(pp) = B(2);
  bt1231clr_trend_Q90_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231cld_cntr,4);
  bt1231cld_trend_cntr(pp) = B(2);
  bt1231cld_trend_cntr_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231cld_mean,4);
  bt1231cld_trend_mean(pp) = B(2);
  bt1231cld_trend_mean_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,bt1231cld_Q90,4);
  bt1231cld_trend_Q90(pp) = B(2);
  bt1231cld_trend_Q90_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_mean,4);
  skt_trend_mean(pp) = B(2);
  skt_trend_mean_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_Q90,4);
  skt_trend_Q90(pp) = B(2);
  skt_trend_Q90_err(pp) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(days,stemp_cntr,4);
  skt_trend_cntr(pp) = B(2);
  skt_trend_cntr_err(pp) = err.se(2);

end


whos *trend*
comment1 = 'created using data from   clust_loop_make_monthly_tile_273points.m --> clust_find_hottest_10percent_from_ERA5clearcalcs.m';
comment2 = 'then                       driver_load_in_summary_10percent_from_ERA5clearcalc.m --> find_trends_summary_10percent_from_ERA5clearcalc.m';
saver = ['save trends_summary_10percent_from_ERA5clearcalc.mat *trend* comment* '];
%eval(saver)

do_XX_YY_from_X_Y
load llsmap5
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
