for pp = 1 : length(thepoints)
  if mod(pp,1000) == 0
    fprintf(1,'+')
  elseif mod(pp,100) == 0
    fprintf(1,'.')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% use raw data
  disp('raw T data')
  for ll = 1 : 101  
    for ii = 1 : 240    
      junk = squeeze(gas_1(ii,:,:,ll,pp));
      junk = sort(junk(:));
      gas_1_mean(ii) = mean(junk);
      gas_1_Q90(ii)  = mean(junk(indQ90));
  
      junk = squeeze(gas_1(ii,11,6,ll,pp));
      gas_1_cntr(ii) = junk;
    end

    %%%%%

    gas_1_cntr = gas_1_cntr/(eps+nanmean(gas_1_cntr));
    gas_1_mean = gas_1_mean/(eps+nanmean(gas_1_mean));
    gas_1_Q90  = gas_1_Q90/(eps+nanmean(gas_1_Q90));

    [B, err, stats] = Math_tsfit_lin_robust(days,gas_1_cntr,4);
    gas_1_trend_cntr(pp,ll) = B(2);
    gas_1_trend_cntr_err(pp,ll) = err.se(2);
  
    [B, err, stats] = Math_tsfit_lin_robust(days,gas_1_mean,4);
    gas_1_trend_mean(pp,ll) = B(2);
    gas_1_trend_mean_err(pp,ll) = err.se(2);
  
    [B, err, stats] = Math_tsfit_lin_robust(days,gas_1_Q90,4);
    gas_1_trend_Q90(pp,ll) = B(2);
    gas_1_trend_Q90_err(pp,ll) = err.se(2);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  disp('Q90 clr based')
  for ll = 1 : 101  
    for ii = 1 : 240    
      junkBT0 = squeeze(bt1231clr(ii,:,:,pp));
      [junkBT,IBT] = sort(junkBT0(:));

      junk = squeeze(gas_1(ii,:,:,ll,pp));
      gas_1_Q90(ii)  = mean(junk(IBT(indQ90)));
    end

    %%%%%
    gas_1_Q90  = gas_1_Q90/(eps+nanmean(gas_1_Q90));

    [B, err, stats] = Math_tsfit_lin_robust(days,gas_1_Q90,4);
    gas_1_trend_Q90_clr(pp,ll) = B(2);
    gas_1_trend_Q90_clr_err(pp,ll) = err.se(2);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  disp('Q90 cld based')
  for ll = 1 : 101  
    junkCT0 = squeeze(bt1231cld(ii,:,:,pp));
    [junkCT,ICT] = sort(junkCT0(:));

    for ii = 1 : 240    
      junk = squeeze(gas_1(ii,:,:,ll,pp));
      gas_1_Q90(ii)  = mean(junk(ICT(indQ90)));
    end

    %%%%%
    gas_1_Q90  = gas_1_Q90/(eps+nanmean(gas_1_Q90));

    [B, err, stats] = Math_tsfit_lin_robust(days,gas_1_Q90,4);
    gas_1_trend_Q90_cld(pp,ll) = B(2);
    gas_1_trend_Q90_cld_err(pp,ll) = err.se(2);
  end

end     %% loop over grid points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plays = load('/home/sergio/MATLABCODE/airslevels.dat');
plays = [1100*1.001; plays];
plays = flipud(meanvaluebin(plays));

figure(1); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_mean,   72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : mean');
figure(2); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_Q90,    72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : Q90');
figure(3); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_cntr,   72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : center');
figure(4); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_Q90_clr,72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : Q90 BT1231 Clr');
figure(5); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_Q90_cld,72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.1); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : Q90 BT1231 Cld');

figure(1); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_mean - gas_1_trend_mean,   72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : mean - mean');
figure(2); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_mean - gas_1_trend_Q90,    72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([ -1 +1]*0.01); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : mean-Q90');
figure(3); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_mean - gas_1_trend_cntr,   72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : mean - center');
figure(4); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_mean - gas_1_trend_Q90_clr,72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([ -1 +1]*0.01); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : mean-Q90 BT1231 Clr');
figure(5); clf; pcolor(Ylat,plays,squeeze(nanmean(reshape(gas_1_trend_mean - gas_1_trend_Q90_cld,72,length(thepoints)/72,101),1))'); 
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); shading interp;  caxis([-1 +1]*0.01); colorbar; colormap(llsmap5); ylim([10 1000]); title('GAS_1 trend : mean - Q90 BT1231 Cld');

