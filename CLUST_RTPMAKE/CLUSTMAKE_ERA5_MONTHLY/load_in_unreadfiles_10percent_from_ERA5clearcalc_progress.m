ddd = 1;

for JOB = JOBB(1) : JOBB(end)  %% 20 years
  JOBind = JOB - JOBB(1) + 1;

  fprintf(1,'JOB = %3i \n',JOB);
  eeeXY_cnt = 0;
  for  eeeX = 1 : 21   %% 7X ERA grid points
    for eeeY = 1 : 12  %% 4Y ERA grid points
      eeeXY_cnt =   eeeXY_cnt + 1;
      %% dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_16day/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      %% fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_day_' num2str(ddd,'%02d') '_individual_timestep_' num2str(JOB,'%03d') '.mat'];
      if  run_sarta.cumsum == -1
        fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_individual_timestep_' num2str(JOB,'%03d') '.mat'];
      elseif run_sarta.cumsum == 9999
        fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_individual_timestep_' num2str(JOB,'%03d') '_9999.mat'];
      end

      if existence(JOBind,eeeX,eeeY) == 1 & loadedfiles(JOBind,eeeX,eeeY) == 0
        loader = ['a = load(''' fout ''');'];
        eval(loader);
        bt1231cld(JOBind,eeeX,eeeY,:) = rad2bt(1231,a.pnew_op.rcalc);
        bt1231clr(JOBind,eeeX,eeeY,:) = rad2bt(1231,a.pnew_op.sarta_rclearcalc);
        stemp(JOBind,eeeX,eeeY,:)     = a.pnew_op.stemp;
        landfrac(JOBind,eeeX,eeeY,:)  = a.pnew_op.landfrac;
        loadedfiles(JOBind,eeeX,eeeY) = 1;
      end

    end
  end

  ind = 1000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(2); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 
  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOB = ' num2str(JOB) ' BT1231 clr/cld for Tile ' num2str(ind)])

  ind = 2000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(3); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 
  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOB = ' num2str(JOB) ' BT1231 clr/cld for Tile ' num2str(ind)])

  ind = 3000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(4); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 
  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOB = ' num2str(JOB) ' BT1231 clr/cld for Tile ' num2str(ind)])

  ind = 4000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(5); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 
  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOB = ' num2str(JOB) ' BT1231 clr/cld for Tile ' num2str(ind)])

  if run_sarta.cumsum == 9999
    moo = squeeze(bt1231cld);
    dbt = 200:340;
    for jj = 1 : 64
      jjind = (1:72) + (jj-1)*72;
      junk = moo(:,:,jjind);
      junk = junk(:);
      histcld(jj,:) = histc(junk,dbt)/length(junk);
    ende
  figure(6); clf; pcolor(dbt,1:64,histcld); colormap(jet); shading interp; colorbar; title('BT1231 cld')
  caxis([0 0.08]); jett = jet(256); jett(1,:) = 1; colormap(jett); xlim([210 310])
  end

  pause(0.1)
end
