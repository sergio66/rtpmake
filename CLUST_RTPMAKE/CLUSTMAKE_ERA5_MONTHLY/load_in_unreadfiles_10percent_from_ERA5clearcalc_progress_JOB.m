ddd = 1;

JOBB = (1:numJOBS) + (JOB-1)*numJOBS;

iFailSave = -1;
iFailCnt = 0;

%% this is new, March 2025
iSaveProfile = +1;

for JOBind = 1 : length(JOBB)    
  JOBX = JOBB(JOBind);
  fprintf(1,'JOBind = %3i --> JOBX = %3i     , out of %3i \n',JOBind,JOBX,length(JOBB));
  eeeXY_cnt = 0;
  for  eeeX = 1 : 21   %% 7X ERA grid points
    for eeeY = 1 : 12  %% 4Y ERA grid points
      iFail = -1;
      eeeXY_cnt =   eeeXY_cnt + 1;
      %% dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_16day/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      %% fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_day_' num2str(ddd,'%02d') '_individual_timestep_' num2str(JOBX,'%03d') '.mat'];
      fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_individual_timestep_' num2str(JOBX,'%03d') '.mat'];

      if existence(JOBind,eeeX,eeeY) == 1 & loadedfiles(JOBind,eeeX,eeeY) == 0                
        loader = ['a = load(''' fout ''');'];
        try
          eval(loader);
        catch
          disp('oh oh bad file, trying splitmat')
          foutjunk = [fout(1:end-4)];
          junkN = splitmat(foutjunk);
          fprintf(1,'expecting 10 splits, got %2i splits \n',junkN)
          try
            eval(loader);       
          catch
            fprintf(1,'even after splitaprt %s is still messed!!!!!! \n',fout);
            iFail = +1;
            iFailSave = +1;
            iFailCnt = iFailCnt + 1;
            badfile(iFailCnt,:) = fout;
          end
        end

        if ~isfield(a,'hnew_op') | ~isfield(a,'hnew_op') ~isfield(a,'pnew_op') | ~isfield(a,'pnew_op')
          fprintf(1,'oh oh %s does not have hnew_ip/op or pnew_ip/op \n',fout)
          iFailSave = +1;
          iFailCnt = iFailCnt + 1;
          badfile(iFailCnt,:) = fout;
        end
      
        if iFail < 0 
          bt1231cld(JOBind,eeeX,eeeY,:) = rad2bt(1231,a.pnew_op.rcalc);
          bt1231clr(JOBind,eeeX,eeeY,:) = rad2bt(1231,a.pnew_op.sarta_rclearcalc);
          stemp(JOBind,eeeX,eeeY,:)     = a.pnew_op.stemp;
          landfrac(JOBind,eeeX,eeeY,:)  = a.pnew_op.landfrac;
          loadedfiles(JOBind,eeeX,eeeY) = 1;

          %% this is new, March 2025
          if iSaveProfile > 0
            ptemp(JOBind,eeeX,eeeY,:,:)     = a.pnew_op.ptemp;
            gas_1(JOBind,eeeX,eeeY,:,:)     = a.pnew_op.gas_1;
            %gas_3(JOBind,eeeX,eeeY,:,:)     = a.pnew_op.gas_3;
          end

        end

      end

    end
  end

  ind = 1000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(2); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 
  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOBind = ' num2str(JOBind) ' BT1231 clr/cld for Tile ' num2str(ind)])

  ind = 2000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(3); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 
  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOBind = ' num2str(JOBind) ' BT1231 clr/cld for Tile ' num2str(ind)])

  ind = 3000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(4); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 
  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOBind = ' num2str(JOBind) ' BT1231 clr/cld for Tile ' num2str(ind)])

  ind = 4000;
  moo = squeeze(bt1231clr(JOBind,:,:,ind)); moo = moo(:);  oo = find(moo > 0.1 & isfinite(moo)); moo = moo(oo); 
  zoo = squeeze(bt1231cld(JOBind,:,:,ind)); zoo = zoo(:);  oo = find(zoo > 0.1 & isfinite(zoo)); zoo = zoo(oo); 
  soo = squeeze(stemp(JOBind,:,:,ind)); soo = soo(:);      oo = find(soo > 0.1 & isfinite(soo)); soo = soo(oo); 
  figure(5); plot(1:length(moo),moo,'b',1:length(moo),zoo,'k',1:length(moo),soo,'r'); 

  hl = legend('CLR','CLD','ST','location','best','fontsize',10); title(['JOBind = ' num2str(JOBind) ' BT1231 clr/cld for Tile ' num2str(ind)])

  pause(0.1)
end

if iFailSave < 0 & iSaveProfile < 0
  filesave = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_' num2str(JOB,'%02d') '.mat'];
  saver = ['save ' filesave ' existence loadedfiles bt1231clr bt1231cld stemp landfrac JOBB'];
  if exist(filesave)
    fprintf(1,'%s already exists, not saving \n',filesave)
  else
    fprintf(1,'%s DNE, saving \n',filesave)
    eval(saver);
  end

elseif iFailSave < 0 & iSaveProfile > 0
  filesave = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_with_Tprofile_' num2str(JOB,'%02d') '.mat'];
  saver = ['save -v7.3 ' filesave ' existence loadedfiles bt1231clr bt1231cld landfrac stemp ptemp gas_1 gas_3 JOBB'];
  saver = ['save -v7.3 ' filesave ' existence loadedfiles bt1231clr bt1231cld landfrac stemp ptemp JOBB'];
  if exist(filesave)
    fprintf(1,'%s already exists, not saving \n',filesave)
  else
    fprintf(1,'%s DNE, saving \n',filesave)
    eval(saver);
  end

  filesave = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_with_WVprofile_' num2str(JOB,'%02d') '.mat'];
  saver = ['save -v7.3 ' filesave ' existence loadedfiles bt1231clr bt1231cld landfrac stemp ptemp gas_1 gas_3 JOBB'];
  saver = ['save -v7.3 ' filesave ' existence loadedfiles bt1231clr bt1231cld landfrac stemp gas_1 JOBB'];
  if exist(filesave)
    fprintf(1,'%s already exists, not saving \n',filesave)
  else
    fprintf(1,'%s DNE, saving \n',filesave)
    eval(saver);
  end

  %filesave = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_with_OZprofile_' num2str(JOB,'%02d') '.mat'];
  %saver = ['save -v7.3 ' filesave ' existence loadedfiles bt1231clr bt1231cld landfrac stemp ptemp gas_1 gas_3 JOBB'];
  %saver = ['save -v7.3 ' filesave ' existence loadedfiles bt1231clr bt1231cld landfrac stemp gas_3 JOBB'];
  %if exist(filesave)
  %  fprintf(1,'%s already exists, not saving \n',filesave)
  %else
  %  fprintf(1,'%s DNE, saving \n',filesave)
  %  eval(saver);
  %end

else
  disp('iFailSave > 0 so not saving')
  for ii = 1 : iFailCnt
    fprintf(1,'%3i failed %s \n',ii,badfile(ii,:))
  end
end
