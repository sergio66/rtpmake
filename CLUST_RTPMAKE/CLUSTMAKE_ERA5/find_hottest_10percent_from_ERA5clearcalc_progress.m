if run_sarta.cumsum == -1
  JOBB = 1 : 240
elseif run_sarta.cumsum == 9999
  JOBB = 120
end

ddd = 1;
for JOB = JOBB(1) : JOBB(end)    %% 20 years
  JOBind = JOB - JOBB(1) + 1;

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

      ee = exist(fout);
      if eeeX == 1 & eeeY == 1
        fprintf(1,'JOB = %3i xscan = %2i of 21 (5 deg lon @ 0.25 deg res)     yscan = %2i of 12 (3 deg lat @0.25 deg res)   exist = %2i \n',JOB,eeeX,eeeY,ee)
      end

      if ee > 0
        %% rmer = ['!/bin/rm ' fout]; eval(rmer); %%%%%%%        
        existence(JOBind,eeeX,eeeY) = 1;
      else
        existence(JOBind,eeeX,eeeY) = 0;
      end
    end
  end
end
fprintf(1,' have grand total of %6i files done out of out of 240 months (20 yrs) x (21 X pts x 12 Y pts) = %6i \n',sum(existence(:)),12*20*21*12);

disp(' ')
for JOB = JOBB(1) : JOBB(end)
  JOBind = JOB - JOBB(1) + 1;
  junk = squeeze(existence(JOBind,:,:));
  ind_existence(JOBind) = sum(junk(:));
end
disp('this is breakdown of files made by month 1-240, we expect max (21 x 12 = 252)')
printarray(ind_existence)
figure(1); plot(ind_existence,'o-'); title('Expect 21x12=252 files per month')
pause(0.1)

bad = find(ind_existence == 0); 
if length(bad) > 0
  printarray(bad,'these are the months that so far have ZERO files made');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
baddy = find(ind_existence < 252); 
if length(baddy) > 0
  printarray(baddy,'these are the months that so far have less than 252 files made');

  %for bb = 1 : length(baddy)
  %  str = ['sbatch -exclude=cnode203,cnode204,cnode260,cnode267 --array=' num2str(baddy(bb)) ' sergio_matlab_jobB.sbatch'];
  %  fprintf(fid,'%s \n',str);
  %end

  disp('still need to get the cluster to do these : look at jobs_not_done.sc in current working dir')
  disp('                                          : XYZ is what YOU need to edit and set')
  disp('                                          : also look at notdone.txt for the list');
  fid = fopen('notdone.txt','w');
  fprintf(fid,'%6i \n',baddy);
  fclose(fid);

  excludelist = 'cnode014';
  excludelist = ' ';

  partition = 'high_mem';
  partition = 'cpu2021';

  fid = fopen('jobs_not_done.sc','w');
  fprintf(1,'found that %4i of %4i timesteps did not finish : see jobs_not_done.sc \n',length(baddy),240)
  str = ['sbatch -p ' partition ' --exclude=' excludelist ' --array='];
  fprintf(fid,'%s',str);
  iX = nice_output(fid,baddy);   %% now put in continuous strips, see /home/sergio/MATLABCODE/nice_output.m
  fprintf(1,'length(jobs_not_done) = %4i Num Continuous Strips = %4i \n',length(baddy),iX)
  str = [' sergio_matlab_jobB.sbatchX XYZ'];
  str = [' sergio_matlab_jobB.sbatch 5'];
  fprintf(fid,'%s \n',str);
  fclose(fid);

  %disp('last part of jobs_not_done.sc is "sergio_matlab_jobB.sbatchX XYZ" which is wrong : edit as needed!!!!')
  %disp('last part of jobs_not_done.sc is "sergio_matlab_jobB.sbatchX XYZ" which is wrong : edit as needed!!!!')
  %disp('last part of jobs_not_done.sc is "sergio_matlab_jobB.sbatchX XYZ" which is wrong : edit as needed!!!!')

  chmodder = ['!chmod +x jobs_not_done.sc'];
  eval(chmodder)
end


%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('directly doing ls -lt blah | wc -l')
if run_sarta.cumsum == -1
  junk = ['!ls -lt /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex*/*.mat | wc -l'];
elseif run_sarta.cumsum == 9999
  junk = ['!ls -lt /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex*/*_9999.mat | wc -l'];
end

eval(junk);

pause(0.1)
