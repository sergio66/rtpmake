ddd = 1;

JOBB = (1:numJOBS) + (JOB-1)*numJOBS;

for JOBind = 1 : length(JOBB)
  JOBX = JOBB(JOBind);
  fprintf(1,'JOBX = %3i \n',JOBX);
  eeeXY_cnt = 0;
  for  eeeX = 1 : 21   %% 7X ERA grid points
    for eeeY = 1 : 12  %% 4Y ERA grid points
      eeeXY_cnt =   eeeXY_cnt + 1;
      %% dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_16day/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      %% fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_day_' num2str(ddd,'%02d') '_individual_timestep_' num2str(JOB,'%03d') '.mat'];
      fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_individual_timestep_' num2str(JOB,'%03d') '.mat'];

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
fprintf(1,' for %3i months have total of %5i files made from expected (%2i months x 21X pts x 12Y pts) = %6i \n',length(JOBB),sum(existence(:)),length(JOBB),length(JOBB)*21*12);

disp(' ')
for JOBind = 1 : length(JOBB)
  junk = squeeze(existence(JOBind,:,:));
  ind_existence(JOBind) = sum(junk(:));
end
fprintf(1,'this is breakdown of files made by month %3i to %3i \n',JOBB(1),JOBB(end))
printarray(ind_existence)
figure(1); plot(ind_existence,'o-'); title('Expect 21x12=252 files per month')
pause(0.1)
