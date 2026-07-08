%% this is basically copied from /home/sergio/MATLABCODE/oem_pkg_run/CRIS_new_clear_scan_January2020/read_the_anom_T_WV_retrievals.m

hpc_usage(-1)

simulateYear = 2012;

iCnt = 0;
for oo = 1 : 999
  fin  = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(oo,'%03d') '.mat'];
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(oo,'%04d') '.mat'];
  if exist(fin)
    iCnt = iCnt + 1;
    mver = ['!mv ' fin ' '  fout];
    eval(mver);
  end
end
if iCnt > 0
  fprintf(1,'had to move %4i files from ABC to WXYZ ie 3 digit to 4 digit \n',iCnt)
end

iaBad = zeros(1,4608);
for oo = 1 : 4608
  fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(oo,'%04d') '.mat'];
  if exist(fin)
    iaBad(oo) = 1;
  end
end

baddy = find(iaBad == 0);
if length(baddy) > 0
  disp('still need to get the cluster to do these timesteps')
  %for bb = 1 : length(baddy)
  %  str = ['sbatch -exclude=cnode203,cnode204,cnode260,cnode267 --array=' num2str(baddy(bb)) ' sergio_matlab_jobB.sbatch'];
  %  fprintf(fid,'%s \n',str);
  %end

  excludelist = ' ';
  excludelist = 'cnode014';

  partition = 'high_mem';
  partition = 'cpu2021';

  fid = fopen('badanom.sc','w');
  fprintf(1,'found that %4i of 4608 tiles  timesteps did not finish : see badanom.sc \n',length(baddy))
  str = ['sbatch -p ' partition ' --exclude=' excludelist ' --array='];
  fprintf(fid,'%s',str);
  iX = nice_output(fid,baddy);   %% now put in continuous strips, see /home/sergio/MATLABCODE/nice_output.m
  fprintf(1,'length(badanom) = %4i Num Continuous Strips = %4i \n',length(baddy),iX)
  str = [' sergio_matlab_jobB.sbatch 7'];
  fprintf(fid,'%s \n',str);
  fclose(fid);

end

figure(4);
clf; colormap jet
pcolor(reshape(iaBad,72,64)'); colorbar
hold on; plot(sum(reshape(iaBad,72,64),1),1:64,'y+-','linewidth',2); hold off
caxis([0 1]); xlabel('LonBin (1-72)');   ylabel('LatBin (1-64)'); title('Colorbar 0/1 and yellow=sum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
