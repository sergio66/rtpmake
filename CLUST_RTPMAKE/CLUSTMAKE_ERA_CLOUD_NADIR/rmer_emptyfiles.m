filelist = '/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA_CLOUD_NADIR/eracloudrtp_nadir.txt'; %% 20020901 to 20140831

%% /asl/data/rtprod_airs/2015/01/07/test_new_fillera_nadir_cloudy_airs_l1b_era_sarta_baum_ice.2015.01.07.00.rtp

thefilelist0 = load(filelist);

iX = 0;
for JOB = 4370:4780 
  thefilelist = thefilelist0(JOB,:);
  yymmdd0  = thefilelist(1:3); %% YY MM DD
  dir0 = ['/asl/data/rtprod_airs/' num2str(yymmdd0(1)) '/' num2str(yymmdd0(2),'%02d') '/' num2str(yymmdd0(3),'%02d') '/'];  
  thedir = dir([dir0 '/test_new_fillera_nadir_cloudy_airs_l1b_era_sarta_baum_ice*.rtp']);
  thedir = dir([dir0 '/sergio_nadir_cloudy_airs_l1b_era_sarta_baum_ice*.rtp']);
  for ii = 1 : length(thedir)
    fsz = thedir(ii).bytes;
    fnm = thedir(ii).name;
    if fsz == 0
      iX = iX + 1;
      rmer = ['!/bin/rm ' dir0 fnm];
      eval(rmer);
    end
  end
end
iX
