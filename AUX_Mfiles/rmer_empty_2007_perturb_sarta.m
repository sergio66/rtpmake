addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

iMakeDir = +1;  %% do    make dirs
iMakeDir = -1;  %% do not make dirs

iMvFiles = +1;  %% do    move file
iMvFiles = -1;  %% do not move files

iRemove = +1; %% do     remove
iRemove = -1; %% do not remove

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pertSTR = '10010';   %% water cloud, amt
pertSTR = '20010';   %% ice cloud,   amt
pertSTR = '11000';   %% water cloud, perturb cprtop by 10%
pertSTR = '21000';   %% water cloud, perturb cprtop by 10%
pertSTR = '12000';   %% water cloud, perturb cprtop by 20%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yy0 = 2007; mm0 = 01; dd0 = 06;

%% first make the subdirs
if iMakeDir > 0
  disp('Making subdirs ...')
  while yy0 == 2007
    dir0 = ['/asl/data/rtprod_airs/2007/' num2str(mm0,'%02d') '/' num2str(dd0,'%02d') '/'];
    cder = ['cd ' dir0];
    eval(cder);
    pertdir = ['PERT_' pertSTR];
    ee = exist(pertdir,'dir');
    if ee == 0
      mker = ['!mkdir ' pertdir];
      disp('making dir')
      eval(mker)
    end
    [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,16);
    yy0 = yy1; dd0 = dd1; mm0 = mm1;
  end
end

%% then move the files
if iMvFiles > 0
  disp('Moving files to subdirs ...')
  while yy0 == 2007
    dir0 = ['/asl/data/rtprod_airs/2007/' num2str(mm0,'%02d') '/' num2str(dd0,'%02d') '/'];
    cder = ['cd ' dir0];
    eval(cder);
    pertdir = ['PERT_' pertSTR];
    dirx = dir(['cloudy_airs_l1b_ecm_sarta_baum_ice_pert_' pertSTR '*.rtp']);
    if length(dirx) > 0
      mver = ...
        ['!/bin/mv cloudy_airs_l1b_ecm_sarta_baum_ice_pert_' pertSTR '*.rtp ' pertdir '/.'];
      fprintf(1,'%s \n',mver)
      eval(mver);
    end
    [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,16);
    yy0 = yy1; dd0 = dd1; mm0 = mm1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yy0 = 2007; mm0 = 01; dd0 = 06;

iDays = 0;
iOverallBad = 0;
iOverallGood = 0;
while yy0 == 2007
  dir0 = ['/asl/data/rtprod_airs/2007/' num2str(mm0,'%02d') '/' num2str(dd0,'%02d') '/'];
  cder = ['cd ' dir0];
  eval(cder);
  dirx = dir(['cloudy_airs_l1b_ecm_sarta_baum_ice_pert_' pertSTR '*.rtp']);
  iCnt = 0;
  iaBad = [];
  iDays = iDays + 1;
  for iii = 1 : length(dirx)
    if dirx(iii).bytes == 0
      iCnt = iCnt + 1;
      iaBad(iCnt) = iii;
    end
  end
  totalnum(iDays) = length(dirx);
  badnum(iDays)   = length(iaBad);
  iOverallGood = iOverallGood + length(dirx);
  if length(iaBad) > 0
    iOverallBad  = iOverallBad  + length(iaBad);
    iOverallGood = iOverallGood - length(iaBad);
    fprintf(1,'found %4i empty files for 2007/%2i/%2i for total bad %4i \n',length(iaBad),mm0,dd0,iOverallBad);
    if iRemove > 0
      for ii = 1 : length(iaBad)
        fname = dirx(iaBad(ii)).name;
        rmer = ['!/bin/rm ' fname];  
        fprintf(1,'   rming %s \n',fname);
        eval(rmer);
      end
    end
  end
  [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,16);
  yy0 = yy1; dd0 = dd1; mm0 = mm1;
end

badnum
totalnum
fprintf(1,'found %4i good files out of %4i; expecting 23*240 = %4i \n',iOverallGood,sum(totalnum),23*240);

cd /home/sergio/MATLABCODE/RTPMAKE/