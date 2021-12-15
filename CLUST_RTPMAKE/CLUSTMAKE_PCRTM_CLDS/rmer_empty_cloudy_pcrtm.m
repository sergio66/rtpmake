addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

yy0 = 2007; mm0 = 01; dd0 = 06;

iDays = 0;
iOverallBad = 0;
iOverallGood = 0;
while yy0 == 2007
  dir0 = ['/asl/data/rtprod_airs/2007/' num2str(mm0,'%02d') '/' num2str(dd0,'%02d') '/'];
  cder = ['cd ' dir0];
  eval(cder);
  dirx = dir('cloudy_airs_l1b_ecm_centerfovs_pcrtm_only*.rtp');
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
  iOverallGood = iOverallGood + length(dirx);
  if length(iaBad) > 0
    iOverallBad  = iOverallBad  + length(iaBad);
    iOverallGood = iOverallGood - length(iaBad);
    fprintf(1,'found %4i empty files for 2007/%2i/%2i for total bad %4i \n',length(iaBad),mm0,dd0,iOverallBad);
    for ii = 1 : length(iaBad)
      fname = dirx(iaBad(ii)).name;
      rmer = ['!/bin/rm ' fname];  
      fprintf(1,'   rming %s \n',fname);
      eval(rmer);
    end
  end
  [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,16);
  yy0 = yy1; dd0 = dd1; mm0 = mm1;
end

totalnum
fprintf(1,'found %4i good files \n',iOverallGood);

cd /home/sergio/MATLABCODE/RTPMAKE/