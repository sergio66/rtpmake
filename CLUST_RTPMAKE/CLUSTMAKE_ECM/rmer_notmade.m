findstr = ['!find /asl/s1/sergio/rtp/rtp_airibrad_v5 -type f -empty']; eval(findstr)
findstr = ['!find /asl/s1/sergio/rtp/rtp_airibrad_v5 -type f -empty | wc -l']; eval(findstr)
findstr = ['!find /asl/s1/sergio/rtp/rtp_airibrad_v5 -type f -empty -delete']; %eval(findstr)
disp('ret'); pause

thedir = dir(['/asl/s1/sergio/rtp/rtp_airibrad_v5/20*/*/*/*.rtp']);
iCnt = 0;
for ii = 1 : length(thedir)
  xname = thedir(ii).name;
  xsize = thedir(ii).bytes;
  if xsize == 0
    iCnt = iCnt + 1;
    info = xname(36:50);
    ystr = info(1:4); mstr = info(6:7); dstr = info(9:10); gstr = info(12:14);  
    fullname = ['/asl/s1/sergio/rtp/rtp_airibrad_v5/' ystr '/' mstr '/' dstr '/cloudy_airs_l1b_era_sarta_baum_ice.' info 'rtp'];
  
    lser = ['!ls -lt ' fullname];
    eval(lser);
    rmer = ['!/bin/rm  ' fullname];
    eval(rmer);
  end
end
fprintf(1,'found %5i files of which %5i were size zero McLaren Honda Renault failures \n',length(thedir),iCnt);