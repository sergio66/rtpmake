filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/MAGIC/magicsonde_matchup_airsfile_list.txt';
filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/TWP_quicksonde/twp_stuff.txt';
filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/TWP_quicksonde/twp_stuff_NEW.txt';
filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/TWP_quicksonde/twp_stuff_FIXED.txt';

filelist = load(filelist);
jjF = 0;
jjR = 0;
for ii = 1 : length(filelist)
  dada = filelist(ii,:);
  yy = dada(1);
  mm = dada(2);
  dd = dada(3);
  gg = dada(4);
  filename = ['/asl/data/rtprod_airs/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
  filename = [filename 'cloudy_airs_l1b_ecm_sarta_baum_ice.'];
  %filename = [filename 'cloudy_airs_l1b_era_sarta_baum_ice.'];
  filename = [filename num2str(yy,'%04d') '.' num2str(mm,'%02d') '.' num2str(dd,'%02d') '.' num2str(gg,'%03d') '.rtp'];
  if exist(filename)
    jjF = jjF + 1;
    thedir = dir(filename);
    %fprintf(1,'%3i %s %8i \n',ii,filename,thedir.bytes);
    if thedir.bytes == 0
      rmer = ['!/bin/rm ' filename];
      eval(rmer)
      fprintf(1,'%3i empty %s \n',ii,filename);
      jjR = jjR + 1;      
    end
  else
    fprintf(1,'%3i missing %s \n',ii,filename);
  end
end

fprintf(1,'there were %4i files in list, found %4i out of which %4i were empty \n',length(filelist),jjF,jjR)