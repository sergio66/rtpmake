%% run with
%%   all slurm output to one file, Paul does NOT recommend this
%%     sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%%   all slurm output to individual files, messy, but Paul recommends this
%%     sbatch --array=N1-N2 sergio_matlab_jobB.sbatch
%% N1 = 1, N2 = number of files to be processed

%% see breno_vs_sergio_cloud.m
%% see look_for_files_make_eracloudrtp_nadir_sarta_filelist.m

%% this goes after the two nadir tracks, after which you need to post-process
%% one day worth of data to get the nadir cosine sampling
%% see eg
%%     /home/sergio/MATLABCODE/ROSES_2013_GRANT/GlobalAvg/cluster_driver_makeavg.m
%%     /home/sergio/MATLABCODE/ROSES_2013_GRANT/GlobalAvg/final_pass.m
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% specify text file which has YY MM DD GG lst that needs to be processed
%{
yS = 2002;
yE = 2014; yE = 2015;
fid = fopen('eracloudrtp_nadir.txt','w');
for yy = yS : yE
  mmS = 1; mmE = 12;
  if yy == yS
    mmS = 09;
  elseif yy == yE
    mmE = 08;
  end
  if mod(yy,4) == 0
    mos = [31 29 31 30 31 30 31 31 30 31 30 31];  %% leap year
  else
    mos = [31 28 31 30 31 30 31 31 30 31 30 31];  %% normal year
  end
  for mm = mmS : mmE
    for dd = 1 : mos(mm)
      fprintf(fid,'%4i %2i %2i \n',yy,mm,dd);
    end
  end
end
fclose(fid);
%}

filelist = '/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA_CLOUD_NADIR/eracloudrtp_nadir.txt'; %% 20020901 to 20140831

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
% JOB=1585
% JOB=3798
% JOB=1
% JOB = 4470

fprintf(1,'JOB = %5i \n',JOB)

thefilelist = load(filelist);
thefilelist = thefilelist(JOB,:);

yymmdd0  = thefilelist(1:3); %% YY MM DD
fprintf(1,'JOB = %5i yy/mm/dd = %4i/%2i/%2i \n',JOB,yymmdd0(1),yymmdd0(2),yymmdd0(3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% creates an rtp file for ONE HOUR of data for that day
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlab2012/airs/readers
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

%addpath /home/strow/cress/Work/Rtp
%addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
%addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm

theinds = (1 : 2378)';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = 9999;

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

if codeX == 0
  icestr = '_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end

icestr = ['sergio_nadir_cloudy_airs_l1b_era' icestr '.'];              %% this is what was originally done
%% icestr = ['test_new_fillera_nadir_cloudy_airs_l1b_era' icestr '.'];    %% LLS does not know if he changed fill_era ??????

%%%%%%%%%%%%%%%%%%%%%%%%%

for hourloop = 0 : 23
  clear h p
  if hourloop < 23
    iaGlist = (1:10) + (hourloop)*10;  %% 1..10, 11..20, 21..30 etc 
  elseif hourloop == 23
    iaGlist = (1:11) + (hourloop)*10;  %% 231..241
  end

  year  = yymmdd0(1);
  month = yymmdd0(2);
  day   = yymmdd0(3);

  ystr = num2str(yymmdd0(1));
  mstr = num2str(yymmdd0(2),'%02d');
  dstr = num2str(yymmdd0(3),'%02d');
  hstr = num2str(hourloop,'%02d');

  if mod(year,4) == 0
    mos = [31 29 31 30 31 30 31 31 30 31 30 31];  %% leap year
  else
    mos = [31 28 31 30 31 30 31 31 30 31 30 31];  %% normal year
  end
  days_so_far = 0;
  if month > 1
    days_so_far = sum(mos(1:month-1));
  end
  days_so_far = days_so_far + day;
  filenameX = ['/asl/data/airs/AIRIBRAD/' ystr '/' num2str(days_so_far,'%03d') '/'];
  dir0 = filenameX;
  filenameX = [filenameX 'AIRS.' ystr '.' mstr '.' dstr '.*' '.L1B.AIRS_Rad.v5*.hdf'];
  eeX = dir(filenameX);
  
  fOUTDIR = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fnameOUT= [fOUTDIR icestr ystr '.' mstr '.' dstr '.' hstr '.rtp'];

  eeP = exist(fnameOUT);

  if eeP == 0 & length(eeX) > 0
    if ~ exist(fOUTDIR)
      mker = ['!mkdir -p ' fOUTDIR];
      eval(mker);
    end
    
    fprintf(1,' making %s \n',fnameOUT);
    toucher = ['!touch ' fnameOUT];
    eval(toucher)

    clear pTotal hTotal ha pa

    iFound = 0;
    for ixx = 1 : length(iaGlist)
      ix = iaGlist(ixx);
      yymmddgg = [yymmdd0 ix];

      gran  = yymmddgg(4);
      gstr  = num2str(gran,'%03d');

      filename = ['/asl/data/airs/AIRIBRAD/' ystr '/' num2str(days_so_far,'%03d') '/'];
      dir0 = filename;
      filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
      filename = [filename '.L1B.AIRS_Rad.v5*.hdf'];

      thedir = dir(filename);
      if length(thedir) == 1
        fname = [dir0 thedir.name];
	    
        %[meantime, f, prof] = readl1b_all(fname);  %% has the old rtime 1993
        [meantime, f, prof] = xreadl1b_all(fname);  %% has the new rtime 1958
	if length(f) > 0
          f0 = f;
	end	
	if length(meantime) == 0 & length(f) == 0 & length(prof) == 0
	  fprintf(1,'looks like NO data in L1B for %s \n',fname);
	else
          iFound = iFound + 1;	
          prof.upwell = ones(size(prof.solazi));
          prof.findex = ones(size(prof.solazi)) * gran;

          [a,b,c] = sdload(fname);
	  prof.DustScore = a.dust_score(:)';

          head.pfields = 4;  %% ir obs
          head.nchan = 2378;
 	  head.ichan = (1:2378)';
	
          nadir = find(prof.xtrack == 45 | prof.xtrack == 46);
          %% nadir = find((prof.xtrack >= 44 & prof.xtrack <= 47) & mod(prof.atrack,2) == 1);
	  prof = xyz_subset_rtp(prof,nadir);
	
          if iFound == 1
            hTotal.pfields = 4;  %% ir obs	
	    hTotal.nchan = 2378;
 	    hTotal.ichan = (1:2378)';
	    pTotal = prof;
         else
            [hTotal,pTotal] = cat_rtp(hTotal,pTotal,head,prof);
	  end
	end
      else
        fprintf(1,'%s DNE \n',filename);
      end       %% does AIRS L1B file exist
    end         %% loop over the 10 (or 11) files

    if (~exist('hTotal') & ~exist('pTotal'))
      fprintf(1,'oops does not look like AIRS L1B data was found to make %s \n',fnameOUT);
      rmer = ['!/bin/rm ' fnameOUT];
      eval(rmer);
    else
      %% the last call to xreadl1b could have been for a bad file, in which case f is [] and then
      %% h.vchan = f(h.ichan) few lines below will fail
      if length(f) == 0
        f = f0;
      end
      h = hTotal;
      p = pTotal;

      scatter(p.rlon,p.rlat,20,p.DustScore); colorbar
      p.pobs = zeros(size(p.solazi));

      pa = {{'profiles','rtime','seconds since 1958'}};
      ha = {{'header','hdf file',filename}};

      h.pfields=5; % (1=prof + 4=IRobs);
      h.nchan = length(theinds);
      h.ichan = theinds;;
      h.vchan = f(h.ichan);;

      %%% this is NEW
      p.landfrac_fromL1B = p.landfrac;
      p.salti_fromL1B    = p.salti;
      [salti, landfrac]  = usgs_deg10_dem(p.rlat, p.rlon);
      p.landfrac = landfrac;
      p.salti    = salti;

      %clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
      %cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
      %           'CC','CIWC','CLWC'};
      %[h,ha,p,pa] = rtpadd_era_data(h,ha,p,pa,cldfields); %%% add on era
      [p,h] = fill_era(p,h);
      if ~isfield(p,'tcc') & isfield(p,'cfrac')
        disp('set tcc = cfrac')
        p.tcc = p.cfrac;
      end
      
      p0 = p;

      %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  broken
      %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);       works
      addpath /asl/packages/time
      addpath /home/sbuczko1/git/rtp_prod2/emis/
      addpath /home/sbuczko1/git/rtp_prod2/util
      
      % [p,pa] = rtp_add_emis(p,pa);
      
      px = p;
      px.rlat = single(px.rlat);
      px.rlon = single(px.rlon);
      px.rtime = single(px.rtime);      
      [p,pa] = rtp_add_emis(px,pa);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      badcc0 = find(p.cc < 0); p.cc(badcc0) = 0;
      badcc1 = find(p.cc > 1); p.cc(badcc1) = 1;     
      badciwc = find(p.ciwc < 0); p.ciwc(badciwc) = 0;
      badclwc = find(p.clwc < 0); p.clwc(badclwc) = 0;
      junk = [length(badcc0) length(badcc1) length(badciwc) length(badclwc)];
      fprintf(1,'found %4i/%4i cc <0/1> %4i/%4i negative ciwc/clwc \n',junk);
      p.cfrac_ERA_ECM_seed = p.cfrac;      
      [p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

      fnamex = fnameOUT;
      [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
      rtpwrite(fnamex,h,ha,p2x,pa)
    end    %% hTotal,pTotal exist
  else
    fprintf(1,' %s already exists \n',fnameOUT)
  end      %% outfile exists
end        %% hourloop
