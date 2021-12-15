addpath /asl/matlab/science
addpath /asl/matlab/rtptoolsV201
addpath /asl/matlab/h4tools/
addpath /asl/matlab/rtptools/
addpath /asl/matlab/gribtools/
addpath /asl/matlab/airs/readers/
addpath /home/sergio/MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE
addpath /home/sergio/MATLABCODE/PLOTTER

iaChanList = [1:2378];
iaChanList = [359 445 903 1291 1557 1614];  %% for MATLABCODE/RATES_CLOUD

%%%%%%%%%%%%%%%%%%%%%%%%%
rCNGWATiceMult = 1.0;
rCPSIZEiceMult = 1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%

rCfracType = 0.5; %% use constant |x| for cfrac1,cfrac2, cfrac12
rCfracType = 0.99; %% use constant |x| for cfrac1,cfrac2, cfrac12
rCfracType = -1;  %% use wghted average      of sum(CC) in cfrac/cfrac2/cfrac12
                  %%                         almost rCfracType = 0, but uses CC (91xn)
rCfracType = -2;  %% use wghted average of sum(CC) (91xn), then SAME cfrac/cfrac2/cfrac12
rCfracType = 0;   %% default, use Scott's methods which uses TCC (1xn)

rCfracType = -2;  %% use wghted average of sum(CC) (91xn), then SAME cfrac/cfrac2/cfrac12
rCfracType = -1;  %% use wghted average      of sum(CC) in cfrac/cfrac2/cfrac12
                  %%                         almost rCfracType = 0, but uses CC (91xn)

%%%%%%%%%%%%%%%%%%%%%%%%%
rCNGWATiceMult = 5.0; rCPSIZEiceMult = 1.0;
rCfracType = 0;   %% default, use Scott's methods which uses TCC (1xn)
rCfracType = -1;  %% use wghted average      of ccfrac in cfrac/cfrac2/cfrac12
                  %%                         almost rCfracType = 0, but uses CC (91xn)
rCNGWATiceMult = 10.0; rCPSIZEiceMult = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE

iCnt = 0;
yymmddgg = [2011 03 10];
for gg = 1 : 240
  disp(' ')
  %clear p h hattr pattr prof yymmddgg
  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(gg,'%03d');

  fnamex = ['/strowdataN/s3/robinso2/PUYEHUE/CLOUDS_6chans/'];

  if (abs(rCNGWATiceMult-1) < eps) & (abs(rCPSIZEiceMult-1) < eps) 
    %% DEFAULT : no ice amount or ice size multiplier
    if rCfracType == 0
      fnamex = [fnamex 'TYPE0/eracld_type0_allfov_'];
    elseif rCfracType == -1
     fnamex = [fnamex 'TYPE1/eracld_type1_allfov_'];
    elseif rCfracType == -2
      fnamex = [fnamex 'TYPE2/eracld_type2_allfov_'];
    elseif rCfracType > 0  & rCfracType <= 1
      fnamex = [fnamex 'TYPEX/eracld_typeX_allfov_'];
    end
  elseif (abs(rCNGWATiceMult-1) > eps) & (abs(rCPSIZEiceMult-1) < eps) 
    %% have put ice amount multiplier, but no ice size multiplier
    if rCfracType == 0
      fnamex = [fnamex 'TYPEI_rCfracType=0_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/eracld_type0_allfov_'];
    elseif rCfracType == -1
      fnamex = [fnamex 'TYPEI_rCfracType=-1_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/eracld_type0_allfov_'];
    elseif rCfracType == -2
      fnamex = [fnamex 'TYPEI_rCfracType=-2_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/eracld_type0_allfov_'];
    end
  elseif (abs(rCNGWATiceMult-1) > eps) & (abs(rCPSIZEiceMult-1) > eps) 
    %% have put ice amount multiplier and ice size multiplier
    if rCfracType == 0
      fnamex = [fnamex 'TYPEIS_rCfracType=0_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/eracld_type0_allfov_'];
    elseif rCfracType == -1
      fnamex = [fnamex 'TYPEIS_rCfracType=-1_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/eracld_type0_allfov_'];
    elseif rCfracType == -2
      fnamex = [fnamex 'TYPEIS_rCfracType=-2_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/eracld_type0_allfov_'];
    else
      error('dude!!l which dir to put stuff into????')
    end
  end

  fnamex = ['/strowdataN/s3/robinso2/PUYEHUE/CLOUDS_6chans/TYPEIS_rCfracType=-1_x10/eracld_type0_allfov_'];
  fnamex = [fnamex ystr '_' mstr '_' dstr '_' gstr '.rtp'];
  ee = exist(fnamex);
  if ee > 0
    iCnt = iCnt + 1;
    [h,ha,p,pa] = oldrtpread(fnamex);
    fprintf(1,'%s \n',fnamex);
    if iCnt == 1
      hall = h; pall = p;
    else
     if (sum(abs(h.ichan - hall.ichan)) ~= 0)
        error('kkkkkkk')
      else
        h.vchan = hall.vchan;
        [hall,pall] = cat_rtp(hall,pall,h,p);
      end
    end
  end
end

jett = jet; jett(1,:) = 1;

water= find(pall.ctype == 101 | pall.ctype2 == 101);
simplemap(pall.rlat(water),pall.rlon(water),log10(pall.cngwat(water))); 
caxis([-2 3]); colormap(jett); colorbar
title('log10(wateramt)')                                                        

ice= find(pall.ctype == 201 | pall.ctype2 == 201);                              
simplemap(pall.rlat(ice),pall.rlon(ice),log10(pall.cngwat(ice))); colorbar
caxis([-2 3]); colormap(jett); colorbar
title('log10(iceamt)')                                                    
