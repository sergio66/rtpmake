%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

klayers  = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_80km/klayersV205_0_80km/BinV221/klayers_airs';
sartaCld = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
sartaCld = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORIG, before Nov 2025
%% addpath /home/sergio/MATLABCODE
%% addpath /asl/matlab2012/airs/readers
%% addpath /asl/matlib/aslutil
%% %addpath /asl/matlib/science
%% addpath /home/sergio/MATLABCODE/matlib/science/
%% addpath /asl/matlib/rtptools
%% addpath /asl/matlib/h4tools/
%% addpath /asl/matlib/rtptools/
%% addpath /asl/matlib/gribtools/
%% addpath /asl/matlib/time
%% addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
%% addpath /home/sergio/MATLABCODE
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis
%% 
%% % addpath /home/strow/cress/Work/Rtp
%% % addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
%% % addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
%% % addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
%% % addpath /asl/packages/rtp_prod2/grib
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/grib
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util
%% addpath /home/sergio/MATLABCODE/NANROUTINES

%% addpath /home/sergio/MATLABCODE/TIME
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
%% 
%% %addpath /asl/rtp_prod2/emis/
%% %addpath /asl/rtp_prod2/util/    
%% %addpath /asl/packages/rtp_prod2/emis/
%% %addpath /asl/packages/rtp_prod2/util/

%    %addpath /asl/rtp_prod2/emis/
%    %addpath /asl/rtp_prod2/util/    
%    %addpath /asl/packages/rtp_prod2/emis/
%    %addpath /asl/packages/rtp_prod2/util/
%    %addpath /asl/rtp_prod2/emis/
%    %addpath /asl/rtp_prod2/util/

%    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
%    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jan 2026

addpath /home/sergio/git/matlabcode/
addpath /home/sergio/git/matlabcode/PLOTTER
addpath /home/sergio/git/matlabcode/NANROUTINES

addpath /home/sergio/git/matlabcode
addpath /home/sergio/git/matlabcode/matlibSergio/matlab2012/airs/readers
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/science/

%% these two are now the same
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/clouds/sarta
addpath /home/sergio/git/sergio_matlib/matlib/clouds/sarta/

%% addpath /home/strow/cress/Work/Rtp
%% addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
%% addpath /home/sergio/git/matlabcode/CRIS_HiRes             %% for sergio_fill_ecmwf
%% addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
%% addpath /asl/packages/rtp_prod2/grib

%addpath /asl/rtp_prod2/emis/
%addpath /asl/rtp_prod2/util/    
%addpath /asl/packages/rtp_prod2/emis/
%addpath /asl/packages/rtp_prod2/util/

addpath /home/sergio/git/rtp_prod2/emis

%%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2/emis    %%%% >>>>>>>>
%%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2/grib    %%%% >>>>>>>>
%%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2/util    %%%% >>>>>>>>
%% could try these???
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2_sbuczkowski_nogit/grib
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2_Aug11_2020/emis
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2_Aug11_2020/grib
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2_Aug11_2020/util

addpath /home/sergio/git/matlabcode/TIME
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/aslutil
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/science
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/h4tools
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtptools
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/gribtools
%% addpath /home/sergio/git/matlabcode/matlibSergio/matlib/TIME

addpath ../GRIB
addpath /home/sergio/git/matlabcode/PLOTTER
addpath /home/sergio/git/matlibcode/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/
