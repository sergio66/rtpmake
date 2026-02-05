%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this was a test run, Jan 2023 +/- 3 months
setstr = '/SET1/';
clear dir0
dir0{1}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s456/N68p25/';
dir0{2}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s457/N68p25/';
dir0{3}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s458/N68p25/';
dir0{4}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s459/N68p25/';
dir0{5}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s460/N68p25/';
dir0{6}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s461/N68p25/';
dir0{7}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s462/N68p25/';
dir0{8}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s463/N68p25/';
dir0{9}  = '/umbc/rs/strow/asl/airs/tile_test7/2022_s464/N68p25/';
dir0{10} = '/umbc/rs/strow/asl/airs/tile_test7/2022_s465/N68p25/';
dir0{11} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s466/N68p25/';
dir0{12} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s467/N68p25/';
dir0{13} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s468/N68p25/';
dir0{14} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s469/N68p25/';
dir0{15} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s470/N68p25/';
dir0{16} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s471/N68p25/';
dir0{17} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s472/N68p25/';
dir0{18} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s473/N68p25/';
dir0{19} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s474/N68p25/';
dir0{20} = '/umbc/rs/strow/asl/airs/tile_test7/2023_s475/N68p25/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode
%% addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/TIME
%% a = read_netcdf_lls('/umbc/rs/strow/asl/airs/tile_test7/2020_s401/N57p75/tile_2020_s401_N57p75_W105p00.nc'); %% Plot A2 of JGR 2025 paper, Hudson Bay, 59N,-102W on 3/18/2020
%% [yy,mm,dd] = tai2utcSergio(a.tai93(1)+offset1958_to_1993);
setstr = '/SET2/';
clear dir0
dir0{1}  = '/umbc/rs/strow/asl/airs/tile_test7/2019_s392/N57p75/';
dir0{2}  = '/umbc/rs/strow/asl/airs/tile_test7/2019_s393/N57p75/';
dir0{3}  = '/umbc/rs/strow/asl/airs/tile_test7/2019_s394/N57p75/';
dir0{4}  = '/umbc/rs/strow/asl/airs/tile_test7/2019_s395/N57p75/';
dir0{5}  = '/umbc/rs/strow/asl/airs/tile_test7/2019_s396/N57p75/';
dir0{6}  = '/umbc/rs/strow/asl/airs/tile_test7/2020_s397/N57p75/';
dir0{7}  = '/umbc/rs/strow/asl/airs/tile_test7/2020_s398/N57p75/';
dir0{8}  = '/umbc/rs/strow/asl/airs/tile_test7/2020_s399/N57p75/';
dir0{9}  = '/umbc/rs/strow/asl/airs/tile_test7/2020_s400/N57p75/';
dir0{10} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s401/N57p75/'; %%% --->>> this is the event
dir0{11} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s402/N57p75/';
dir0{12} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s403/N57p75/';
dir0{13} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s404/N57p75/';
dir0{14} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s405/N57p75/';
dir0{15} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s406/N57p75/';
dir0{16} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s407/N57p75/';
dir0{17} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s408/N57p75/';
dir0{18} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s409/N57p75/';
dir0{19} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s410/N57p75/';
%dir0{20} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s411/N57p75/'; %% DNE 
dir0{20} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s412/N57p75/';
dir0{21} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s413/N57p75/';
dir0{22} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s414/N57p75/';
dir0{23} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s415/N57p75/';
dir0{24} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s416/N57p75/';
dir0{25} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s417/N57p75/';
dir0{26} = '/umbc/rs/strow/asl/airs/tile_test7/2020_s418/N57p75/'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
