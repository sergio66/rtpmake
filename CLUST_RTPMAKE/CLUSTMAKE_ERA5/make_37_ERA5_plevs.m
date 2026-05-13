%{
%% https://climate.copernicus.eu/climate-reanalysis
ERA5 and ERA5-Land monthly mean 2m temperature Jan 2016

ERA5 is the latest climate reanalysis produced by ECMWF, providing
hourly data on many atmospheric, land-surface and sea-state parameters
together with estimates of uncertainty.

ERA5 data are available in the Climate Data Store on regular
latitude-longitude grids at 0.25o x 0.25o resolution, with atmospheric
parameters on 37 pressure levels.

ERA5 is available from 1940 and continues to be extended forward in
time, with daily updates being made available 5 days behind real time

Initial release data, i.e., data no more than three months behind real
time, are called ERA5T.

The products of the reanalysis are available to the public through the
Climate Data Store.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
https://ams.confex.com/ams/104ANNUAL/meetingapp.cgi/Paper/439790

AMS 103rd Annual Meeting
About the Meeting

Carlos Mario Cuervo López (Presenter)
Central Michigan University
Mount Pleasant, MI
USA

John T. Allen
Central Michigan Univ.
Mount Pleasant, MI
USA

Mateusz Taszarek
NSSL
Norman, OK
USA

Visit Ams Website

70 - Comparative Analysis of ERA5 Model Levels and Pressure Levels Over the Continental U.S.

Reanalyses data is recurrently used as a surrogate for observed
atmospheric data in numerous scientific applications, representing
robust historical environmental profiles derived through data
assimilation. Nevertheless, its performance and suitability for
different atmospheric phenomena have been seldom assessed. Among all
reanalyses, the fifth generation of the European Centre for
Medium-Range Weather Forecasts (ECMWF) reanalysis (ERA5) stands out as
a favored choice for downscaling, case studies, and convective
analysis, in part due to its high temporal and spatial
resolution. ERA5 is used frequently in two vertical level versions:
model level (ML) and pressure level (PL) data. The ML consists of 137
hybrid sigma-pressure vertical levels, 20 of which are on the lowest
kilometer, while the PL consists of 37 pressure levels ranging from
1000 hPa (near surface) to 1 hPa (about 80 km). Additionally, the
accessibility to PL data drives its preference among users rather than
higher-resolution ML products. This is particularly relevant in
applications for dynamic downscaling, case applications, or studies
using convective parameters. Comparing the ML and PL data allows us to
assess these two versions’ relative performance and biases, providing
essential insights into the accuracy and applicability of each ERA5
version.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% https://forums.jcsda.org/t/about-era5-pressurelevels-and-pressurelayers/91
era5plevs = [1 2 3 5 7 10 20 30 50 70 100 125 150 175 200 225 250 300 350 400 450 500 550 600 650 700 750 775 800 825 850 875 900 925 950 975 1000];

comment = 'see make_37_ERA5_plevs.m';
save era5plevs.mat era5plevs comment
