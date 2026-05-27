% Select a critical Richardson number. While 0.25 is common, studies
% suggest using 0.24 for strong stable layers, 0.31 for weak stable,
% and 0.39 for unstable conditions for better accuracy.

%% iVers_Ri = 0;
%% 
%% see J.Clim https://journals.ametsoc.org/view/journals/clim/31/22/jcli-d-17-0498.1.xml?tab_body=pdf
%%   The Climatology of the Atmospheric Boundary Layer in Contemporary Global Climate Models
%%   Richard Davy
%%   Print Publication: 15 Nov 2018
%%   DOI: https://doi.org/10.1175/JCLI-D-17-0498.1
%%   Page(s): 9151–9173
%%
%% also see Climatology of the planetary boundary layer over the
%%   continentalUnited States and EuropeDian J. Seidel, 1 Yehui Zhang, 2
%%   Anton Beljaars, 3 Jean-Christophe Golaz, 4Andrew R. Jacobson, 5 and
%%   Brian Medeiros
%% JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 117, D17106, doi:10.1029/2012JD018143, 2012
%%
%% convert to Ri(z) =  g         (Tpot(z)-Tpot0)*(z-z0)
%%                    ---   -------------------------------
%%                    Tpot0   (horizspeed(z) - surfspeed)^2
%%
%% where Tpot0,zpot0 = potential temp at surface, surface altitude
%% units m/s2 K m / (K m2/s2) = m2/s2/(m2/s2) = []
%%
%%
%%      In this method the PBL depth is defined to be
%%      the height at which the bulk Richardson number first
%%      exceeds some critical value, which we took to be 0.25
%%      based upon observational evidence of the threshold
%%      for Kolmogorov turbulence (Grachev et al. 2013). We
%%      tested the sensitivity of the PBL depth to this critical
%%      value and found variations in climatological-mean PBL
%%      depths in individual models of less than 3% in response
%%      to 20% variations in the critical value from 0.2 to
%%      0.3, which is the range of likely values for the critical
%%      Richardson number from analysis of observations (Cheng
%%      et al. 2002; Zilitinkevich and Esau 2007).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% iVers_Ri = 1;
%% 
%% Evaluation of the Planetary Boundary Layer Height From ERA5 Reanalysis With MOSAiC Observations Over the Arctic Ocean
%% Xingya Xi, Qinghua Yang, Changwei Liu, Matthew D. Shupe, Bo Han, Shijie Peng, Shaohui Zhou, Dake Chen
%% JGR Atmospheres : 22 June 2024 https://doi.org/10.1029/2024JD040779
%% https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024JD040779
%% 
%% In the ECMWF Integrated Forecasting System (IFS), which generates the
%% ERA5 reanalysis, the bulk Richardson number Ri_{b} profile
%% calculation depends on the height z above the surface.
%% 
%% It relates the buoyant destruction or production of turbulence to the
%% mechanical generation of turbulence by wind shear.The exact equation
%% used is defined as:
%%   Ri_b(z) = g        S_v(z) - S_v(z_s)
%%            --      --------------------
%%          Svlayer       U(z)^2
%% 
%% Where the parameters represent:
%% g      : Acceleration due to gravity ~ 9.80665 m/s2
%% S_v(z) : Virtual dry static energy at height z\ calculated as S_v = c_p T_v + gz,
%% T_{v}  : the virtual temperature
%% c_{p}  : pecific heat of dry air,
%% z_{s}  : Height of the model's lowest prognostic level (typically around 10 m)
%% S_{v_{layer}} : average virtual dry static energy across ayer between z and z_{s}
%% U(z)   : Magnitude of the horizontal wind vector at height z relative to the lowest model
%%            level (i.e., \Delta U = |f{U}(z) - {U}(z_s)|.
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RiCritical = 0.20;  %% set this empirically, did not work

RiCritical = 0.27;   %% I used this initially, looks like it is a little too high

RiCritical = 0.20;   %% try this since I am still 25% too high (UMBC_PBLH = 1.25 ERA5_PBLH, but I still think this is not low enough

RiCritical = 0.25;   %% ERA5 uses this, but they also have a surface friction velocity term in denom   u^2 + v^2 + b ustar^2
                     %% where ustar = surfae friction velocity, b is a coeff accoounting for this (so effectively I have b=0)
                     %% results in 0.06 km lowering of PBLH over ocean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iVers_Ri = 0;    %% simple one,  from Richard Davy, J. Clim, 2018 and Dian Siedel,  JGR-Atmos, 2012, uses potential temp
iVers_Ri = 1;    %% fancier one, from Xingya Xi , JGR_Atmos 2024 (and ERA5)                          uses staticE
%iVers_Ri = 2;    %% fancier one, from Xingya Xi , JGR_Atmos 2024 (and ERA5)                          uses staticE, and with diff(numer)/diff(denom) rather than bulk

