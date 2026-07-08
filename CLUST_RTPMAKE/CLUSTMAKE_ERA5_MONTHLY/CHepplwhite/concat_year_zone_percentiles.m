function [] = concat_year_zone_percentiles()

% function [] = concat_year_zone_percentiles()
%
%   concatenate a given zone percentile file into single
% time series
%   pctr[23 x 72 x 2645 x 100]
%   dnum[23 x 72]
%


% Check zone string
zn_str = 'S90p00';

d.home = '/home/chepplew/data/rates_anomalies/tiling/percentiles/';
d.dir  = dir([d.home '2*_' zn_str '_percentiles.mat']);

nt  = 1;      % time step counter

ifn = 1;
  x = load([d.dir(ifn).folder,'/' d.dir(ifn).name]);

  nsz = size(x.pctr);
  all_pctr(nt:nt+nsz(1)-1,:,:,:) = x.pctr;
  all_dnum(nt:nt+nsz(1)-1,:)     = x.dnum;
  nt = nt+nsz(1);

end
