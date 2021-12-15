function [fits] = get_era5_tile_rates(era, scen, smean)

% INPUTS:
%       era:   structure of loaded era5 monthly values
%       scen:  structure of center values
%       smean: structure of mean values
%
% OUTPUT: fits: structure with two members
%   fits.cen;  for center value of tile
%   fits.mean; for mean of all era values in the tile
%   each member has fields for variables, skt, T(z), O3(z), Q(z)
%        clwc(z), ciwc(z), cc(z)
%   each variable includes the fir coefficients, lag-1 error
%        and anomaly
%

% check input arguments
if(nargin ~= 3)
  error('Need 3 inputs')
end

% check required fields are present
% 1 - center values from structure cen
cenflds = {'skt','tz','o3','q','clwc','ciwc','cc'};
if( ~all(ismember(cenflds, fieldnames(scen))) )
  disp('Error: not all center fields are present')
end

% 2 - mean values from struct: mean
mnflds = cenflds;
if( ~all(ismember(mnflds, fieldnames(smean))) )
  disp('Error: not all average fields are present')
end

% subset time period
% not used
indt = find(era.sdnum);

  warning('off','all')

disp('Fiting for center values in tiles')


% 1-i) cen skt trend
if(isfield(scen,'skt'))

  disp('fitting skt')

  xdat  = era.sdnum(indt);
  ydat = squeeze(scen.skt(indt)) - nanmean(scen.skt(indt));
    try
      [B, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,2);
      fits.cen.skt_bcoef = B;
      %yhat = B(2)*(days(ndt1:ndt2)-days(ndt1))./365 + B(1);
      %linr(i, j,:) = B(2);
% Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      fits.cen.skt_rerr   = Bstats.se(2).*lag1;
% to compare with SM anom retrievals add back linear term
      dr = (xdat - xdat(1)+1)./365. *B(2);
      fits.cen.skt_anom = Bstats.resid + dr;
    catch ME
      fprintf(1, '%s\n', ME.message);
      %continue;
    end

end   % end if isfield(skt)

% ---------------------------
% 1-ii) cen levels field

warning('off','all')
%olinfit = zeros(1440,721,37,1);
%anom    = [];
%fiterr  = zeros(1440,721,37,1);
indt    = find(era.ednum);

for i = 2:length(cenflds)
  disp(['fitting ' str2mat(cenflds(i))])
  cfldname = [(cell2mat(cenflds(i))) '_coef'];
  afldname = [(cell2mat(cenflds(i))) '_anom'];
  rfldname = [(cell2mat(cenflds(i))) '_rerr'];
  tfldn    = str2mat(cenflds(i));

  xdat = era.ednum(indt);
  
  for ilev = 1:length(era.plevs)
    ydat = squeeze(scen.(tfldn)(ilev,indt)) - nanmean(scen.(tfldn)(ilev,indt),2);
    try
      [B, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,2);
      fits.cen.(cfldname)(ilev,:) = B;
      %olinfit(i, j, ilev,:) = B(2);
      % to compare with SM anom retrievals add back linear term
      % Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      fits.cen.(rfldname)(ilev,:) = Bstats.se(2).*lag1;
      %
      dr = (xdat - xdat(1)+1)./365. *B(2);
      fits.cen.(afldname)(ilev,:) = Bstats.resid + dr;
    catch ME
      fprintf(1, '%s\n', ME.message);
      continue;
    end

  end   % end plevs
end     % end cen.fields      

% ---------------------------------------------
%  Fit Mean of era values in the selected tile
% ---------------------------------------------

disp('Fiting for means of values in tiles')

% 1-i) mean skt trend
if(isfield(smean,'skt'))

  disp('fitting skt')
  xdat  = era.sdnum(indt);
  warning('off','all')
  ydat = squeeze(smean.skt(indt)) - nanmean(smean.skt(indt));
    try
      [B, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,2);
      fits.mean.skt_bcoef = B;
      %yhat = B(2)*(days(ndt1:ndt2)-days(ndt1))./365 + B(1);
      %linr(i, j,:) = B(2);
% Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      fits.mean.skt_rerr   = Bstats.se(2).*lag1;
% to compare with SM anom retrievals add back linear term
      dr = (xdat - xdat(1)+1)./365. *B(2);
      fits.mean.skt_anom = Bstats.resid + dr;
    catch ME
      fprintf(1, '%s\n', ME.message);
      %continue;
    end
end   % end if isfield(skt)

% ---------------------------
% 1-ii) cen levels field

warning('off','all')
%olinfit = zeros(1440,721,37,1);
%anom    = [];
%fiterr  = zeros(1440,721,37,1);
indt    = find(era.ednum);

for i = 2:length(mnflds)
  disp(['fitting ' str2mat(mnflds(i))])
  cfldname = [(cell2mat(mnflds(i))) '_coef'];
  afldname = [(cell2mat(mnflds(i))) '_anom'];
  rfldname = [(cell2mat(mnflds(i))) '_rerr'];
  tfldn    = str2mat(mnflds(i));

  xdat = era.ednum(indt);
  
  for ilev = 1:length(era.plevs)
    ydat = squeeze(smean.(tfldn)(ilev,indt)) - nanmean(smean.(tfldn)(ilev,indt),2);
    try
      [B, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,2);
      fits.mean.(cfldname)(ilev,:) = B;
      %olinfit(i, j, ilev,:) = B(2);
      % to compare with SM anom retrievals add back linear term
      % Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      fits.mean.(rfldname)(ilev,:) = Bstats.se(2).*lag1;
       %
      dr = (xdat - xdat(1)+1)./365. *B(2);
      fits.mean.(afldname)(ilev,:) = Bstats.resid + dr;
    catch ME
      fprintf(1, '%s\n', ME.message);
      continue;
    end

  end   % end plevs
end     % end mean.fields      

  warning('on','all')

% prep final structure to output

fits.cen.type  = 'center';
fits.mean.type = 'average of era over airs tile';


% ------ saving -------

savdr = '/home/chepplew/data/rates_anomalies/tiled/era5_mon/fits/';
savfn = ['era5_tile_fits_lonbin'  ...
         sprintf('%02d_latbin%02d',era.lonbin,era.latbin) '.mat'];

if(~exist(savdr))
  disp(['Creating: ' savdr]);
  mkdir(savdr);
end

disp(['saving fits to: ' [savdr savfn]]);
save([savdr savfn],'fits','xdat');

