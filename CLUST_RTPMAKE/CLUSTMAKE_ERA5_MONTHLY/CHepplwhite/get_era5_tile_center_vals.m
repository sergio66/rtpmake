function [scen, smean, stdv] = get_era5_tile_center_vals(era)

%
% first run: era5 = load_era5_monthly()
%
%

% check presence of fields
allfields = {'inda','indn','tz','q','clwc','ciwc','cc','skt',...
              'ednum','sdnum','lat','lon'};
if( ~all(ismember(allfields,fieldnames(era))) )
  disp('Error: not all fields present');
end

% center bin (original era and local array). lat=a, lon=n.
scen.inda = floor(0.5*(era.inda(1) + era.inda(end) ));
scen.bina = floor(length(era.inda)/2);
scen.indn = ceil(0.5*(era.indn(1) + era.indn(end) ));
scen.binn = ceil(length(era.indn)/2);

% center mean and std.dev for each field
scen.tz      = squeeze(era.tz(scen.binn, scen.bina, :, :));
scen.q       = squeeze(era.q(scen.binn, scen.bina, :, :));
scen.o3      = squeeze(era.o3(scen.binn, scen.bina, :, :));
scen.clwc    = squeeze(era.clwc(scen.binn, scen.bina, :, :));
scen.ciwc    = squeeze(era.ciwc(scen.binn, scen.bina, :, :));
scen.cc      = squeeze(era.cc(scen.binn, scen.bina, :, :));
scen.skt     = squeeze(era.skt(scen.binn, scen.bina, :));
%
smean.tz     = squeeze(nanmean(era.tz(:, :, :, :),[1 2]));
smean.q      = squeeze(nanmean(era.q(:, :, :, :),[1 2]));
smean.o3     = squeeze(nanmean(era.o3(:, :, :, :),[1 2]));
smean.clwc   = squeeze(nanmean(era.clwc(:, :, :, :),[1 2]));
smean.ciwc   = squeeze(nanmean(era.ciwc(:, :, :, :),[1 2]));
smean.cc     = squeeze(nanmean(era.cc(:, :, :, :),[1 2]));
smean.skt    = squeeze(nanmean(era.skt(:, :, :),[1 2]));
%
stdv.tz     = squeeze(nanstd(era.tz(:, :, :, :),0,[1 2]));
stdv.q      = squeeze(nanstd(era.q(:, :, :, :),0,[1 2]));
stdv.o3     = squeeze(nanstd(era.o3(:, :, :, :),0,[1 2]));
stdv.clwc   = squeeze(nanstd(era.clwc(:, :, :, :),0,[1 2]));
stdv.ciwc   = squeeze(nanstd(era.ciwc(:, :, :, :),0,[1 2]));
stdv.cc     = squeeze(nanstd(era.cc(:, :, :, :),0,[1 2]));
stdv.skt    = squeeze(nanstd(era.skt(:, :, :, :),0,[1 2]));
%
scen.lon     = era.lon(scen.indn);
scen.lat     = era.lat(scen.inda);
%
if(~isequal(era.sdnum, era.ednum))
  disp('warning: layer and sfc times differ')
end
