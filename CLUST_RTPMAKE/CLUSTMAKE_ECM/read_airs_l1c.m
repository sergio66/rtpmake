function [l1c] = read_airs_l1c(fn);

% cp ~chepplew/myLib/matlib/read_airs_l1c.m .

% Granule dimensions
nchan=2645;
nxtrack=90;
natrack=135;
nobs=nxtrack*natrack;


l1c.Latitude    = hdfread(fn,'Latitude');
l1c.Longitude   = hdfread(fn,'Longitude');
l1c.Time        = hdfread(fn,'Time');
l1c.radiances   = hdfread(fn,'radiances');
l1c.scanang     = hdfread(fn,'scanang');
l1c.satzen      = hdfread(fn,'satzen');
l1c.satazi      = hdfread(fn,'satazi');
l1c.solzen      = hdfread(fn,'solzen');
l1c.solazi      = hdfread(fn,'solazi');
l1c.sun_glint_distance = hdfread(fn,'sun_glint_distance');
l1c.topog       = hdfread(fn,'topog');
l1c.topog_err   = hdfread(fn,'topog_err');
l1c.landFrac    = hdfread(fn,'landFrac');
l1c.landFrac_err = hdfread(fn,'landFrac_err');
l1c.ftptgeoqa   = hdfread(fn,'ftptgeoqa');
l1c.zengeoqa    = hdfread(fn,'zengeoqa');
l1c.demgeoqa    = hdfread(fn,'demgeoqa');
l1c.state       = hdfread(fn,'state');
l1c.Rdiff_swindow = hdfread(fn,'Rdiff_swindow');
l1c.Rdiff_lwindow = hdfread(fn,'Rdiff_lwindow');
l1c.SceneInhomogeneous = hdfread(fn,'SceneInhomogeneous');
l1c.dust_flag   = hdfread(fn,'dust_flag');
l1c.dust_score  = hdfread(fn,'dust_score');
l1c.spectral_clear_indicator = hdfread(fn,'spectral_clear_indicator');
l1c.BT_diff_SO2 = hdfread(fn,'BT_diff_SO2');
l1c.AB_Weight   = hdfread(fn,'AB_Weight');
l1c.L1cProc     = hdfread(fn,'L1cProc');
l1c.L1cSynthReason = hdfread(fn,'L1cSynthReason');
l1c.NeN         = hdfread(fn,'NeN');
l1c.Inhomo850   = hdfread(fn,'Inhomo850');
l1c.nadirTAI    = cell2mat(hdfread(fn,'nadirTAI'));
l1c.satheight   = cell2mat(hdfread(fn,'satheight'));
l1c.sat_lat     = cell2mat(hdfread(fn, 'sat_lat'));
l1c.sat_lon     = cell2mat(hdfread(fn, 'sat_lon'));
l1c.node_type   = cell2mat(hdfread(fn, 'node_type'));
l1c.scan_node_type  = cell2mat(hdfread(fn, 'scan_node_type'));

%{
% open swath ID to get more data
file_id    = hdfsw('open',fn,'read');
swath_id   = hdfsw('attach',file_id,'L1C_AIRS_Science');
[junk,s]   = hdfsw('readattr',swath_id,'granule_number');
tmp_granule_id = double(junk(1));

[junk,s]   = hdfsw('readfield',swath_id,'nominal_freq',[],[],[]);
l1c.freq   = junk;

[junk,s]   = hdfsw('readfield',swath_id,'ChanID',[],[],[]);
l1c.chanID = junk;

[junk,s]   = hdfsw('readattr',swath_id,'eq_x_tai');
l1c.eq_x_tai = junk;

% close hdf file
s = hdfsw('detach',swath_id);
if s == -1; disp('Swatch detach error: L1b');end;   
s = hdfsw('close',file_id);
if s == -1; disp('File close error: L1b');end;

%}

tmp_gran_id  = hdfread(fn, 'granule_number');

l1c.eq_x_tai = hdfread(fn, 'eq_x_tai');

l1c.freq     = hdfread(fn, 'nominal_freq');

l1c.chanID   = hdfread(fn, 'ChanID');

tmp_state = reshape( l1c.state, 1,nobs);
i0=find( tmp_state <= 1);  % Indices of "good" FOVs (Process or Special)
n0=length(i0);


% Declare temporary variables for expansion
tmp_atrack=zeros(natrack, nxtrack);
tmp_xtrack=zeros(natrack, nxtrack);

% Loop over along-track and fill in temporary variables
for ia=1:natrack
  for ix = 1:nxtrack
    tmp_xtrack(ia, ix) = ix;
    tmp_atrack(ia, ix) = ia;
  end
end

% junk = reshape(tmp_atrack,natrack,nxtrack);
l1c.atrack = single(tmp_atrack);   % was tmp_atrack(i0)
l1c.xtrack = single(tmp_xtrack);
l1c.gindex = single(double(tmp_gran_id{1})*ones(natrack, nxtrack));

return
