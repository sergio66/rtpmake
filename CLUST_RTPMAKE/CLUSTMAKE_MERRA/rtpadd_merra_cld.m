function [head hattr prof pattr] = rtpadd_merra(head,hattr,prof,pattr,root,iFix)
% function [head hattr prof pattr] = rtpadd_merra(head,hattr,prof,pattr,root,iFix)
%
% root - usuall '/asl' (optional)
% iFix - if iFix = 1, get rid of NANs in profile (default); if -1, keep NaNs as they are
%      - this is optional arg
%
% Add MERRA Model into the given RTP structure
%
% Paul Schou, Breno Imbiriba  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fetch MERRA Data 
  % Fields 	description	RTP field name
  % t  - air temperature  	ptemp
  % qv - specific humidity	gas_1
  % o3 - ozone mixing ratio	gas_3
  % ps - surface pressure		spres
  % ts - surface temperature	stemp
  % u2m - eastward wind at 2 m above the displacement height
  % v2m - northward wind at 2 m above the displacement height
  % ciwc,clwc - cloud mixing ratio
  %
  % Data is given on 3hr files, except 'ts' which is hourly
  %

  if(nargin()==4)
    root = '/asl';
    iFix = +1;
  end

  if(nargin()==5)
    iFix = +1;
  end

  % Assumptions:
  % 1. All pressure grids are the same for all 3D variables
  % 2. Invalid bottom of atmosphere values happens on all 3D variables


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Basic checks
  if(isfield(head,'ptype') && head.ptype~=0)
    disp('You have an RTP file that is a level file. All previous layer information will be removed');
    while head.ngas>0
      if(isfield(prof,['gas_' num2str(head.glist(1))]))
	prof=rmfield(prof,['gas_' num2str(head.glist(1))]);
	head.ngas  = head.ngas-1;
	head.glist = head.glist(2:end,1);
	head.gunit = head.gunit(2:end,1);
	head.ptype = 0;
      else
	warning(['Non existing field gas_' num2str(head.glist(1)) ' indicated by headers. Fixing']);
	head.ngas  = head.ngas-1;
	head.glist = head.glist(2:end,1);
	head.gunit = head.gunit(2:end,1);
	head.ptype = 0;
      end
    end
  end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get the data times:

  % Compute FoVs matlab times and the reference year for TAI (used later)
  % 
  [mtimes mreftime] = rtpdate(prof,pattr);
  mrefyear = str2num(datestr(mreftime,'yyyy'));


  threehours = round((mtimes-floor(min(mtimes)))*8); % [threehours]=3-hour long units
  u3hours = unique(threehours); % unique list of the used 3-hour intervals
  u3hours(isnan(u3hours))=[]; % Remove NaNs that may come from mtimes==NaN
  n3hours = numel(u3hours);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Loop over 
  
  for i3hours=1:n3hours

    % Select Fovs and subset them
    ifov=find(threehours==u3hours(i3hours));
    nfovs = numel(ifov);
    tprof = ProfSubset2(prof,ifov);

    % Get required 3hr time slot
    reqtime = floor(min(mtimes)) + u3hours(i3hours)/8; % [day]=[day]+[3hr]/8


    %%%%%%%%%%%%%%%%%%%% 
    % Set profile variables

    ptemp=[]; pgas_1=[]; pgas_3=[]; pps=[]; pts=[]; ciwc=[]; clwc=[];


    %%%%%%%%%%%%%%%%%%%% 
    % Interpolate 3D variables for each layer
    % ATTENTION: Fill value is not consistent. Some times it is 1e15 (as
    % advertised, but some times it is -9.99e8. Go figure!

    %% http://wiki.seas.harvard.edu/geos-chem/index.php/List_of_MERRA_met_fields kg kg-1
    %% cc - cloud cover
%    [dat_cc plevs lats lons]= getdata_merra(reqtime, 'cloud',[],root);

    % ciwc  - ice cloud
    [dat_ciwc plevs lats lons]= getdata_merra(reqtime, 'qi',[],root);

    % clwc  - liquid cloud
    [dat_clwc plevs lats lons]= getdata_merra(reqtime, 'ql',[],root);
 
    % t  - air temperature  	ptemp
    [dat_t plevs lats lons]= getdata_merra(reqtime, 't',[],root);


    %%%%%%%%%%%%%%%%%%%%
    % ptime is computed by reverting mattime to TAI using the reference year above
    reqTAI = mattime2tai(reqtime, mrefyear);
    tprof.ptime = ones(1,nfovs).*reqTAI; % [sec]


    %%%%%%%%%%%%%%%%%%%%
    % Field quantities. Use Nearest Match. 

    [glats glons] = meshgrid(lats,lons);
    tprof.plat = single(interp_sphere(lats, lons, glats, tprof.rlat, tprof.rlon, 'nearest'));
    tprof.plon = single(interp_sphere(lats, lons, glons, tprof.rlat, tprof.rlon, 'nearest'));


    nlevs=numel(plevs);

%    dat_cc(dat_cc>1e14 | dat_cc<-1)=NaN;
%    for ilev=1:nlevs
%      cc(ilev,:) = interp_sphere(lats, lons, dat_cc(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
%    end

    dat_ciwc(dat_ciwc>1e14 | dat_ciwc<-1)=NaN;
    for ilev=1:nlevs
      ciwc(ilev,:) = interp_sphere(lats, lons, dat_ciwc(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
    end

    dat_clwc(dat_clwc>1e14 | dat_clwc<-1)=NaN;
    for ilev=1:nlevs
      clwc(ilev,:) = interp_sphere(lats, lons, dat_clwc(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
    end

    dat_t(dat_t>1e14 | dat_t<-1)=NaN;
    for ilev=1:nlevs
      ptemp(ilev,:) = interp_sphere(lats, lons, dat_t(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
    end
    
    % qv - specific humidity	gas_1
    [dat_q plevs lats lons]= getdata_merra(reqtime, 'qv',[],root);
    dat_q(dat_q>1e14 | dat_t<-1)=NaN;
    for ilev=1:nlevs
      pgas_1(ilev,:) = interp_sphere(lats, lons, dat_q(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
    end
   
    % o3 - ozone mixing ratio	gas_3
    [dat_o3 plevs lats lons]= getdata_merra(reqtime, 'o3',[],root);
    dat_o3(dat_o3>1e14 | dat_t<-1)=NaN;
    for ilev=1:nlevs
      pgas_3(ilev,:) = interp_sphere(lats, lons, dat_o3(:,:,nlevs-ilev+1), tprof.rlat, tprof.rlon,'nearest');
    end

    t_n = numel(plevs);
    plevs = reshape(plevs,[t_n,1]);
    plevs = plevs(t_n:-1:1,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Compute Valid Levels
   
    % For each profile, find the last valid level. 
    % 
    tprof.plevs = single(plevs*ones(1, nfovs));
    tprof.nlevs = single(nlevs.*zeros(1, nfovs));
    tprof.nlevs(1,:) = single(nlevs-sum(isnan(ptemp)));

%    tprof.cc    = single(cc);
    tprof.clwc  = single(clwc);
    tprof.ciwc  = single(ciwc);
    tprof.ptemp = single(ptemp);
    tprof.gas_1 = single(pgas_1);
    tprof.gas_3 = single(pgas_3);


    % ps - surface pressure	spres
    [dat_ps tmpx lats lons]= getdata_merra(reqtime, 'ps',[],root); % It is in Pa, convert to mbar -> /100
    dat_ps(dat_ps>1e14)=NaN;
    pps = interp_sphere(lats,lons,dat_ps/100, tprof.rlat, tprof.rlon, 'nearest');
    
    % ts - surface temperature	stemp
    [dat_ts tmpx lats lons merra_str]= getdata_merra(reqtime, 'ts',[],root);
    dat_ts(dat_ts>1e14)=NaN;
    pts = interp_sphere(lats,lons,dat_ts, tprof.rlat, tprof.rlon, 'nearest');
   
    tprof.spres = single(pps);
    tprof.stemp = single(pts);  

    
    % wind speed at 2m
    [dat_u2m tmpx lats lons]= getdata_merra(reqtime, 'u2m',[],root);
    [dat_v2m tmpx lats lons]= getdata_merra(reqtime, 'v2m',[],root);
    dat_w2m = sqrt(dat_u2m.^2 + dat_v2m.^2);
    w2m = interp_sphere(lats,lons,dat_w2m,tprof.rlat, tprof.rlon,'nearest'); 

    tprof.wspeed = single(w2m);


    if(n3hours>1)
      tprof_arr(i3hours) = tprof;
    end

  end

  if(n3hours>1)
    tprof = Prof_join_arr(tprof_arr);
    clear tprof_arr
  end 
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fix Header and Attributes


  % Set that the profile type is level.
  head.ptype=0;

  % Set the profile bit
  [ia ib ic]=pfields2bits(head.pfields);
  head.pfields = bits2pfields(1,ib,ic);

  % Add the new two gases - 1 and 3
  for ig=[1 3]
    if(~isfield(head,'glist'))
      head.glist=[];
      head.gunit=[];
      head.ngas=0;
    end
    ik=find(head.glist==ig); 
    if(numel(ik)==0)
      head.glist(end+1,1)=ig; % gad id
      head.gunit(end+1,1)=21; % gas unit (g/g) or something like that
      head.ngas=head.ngas+1;
    else
      head.gunit(ik)=21;
    end 
  end

  % Set pressures
  head.pmin = min(tprof.plevs(:));
  head.pmax = max(tprof.spres(:));

  hattr = set_attr(hattr,'profile',[merra_str ' Nearest']);

  prof=tprof;

  %% fix up the NaNs over land, by (a) putting constant gas, cloud profiles to the bottom and
  %%                               (b) linear interpolation to surface pressure, temperature
  % iFix = +1;   %% clean up  the NaNs    default
  % iFix = -1;   %% live with the NaNs

  iFix = -1   %% live with the NaNs
   
  if iFix > 0
    disp('removing the NANS')
    for i = 1:length(prof.stemp)
      prof.nlevs(1,i) = sum(~isnan(prof.ptemp(:,i)));
      prof.gas_1(prof.nlevs(1,i):end,i) = prof.gas_1(prof.nlevs(1,i),i);
      prof.gas_3(prof.nlevs(1,i):end,i) = prof.gas_3(prof.nlevs(1,i),i);
      prof.ciwc(prof.nlevs(1,i):end,i)  = prof.ciwc(prof.nlevs(1,i),i);
      prof.clwc(prof.nlevs(1,i):end,i)  = prof.clwc(prof.nlevs(1,i),i);

      %this is to dangerous, will make nlevs === 42, since we reset nlevs a couple of lines later
      %prof.ptemp(prof.nlevs(1,i):end,i) = interp1([prof.plevs(prof.nlevs(1,i),i) prof.spres(1,i)],...
      %                                            [prof.ptemp(prof.nlevs(1,i),i) prof.stemp(1,i)],...
      %                                            prof.plevs(prof.nlevs(1,i):end,i),[],'extrap');

      %much safer and better, wrt nlevs
      prof.ptemp(prof.nlevs(1,i):end,i) = interp1([prof.plevs(prof.nlevs(1,i),i) prof.spres(1,i)],...
                                                  [prof.ptemp(prof.nlevs(1,i),i) prof.stemp(1,i)],...
                                                  prof.plevs(prof.nlevs(1,i):end,i));
      prof.nlevs(1,i) = sum(~isnan(prof.ptemp(:,i)));
    end
  else
    disp('NOT!!!!!!!!!! removing the NANS')
  end

end
