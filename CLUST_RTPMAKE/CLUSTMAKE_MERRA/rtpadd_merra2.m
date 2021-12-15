function [head hattr prof pattr] = rtpadd_merra(head,hattr,prof,pattr,fields)
% function [head hattr prof pattr] = rtpadd_merra(head,hattr,prof,pattr,fields)
%
% root - usuall '/asl' (optional)
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
  %
  % Data is given on 3hr files, except 'ts' which is hourly
  %


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

    ptemp=[]; pgas_1=[]; pgas_3=[]; pciwc = []; pclwc = [];


    %%%%%%%%%%%%%%%%%%%% 
    % Interpolate 3D variables for each layer
 
    % t  - air temperature  	ptemp
    [dat_t plevs lats lons]= getdata_merra2(reqtime, 't',[]);
    msize = size(dat_t);
    dat_t = permute(dat_t,[3 1 2]);


    %%%%%%%%%%%%%%%%%%%%
    % ptime is computed by reverting mattime to TAI using the reference year above
    reqTAI = mattime2tai(reqtime, mrefyear);
    tprof.ptime = ones(1,nfovs).*reqTAI; % [sec]


    %%%%%%%%%%%%%%%%%%%%
    % Field quantities. Use Nearest Match. 

    isphere = single(interp_sphere(lats, lons, reshape(1:prod(msize(1:2)),msize(1:2)), tprof.rlat, tprof.rlon, 'nearest'));
    [ix, iy] = ind2sub(msize(1:2),isphere);

    tprof.plat = (lats(iy))';
    tprof.plon = (lons(ix))';


    nlevs=numel(plevs);


    parfor ilev=1:nlevs
      ptemp(ilev,:) = single(dat_t(nlevs-ilev+1,isphere));
    end
    
    % qv - specific humidity	gas_1
    [dat_q]= permute(getdata_merra2(reqtime, 'qv',[]),[3 1 2]);
    parfor ilev=1:nlevs
      pgas_1(ilev,:) = single(dat_q(nlevs-ilev+1,isphere));
    end
   
    % o3 - ozone mixing ratio	gas_3
    [dat_o3]= permute(getdata_merra2(reqtime, 'o3',[]),[3 1 2]);
    parfor ilev=1:nlevs
      pgas_3(ilev,:) = single(dat_o3(nlevs-ilev+1,isphere));
    end

    % ciwc - ice cloud
    [dat_q]= permute(getdata_merra2(reqtime, 'qi',[]),[3 1 2]);
    parfor ilev=1:nlevs
      pciwc(ilev,:) = single(dat_q(nlevs-ilev+1,isphere));
    end

    % clwc - liquid cloud
    [dat_q]= permute(getdata_merra2(reqtime, 'ql',[]),[3 1 2]);
    parfor ilev=1:nlevs
      pclwc(ilev,:) = single(dat_q(nlevs-ilev+1,isphere));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Compute Valid Levels
   
    % For each profile, find the last valid level. 
    % 
    tprof.plevs = single(flipud(plevs)*ones(1, nfovs));
    tprof.nlevs = single(nlevs.*zeros(1, nfovs));
    tprof.nlevs(1,:) = single(sum(~isnan(ptemp)));

    tprof.ptemp = single(ptemp);
    tprof.gas_1 = pgas_1;
    tprof.gas_3 = pgas_3;
    tprof.ciwc  = pciwc;
    tprof.clwc  = pclwc;

    % ps - surface pressure	spres
    [dat_ps xx lats lons]= getdata_merra2(reqtime, 'ps',[]); % It is in Pa, convert to mbar -> /100
    tprof.spres = single(interp_sphere(lats,lons,dat_ps/100, tprof.rlat, tprof.rlon, 'nearest'));
    
    % ts - surface temperature	stemp
    [dat_ts]= getdata_merra2(reqtime, 'ts',[]);
    tprof.stemp = single(interp_sphere(lats,lons,dat_ts, tprof.rlat, tprof.rlon, 'nearest'));
    
    % wind speed at 2m
    [dat_u2m]= getdata_merra2(reqtime, 'u2m',[]);
    [dat_v2m]= getdata_merra2(reqtime, 'v2m',[]);
    dat_w2m = abs(dat_u2m + dat_v2m*1i);
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

  hattr = set_attr(hattr,'profile','MERRA OpenDap using NetCDF libs, Nearest match');

  prof=tprof;

end
