function [head hattr prof pattr] = rtpadd_emis_DanZhou(head,hattr,prof,pattr)

%function [head hattr prof pattr] = rtpadd_emis_DanZhou(head,hattr,prof,pattr)
%
% This function takes an rtp profile / file and adds/replaces the emissivity data
%
% The idea behind this function is to add all the emissivity needed fields to a profile structure so 
%   that higher order functions may be able to interchangably call the emissivity of different methods

%  Written 17 March 2011 - Paul Schou


addpath /asl/matlib/science        % emissivity function
addpath /asl/matlib/aslutil        % get_attr function

debug = 1;  %indicate that we are in debug an print out a bunch of checks
debug = 0;

%%%%%%  Begin rtpadd_emis_DanZhou

    % convert to matlab time
    rtime_str = get_attr(pattr,'rtime');
    if isempty(rtime_str); rtime_str = get_attr(pattr,'L1bCM rtime'); end  % backwards compatibility
    if debug; disp(['rtime str = ' rtime_str]); end

    st_year = 1993;  % default to 1993 year
    if length(rtime_str) > 4
      st_year = str2num(rtime_str(end-4:end));
      if st_year < 1993
        disp('Warning [rtpadd_ecmwf]: The rtime in pattr is in an invalid format')
        st_year = 1993;
      end
    else
      disp('Warning [rtpadd_ecmwf]: Could not find rtime in pattr for start year')
    end
    npro = length(prof.rtime);

    % select and generate the land emissivity set
    iland = prof.landfrac > 0.01;
    iLandOK = sum(iland);
    if iLandOK > 0
      [land_efreq, land_emis]=emis_DanZhou(prof.rlat(iland),prof.rlon(iland),prof.rtime(iland),st_year);
      land_nemis = length(land_efreq);
      iland(iland) = max(land_emis) >= 0;
      land_emis = land_emis(:,max(land_emis) >= 0);
      end
    iLandOK = sum(iland);

    % select and generate the sea emissivity set
    isea = prof.landfrac < 0.99 | ~iland;
    iSeaOK = sum(isea);
    if iSeaOK > 0
      [sea_nemis, sea_efreq, sea_emis]=cal_seaemis2(prof.satzen(isea),prof.wspeed(isea));
      sea_efreq = sea_efreq(:,1);
      sea_nemis = length(sea_efreq);
      end
  
    % interpolation frequency set
    if iLandOK > 0 & iSeaOK > 0
      interp_efreq = sort([land_efreq]);
      interp_efreq = interp_efreq(interp_efreq >= max(min(land_efreq),min(sea_efreq)) & ...
                                interp_efreq <= min(max(land_efreq),max(sea_efreq)));
      interp_nemis = length(interp_efreq);

      % predeclare the size of arrays
      prof.emis = zeros([max(land_nemis,sea_nemis), npro],'single');
      prof.efreq = ones([max(land_nemis,sea_nemis), npro],'single') * -9999;
      prof.nemis = zeros([1, npro],'uint8');

      % replace the clean stuff, here we put all the sea fovs into the sea emis
      prof.emis(1:sea_nemis,isea & ~iland) = sea_emis(:,~iland(isea));
      prof.efreq(1:sea_nemis,isea & ~iland) = repmat(sea_efreq(:),[1 sum(isea & ~iland)]);
      prof.nemis(1,isea & ~iland) = sea_nemis;

      % replace the land fovs in the emissivity
      prof.emis(1:land_nemis,~isea & iland) = land_emis(:,~isea(iland));
      prof.efreq(1:land_nemis,~isea & iland) = repmat(land_efreq(:),[1 sum(~isea & iland)]);
      prof.nemis(1,~isea & iland) = land_nemis;

      % interpolate the rest
      lf = prof.landfrac(isea & iland); lf(lf < 0) = 0; % landfrac
      prof.emis(1:interp_nemis,isea & iland) = ...
        bsxfun(@times,interp1(land_efreq,land_emis(:,isea(iland)),interp_efreq,'linear'), lf) + ...
        bsxfun(@times,interp1(sea_efreq,sea_emis(:,iland(isea)),interp_efreq,'linear'), (1-lf));
      prof.efreq(1:interp_nemis,isea & iland) = repmat(interp_efreq(:),[1 sum(isea & iland)]);
      prof.nemis(1,isea & iland) = interp_nemis;

      set_attr(pattr,'emis',['land(' which('emis_DanZhou') ')  water(' which('cal_seaemis2') ')']);

   elseif iLandOK > 0 & iSeaOK == 0
     prof.efreq(1:length(land_efreq),iland) = repmat(land_efreq,[1 sum(iland)]);
     prof.emis = land_emis;
     prof.nemis(1,iland) = length(land_efreq);

   elseif iLandOK == 0 & iSeaOK > 0
     prof.efreq(1:length(sea_efreq),isea) = repmat(sea_efreq,[1 sum(isea)]);
     prof.emis = sea_emis;
     prof.nemis(1,isea) = length(sea_efreq);
     end

end % Function end
