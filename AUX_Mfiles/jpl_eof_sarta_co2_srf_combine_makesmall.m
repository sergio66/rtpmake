xstartup

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/RTPMAKE/

cd /asl/data/rtprod_airs/2009/03/01/
for gg = 1 : 240
  fnameIN  = ['cloudy_airs_l1b_ecm_sarta_baum_ice.2009.03.01.' num2str(gg,'%03d') '.rtp'];
  fnameOUTMAT = ['SMALLFILES/xjpl_l1b2834_ecm_sarta_baum_ice_2009_03_01_' num2str(gg,'%03d') '.mat'];
  fnameOUTRTP = ['RTPSMALLFILES/xjpl_l1b2834_ecm_sarta_baum_ice_2009_03_01_' num2str(gg,'%03d') '.rtp'];
  ee = exist(fnameIN);

  iaFound = zeros(1,12150);

  if ee > 0
    [h,ha,p,pa] = rtpread(fnameIN);
    fprintf(1,'read in %s \n',fnameIN)

    px.rcalc0 = ones(2834,12150)*NaN;
    px.rcalc0(1:2378,:) = p.rcalc;      %% this is what was in my orig sim files

    px.rcalc  = ones(2834,12150)*NaN;
    px.rlon   = p.rlon;
    px.rlat   = p.rlat;
    px.satzen = p.satzen;
    px.solzen = p.solzen;
    px.scanang  = p.scanang;
    px.stemp    = p.stemp;
    px.landfrac = p.landfrac;

    rcalc0 = px.rcalc;

    %% now have to put the new calcs in, together with [co2][sarta]
    for ii = 1 : 3
      for jj = 1 : 3
        xyz  = ['JPL_EOF_' num2str(ii) num2str(jj)];
        fxyz = [xyz '/jpl_' num2str(ii) num2str(jj) '_' num2str(gg,'%03d') '.rtp'];
        [hxyz,ha,pxyz,pa] = rtpread(fxyz);
        iaInd = (pxyz.atrack-1)*90 + pxyz.xtrack;
        iaFound(iaInd) = 1;
        px.rcalc(:,iaInd) = pxyz.rcalc;
        px.co2_sarta_model(iaInd) = str2num([num2str(ii) num2str(jj)]);
      end
    end
    if sum(iaFound) ~= 12150
      error('oops : not all fovs found!!')
    end


%    for ii = 1 : length(p.stemp)
%      nemis = 1:p.emis(ii);
%      px.emis1231(ii) = interp1(p.efreq(1:nemis,ii).p.eemis(1:nemis,ii),1231.2);
%    end

%    saver = ['save -v7.3 ' fnameOUT ' px comment'];
%    fprintf(1,'saver = %s \n',saver)
%    eval(saver)

    rtpwrite(fnameOUTRTP,hxyz,ha,px,pa);

    figure(1); scatter_coast(px.rlon,px.rlat,20,px.co2_sarta_model); 
      title(['co2/sarta model for gran ' num2str(gg)])
    figure(2); plot(hxyz.vchan,rad2bt(hxyz.vchan,px.rcalc0) - rad2bt(hxyz.vchan,px.rcalc))
      title(['co2/sarta model for gran ' num2str(gg)])
    pause(0.1); 

  end
end
