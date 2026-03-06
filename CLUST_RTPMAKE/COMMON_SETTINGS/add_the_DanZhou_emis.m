%% default is to use rtp_add_emis for : Masuda over ocean, Dan Zhou emis over land

% [h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
% p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
% p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
% [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);

%% hopefully you have done this already
%% p.rlon = wrapTo180(p.rlon);

if exist([set_path_to_danz '/danz_interpolant.mat'])
  [p,pa] = rtp_add_emis(p,pa);
else  
  disp('no /asl/data/iremis/danz/danz_interpolant.mat so skip rtp_add_emis  ... fix paths in set_path_to_danz.m. ... for now use constant emis')
  disp('no /asl/data/iremis/danz/danz_interpolant.mat so skip rtp_add_emis  ... fix paths in set_path_to_danz.m. ... for now use constant emis')
  disp('no /asl/data/iremis/danz/danz_interpolant.mat so skip rtp_add_emis  ... fix paths in set_path_to_danz.m. ... for now use constant emis')
  error('my codes assumes DanZhou emissivity so bit dangerous to set constant 0.98 emis without the user being told')
  p.nemis = ones(size(p.stemp)) * 2;
  p.efreq = [600 3000]' * ones(1,length(p.stemp));
  p.emis  = [0.98 0.98]' * ones(1,length(p.stemp));
  p.rho = (1-p.emis)/pi;
end
