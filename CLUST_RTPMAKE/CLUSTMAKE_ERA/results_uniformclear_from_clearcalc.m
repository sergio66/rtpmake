rlat = [];
rlon = [];
solzen = [];
landfrac = [];
gran = [];

for JOB = 1 : 240
  gstr = num2str(JOB,'%03d');
  fin = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/09/06/clear_airs_l1c_era.2018.09.06.' gstr '.rtp'];
  fout = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/09/06/uniform_airs_l1c_era.2018.09.06.' gstr '.rtp'];
  if exist(fout)
    [head,hattr,prof,pattr] = rtpread(fout);
    rlat = [rlat prof.rlat];
    rlon = [rlon prof.rlon];
    solzen = [solzen prof.solzen];
    landfrac = [landfrac prof.landfrac];    
    gran = [gran ones(size(prof.rlat))*JOB];
  end
end

steve = '/asl/rtp/rtp_airicrad_v6/clear/2018/era_airicrad_day249_clear.rtp';
[hs,has,psteve,pas] = rtpread(steve);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oo = find(landfrac == 0);
figure(1); scatter_coast(rlon(oo),rlat(oo),5,ones(size(oo))); title('sergio')

oo = find(psteve.landfrac == 0);
figure(2); scatter_coast(psteve.rlon(oo),psteve.rlat(oo),5,ones(size(oo))); title('steve')
