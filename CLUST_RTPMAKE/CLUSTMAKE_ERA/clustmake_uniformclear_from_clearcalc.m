JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 55  %% most daytime clear, according to Steve
%JOB = 103 %% wierdnees near S. Africa/Mozambique

addpath /asl/matlab2012/airs/readers
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/rtp_prod2/airs/util

gstr = num2str(JOB,'%03d');
fin = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/09/06/clear_airs_l1c_era.2018.09.06.' gstr '.rtp'];
fout = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/09/06/uniform_airs_l1c_era.2018.09.06.' gstr '.rtp'];

figure(1);
if exist(fin)
  [head,hattr,prof,pattr] = rtpread(fin);

  [iclear, dbtun] = airs_uniform_clear(head, hattr, prof, pattr);
  i1231 = find(head.vchan >= 1231,1);
  scatter_coast(prof.rlon,prof.rlat,10,rad2bt(1231,prof.robs1(i1231,:))); colormap jet
  hold on; plot(prof.rlon(iclear),prof.rlat(iclear),'x','linewidth',5,'markersize',5,'color','k'); 
  ax = axis;
  hold off

  whos iclear

  [head,prof] = subset_rtp(head,prof,[],[],iclear);
  rtpwrite(fout,head,hattr,prof,pattr)
end

iDoSteve = -1;
if iDoSteve > 0
  dinm = [31 28 31 30 31 30 31 31 30 31 30 31];
  doy = sum(dinm(1:8)) + 6

  steve = '/asl/rtp/rtp_airicrad_v6/clear/2018/era_airicrad_day249_clear.rtp';
  if ~exist('psteve')
    [hs,has,psteve,pas] = rtpread(steve);
    for ix = 1 : 240
      woo = find(psteve.findex == ix);
      if length(woo) > 0
        steve_uc.number(ix) = length(woo);
        steve_uc.solzen(ix) = mean(psteve.solzen(woo));
      else
        steve_uc.number(ix) = 0;
        steve_uc.solzen(ix) = NaN;
      end
    end
  end
  figure(2); ix = find(steve_uc.solzen < 90); plot(ix,steve_uc.number(ix),'x-')
  
  ix = find(psteve.findex == JOB);
  if length(ix) > 0
    whos ix
    xi1231 = find(hs.vchan >= 1231,1);
    figure(2); 
    scatter_coast(psteve.rlon(ix),psteve.rlat(ix),10,rad2bt(xi1231,psteve.robs1(1291,ix))); colormap jet
    axis(ax)
  end
end

