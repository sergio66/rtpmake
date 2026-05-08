% if iSNPPorJ1orJ2 == 0
%   NONONOdout = ['/asl/rtp/cris/npp_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
%   dout = ['/asl/s1/sergio/rtp/npp_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
% elseif iSNPPorJ1orJ2 == 1
%   NONONOdout = ['/asl/rtp/cris/j1_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
%   dout = ['/asl/s1/sergio/rtp/j1_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
% else
%   error('unknow SNPP, J1 or ... ?')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iSNPPorJ1orJ2 == 0
  % NONONOdout = ['/asl/rtp/cris/npp_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
  dout = ['/asl/s1/sergio/rtp/npp_ccast_hires/'];
  dout = ['/umbc/rs/pi_sergio/WorkDirDec2025/sergio_temp_rtp_files/npp_ccast_hires/'];
  dout = ['/home/sergio/nogit/sergio_temp_rtp_files/npp_ccast_hires/'];    
  dout = [dout  '/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
elseif iSNPPorJ1orJ2 == 1
  % NONONOdout = ['/asl/rtp/cris/j1_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
  dout = ['/asl/s1/sergio/rtp/j1_ccast_hires/'];
  dout = ['/umbc/rs/pi_sergio/WorkDirDec2025/sergio_temp_rtp_files/j1_ccast_hires/'];
  dout = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/'];
  dout = [dout '/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
else
  error('unknow SNPP, J1 or ... ?')
end
