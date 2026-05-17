% fdirOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
% fdirOUT = ['/asl/rtp/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
if iv5or6 == 5
  fdirOUT = ['/asl/s1/sergio/rtp/rtp_airibrad_v5/' ystr '/' mstr '/' dstr '/'];  %% till 2018
  fdirOUT = ['/umbc/rs/pi_sergio/WorkDirDec2025/sergio_temp_rtp_files/rtp_airibrad_v5/' ystr '/' mstr '/' dstr '/'];  %% till 2018
  fdirOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/rtp_airicrad_v5/' ystr '/' mstr '/' dstr '/'];                 %% after 2018
elseif iv5or6 == 6
  fdirOUT = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/' ystr '/' mstr '/' dstr '/'];  %% after 2018
  fdirOUT = ['/umbc/rs/pi_sergio/WorkDirDec2025/sergio_temp_rtp_files/rtp_airibrad_v6/' ystr '/' mstr '/' dstr '/'];  %% till 2018
  fdirOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/rtp_airicrad_v6/' ystr '/' mstr '/' dstr '/'];                 %% after 2018  
end
