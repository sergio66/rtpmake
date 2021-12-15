function [h,ha,p,pa] = driver_make_rtp_airsL3climatology(yIN,mIN,dIN,gIN);

%% depends on climatology in /asl/s1/sergio/AIRS_L3/airs_L3v6_*_Sept2014.mat

addpath /home/sergio/MATLABCODE/AIRS_L3/
[h,ha,p,pa] = rtp_airsL3_climatology(yIN,mIN,dIN,gIN);

