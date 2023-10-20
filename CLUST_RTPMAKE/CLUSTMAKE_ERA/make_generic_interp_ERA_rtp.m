function [h,ha,p,pa] = make_generic_interp_ERA_rtp(hin,hain,pin,pain);

h  = hin;
ha = hain;
p  = pin;
pa = pain;

    %%% this is NEW
    p.landfrac_fromL1B = p.landfrac;
    p.salti_fromL1B = p.salti;
    [salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
    p.landfrac = landfrac;
    p.salti    = salti;

    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};

    %[h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm
    %[p,h] = sergio_fill_ecmwf(p,h,'/asl/data/ecmwf/',-1);
    %save test_2002_09_08_g044_sergio.mat p h

    addpath /home/sergio/MATLABCODE/TIME
    [xyy,xmm,xdd,xhh] = tai2utcSergio(p.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
    time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
    co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
    p.co2ppm = co2ppm;
    fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),p.co2ppm(1));
    fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),p.co2ppm(end));

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% [p,h] = fill_era_interp(p,h);                       %%% add on era the NEW WAY
    new_4_ERAfiles_interp_analysis                          %%% even better Oct 2022
    %%%%%%%%%%%%%%%%%%%%%%%%%

    p0 = p;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
    %addpath /asl/rtp_prod2/emis/
    %addpath /asl/rtp_prod2/util/    
    %addpath /asl/packages/rtp_prod2/emis/
    %addpath /asl/packages/rtp_prod2/util/
    %addpath /asl/rtp_prod2/emis/
    %addpath /asl/rtp_prod2/util/
 
    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

    p.rlon = wrapTo180(p.rlon);
    [p,pa] = rtp_add_emis(p,pa);

    %figure(1)
    %scatter_coast(p.rlon,p.rlat,10,p.nemis); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h = hB2;

    disp('now you run "[p2] = driver_sarta_cloud_rtp(hB2,ha,p,pa,run_sarta);" ');
    disp('now you run "[p2] = driver_sarta_cloud_rtp(hB2,ha,p,pa,run_sarta);" ');
    disp('now you run "[p2] = driver_sarta_cloud_rtp(hB2,ha,p,pa,run_sarta);" ');

