%% creates an rtp file for ONE granule
%% can be modified for more!

%% JOB 01 - 12 = 2002/01/01 to 2002/12/31 
%% JOB 13 - 24 = 2003/01/01 to 2003/12/31 
%% etc

% local running to test
% clustcmd -L clustcmd_loop_make_merrartp_sergio_clear_from_era.m 001:140
%
% otherwise when happy
% clustcmd -q medium -n 16 -p 2 clustcmd_loop_make_merrartp_sergio_clear_from_era.m 001:140
%
% or
% clustcmd -q long_contrib -n 16 -p 4 clustcmd_loop_make_merrartp_sergio_clear_from_era.m 001:140


%for JOB = 0 : 140

%% loops on a monthly basis

for JOBX = JOB
  yy = floor(JOBX/12);   %% assuming 31 days per month = 372 days per year

  mm = JOBX - yy*12 + 1;
  yy = yy + 2002;

  YY = yy;
  MM = mm;
  fprintf(1,'JOBX = %5i   YY = %4i MM = %2i \n',JOB,YY,MM)

  for DD = 1 : 31
    yymmddgg = [YY MM DD -1];
    loop_make_merrartp_sergio_clear_from_era  
  end
end


