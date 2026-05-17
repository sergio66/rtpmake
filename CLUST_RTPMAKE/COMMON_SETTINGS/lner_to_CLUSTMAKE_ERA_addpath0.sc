# lrwxrwxrwx 1 sergio pi_sergio    35 Mar  5 06:16 set_path_to_danz.m -> ../COMMON_SETTINGS/set_path_to_danz.m
# lrwxrwxrwx 1 sergio pi_sergio    27 Mar  5 06:16 addpath0.m -> ../COMMON_SETTINGS/addpath0.m
# lrwxrwxrwx 1 sergio pi_sergio    56 Mar  5 06:15 add_the_paths_and_klayers_sarta_execs.m -> ../COMMON_SETTINGS/add_the_paths_and_klayers_sarta_execs.m
# lrwxrwxrwx 1 sergio pi_sergio    36 Mar  5 06:15 set_path_to_execs.m -> ../COMMON_SETTINGS/set_path_to_execs.m
# lrwxrwxrwx 1 sergio pi_sergio    31 Jan 25 05:10 set_filelist.m -> ../COMMON_SETTINGS/set_filelist.m
# lrwxrwxrwx 1 sergio pi_sergio    52 Jan 25 05:10 test_tcc.m -> /home/sergio/MATLABCODE/matlib/clouds/TCC/test_tcc.m
# lrwxrwxrwx 1 sergio pi_sergio    37 Jan 25 05:10 read_in_L1B_or_L1C.m -> ../COMMON_SETTINGS/read_in_L1B_or_L1C.m
# lrwxrwxrwx 1 sergio pi_sergio    28 Jan 25 05:10 onetime.txt -> ../COMMON_SETTINGS/onetime.txt
# lrwxrwxrwx 1 sergio pi_sergio    39 Mar  5 06:43 add_the_DanZhou_emis.m -> ../COMMON_SETTINGS/add_the_DanZhou_emis.m
# lrwxrwxrwx 1 sergio pi_sergio    40 Mar  5 10:27 sarta_chans_for_l1c.mat -> ../COMMON_SETTINGS/sarta_chans_for_l1c.mat
# lrwxrwxrwx 1 sergio pi_sergio    53 Mar  5 16:40 set_landfrac_using_L1B_L1C_or_usgs.m -> ../COMMON_SETTINGS/set_landfrac_using_L1B_L1C_or_usgs.m
# test_tcc.m --> /home/sergio/git/matlabcode/matlibSergio/matlib/clouds/TCC/test_tcc.m

########################################################################

echo "these are initial files, make sure they are SYMBOLIC LINKS and NOT actual files"
ls -lt set_path_to_danz.m addpath0.m add_the_paths_and_klayers_sarta_execs.m set_path_to_execs.m set_filelist.m \
      set_run_sarta_options.m  set_landfrac_using_L1B_L1C_or_usgs.m test_tcc.m read_in_L1B_or_L1C.m onetime.txt \
      add_the_DanZhou_emis.m sarta_chans_for_l1c.mat \
      plot_richardson_PBLH.m get_richardson_number_levels.m get_richardson_number_layers.m
read -p "Press [Enter] key to continue if they do not exist or they are SYMBOLIC LINKS ... (Ctrl C) if they are actual files"

echo "now setting symbolic links to files in COMMON_SETTINGS"

rm set_path_to_danz.m addpath0.m add_the_paths_and_klayers_sarta_execs.m set_path_to_execs.m set_filelist.m \
   set_run_sarta_options.m set_landfrac_using_L1B_L1C_or_usgs.m test_tcc.m read_in_L1B_or_L1C.m onetime.txt \
   add_the_DanZhou_emis.m sarta_chans_for_l1c.mat

rm get_richardson_number_levels.m get_richardson_number_layers.m plot_richardson_PBLH.m

#########################

ln -s ../COMMON_SETTINGS/onetime.txt                               .
ln -s ../COMMON_SETTINGS/read_in_L1B_or_L1C.m                      .
ln -s ../COMMON_SETTINGS/set_filelist.m                            . 

ln -s ../COMMON_SETTINGS/set_path_to_execs.m                       .
ln -s ../COMMON_SETTINGS/add_the_paths_and_klayers_sarta_execs.m   .
ln -s ../COMMON_SETTINGS/addpath0.m                                .
ln -s ../COMMON_SETTINGS/set_path_to_danz.m                        .
ln -s ../COMMON_SETTINGS/add_the_DanZhou_emis.m                    .
ln -s ../COMMON_SETTINGS/sarta_chans_for_l1c.mat                   .
ln -s ../COMMON_SETTINGS/set_landfrac_using_L1B_L1C_or_usgs.m      .
ln -s ../COMMON_SETTINGS/set_run_sarta_options.m                   .

ln -s ../COMMON_SETTINGS/get_richardson_number_levels.m            .
ln -s ../COMMON_SETTINGS/get_richardson_number_layers.m            .
ln -s ../COMMON_SETTINGS/plot_richardson_PBLH.m                    .

ln -s /home/sergio/git/matlabcode/matlibSergio/matlib/clouds/TCC/test_tcc.m . 

#########################

ls -lt set_path_to_danz.m addpath0.m add_the_paths_and_klayers_sarta_execs.m set_path_to_execs.m set_filelist.m \
      set_run_sarta_options.m  set_landfrac_using_L1B_L1C_or_usgs.m test_tcc.m read_in_L1B_or_L1C.m onetime.txt \
      add_the_DanZhou_emis.m sarta_chans_for_l1c.mat
      get_richardson_number_layers.m get_richardson_number_levels.m plot_richardson_PBLH.m
