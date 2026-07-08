scalar_stuff = load('/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/scalar_2002_09_2022_08.mat');
profile_stuff = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_asc.mat');

st_profile = profile_stuff.all.stemp;
st_scalar  = scalar_stuff.sktsaveA;
%% st_scalar  = scalar_stuff.sktsaveD;

figure(1); clf; pcolor(reshape(nanmean(st_profile - st_scalar,1),72,64)'); shading interp; colormap(usa2); caxis([-1 +1]*2); colorbar
figure(2); clf; ii = 1; pcolor(reshape(st_profile(ii,:) - st_scalar(ii,:),72,64)'); shading interp; colormap(usa2); caxis([-1 +1]*2); colorbar

