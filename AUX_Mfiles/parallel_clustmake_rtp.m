matlabpool open
poolSize = matlabpool('size');
if poolSize == 0
   error('parallel:demo:poolClosed .... This demo needs an open MATLAB pool to run.');
end
fprintf('This demo is running on %d MATLABPOOL workers.\n', matlabpool('size'));


parfor JOB = 1 : 240
  clust_make_AIRSL2_90x135_cloudrtp_YYMMDD_loopGG
end


matlabpool close
