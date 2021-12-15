% rats: Math_tsfit_lin_robust requires daily data

addpath /asl/matlib/aslutil
addpath /home/strow/Matlab/Math


fn='/home/chepplew/data/atmosSci/UEA_CRU/HadCRUT4-gl.dat';
A=importdata(fn);
B=A(1:2:end,2:13);
YR = A(1:2:end,1);

xt1=datetime(1850,1,1)
xt2=datetime(2019,12,1)

nsam      = length(C(:))
numyr     = YR(end) - YR(1) + 1
numdays   = days(xt2 - xt1)
monthlist = 1:numyr*12;

X1 = datetime(YR(1), monthlist, 1);
C=B';

plot(X,C(:),'-')

junk = NaN(63240/31,1);
Cpad = [C(:) junk];

% Interpolate to daily grid
vq = interp1([1:31:63240],C(:),[1:numdays]);

[Bb, Bstats] = Math_tsfit_lin_robust([1:numdays],vq,2);
Ys_bar       = Math_timeseries_2([1:numdays],Bb);

