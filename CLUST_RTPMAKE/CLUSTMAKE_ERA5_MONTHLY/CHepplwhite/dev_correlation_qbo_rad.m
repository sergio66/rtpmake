
% first: load_tile_daily_data()
% second: prep.: using stats_tile_daily_ssw()

iie = find(rdays);


xt1 = datetime(2003,01,01);
xt2 = datetime(2014,12,31);
iid = find(qbo.dtim >= xt1 & qbo.dtim <= xt2);

p70_mon = qbo.p70(iid);
p50_mon = qbo.p50(iid);
p30_mon = qbo.p30(iid);
p10_mon = qbo.p10(iid);

indays  = datenum(xt1):1:datenum(xt2);
p70_day = interp1(datenum(qbo.dtim(iid)), p70_mon, indays,'spline','extrap');
p50_day = interp1(datenum(qbo.dtim(iid)), p50_mon, indays,'spline','extrap');
p30_day = interp1(datenum(qbo.dtim(iid)), p30_mon, indays,'spline','extrap');
p10_day = interp1(datenum(qbo.dtim(iid)), p10_mon, indays,'spline','extrap');

figure;plot(datetime(rdays(iie),'convertfrom','datenum'), btm(ich,:),'-')

figure;plot(qbo.dtim(iid), qbo.p50(iid),'.-')
hold on; plot(datetime(indays,'convertfrom','datenum'), p50_day,'.-')
     
[~,iin] = intersect(indays, rdays(iie));
[~,iim] = intersect(rdays(iie), indays);

clf;plot(datetime(indays(iin),'convertfrom','datenum'), p50_day(iin),'.-')
hold on;
 plot(datetime(rdays(iie(iim)),'convertfrom','datenum'), 60*(btm(ich,iim)-nanmean(btm(ich,iim),2)),'-')

ich=7; 
X = p10_day(iin);
Y = 60*(btm(ich,iim)-nanmean(btm(ich,iim),2));
inan=find(isnan(Y));
Y(inan)=[];
X(inan)=[];
yhat = nanmean(Y);
xhat = nanmean(X);
bb = (X - xhat)*(Y' - yhat)/sum((X - xhat).^2);
R = corrcoef(X,Y); 
