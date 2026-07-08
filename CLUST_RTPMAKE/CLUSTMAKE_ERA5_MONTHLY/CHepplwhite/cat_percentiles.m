
   dnum: [23x72 double]
    lonB: [1x73 double]
    nsam: [23x72 single]
    pctr: [23x72x2645x100 single]
     qar: [23x72x2645x101 single]
    zonB: {1x64 cell}



x.home = '/home/chepplew/data/rates_anomalies/tiling/percentiles/';

x.dir = dir([x.home '20*_S02p75_percentiles.mat']);


dnum = [];
pctr = [];
for ifn = 1:length(x.dir)
  %y1 = load([x.dir(ifn).folder '/' x.dir(ifn).name]);
  y1 = load([x.dir(ifn).folder '/' x.dir(ifn).name],'dnum','pctr');
  dnum = cat(1, dnum, y1.dnum);
  pctr = cat(1, pctr, y1.pctr);

  fprintf(1,'.')
end
