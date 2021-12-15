indexless = [];
indexmore = [];
for ii = 1 : 97
  boo = find(papa.ptemp(ii,:) < (usstd(ii,4)-50) & papa.nlevs-1 >= ii);
  lessthan{ii} = boo;
  lessthancnt(ii) = length(boo);
  if length(boo) > 0
    lessthanmagnitude{ii} = papa.ptemp(ii,boo) - (usstd(ii,4)-50);
    maxlessthanmagnitude(ii) = max(abs(papa.ptemp(ii,boo) - (usstd(ii,4)-50)));
    meanlessthanmagnitude(ii) = mean(abs(papa.ptemp(ii,boo) - (usstd(ii,4)-50)));
    stdlessthanmagnitude(ii) = std(abs(papa.ptemp(ii,boo) - (usstd(ii,4)-50)));
  else
    lessthanmagnitude{ii} = [];
    maxlessthanmagnitude(ii) = 0.0;
    meanlessthanmagnitude(ii) = 0.0;
    stdlessthanmagnitude(ii) = 0.0;
  end
  if length(indexless) == 0 & length(boo) > 0
    indexless = boo;
  elseif length(indexless) > 0 & length(boo) > 0
    indexless = union(indexless,boo);
  end

  boo = find(papa.ptemp(ii,:) > (usstd(ii,4)+50) & papa.nlevs-1 >= ii);
  morethan{ii} = boo;
  morethancnt(ii) = length(boo);
  if length(boo) > 0
    morethanmagnitude{ii} = papa.ptemp(ii,boo) - (usstd(ii,4)-50);
    maxmorethanmagnitude(ii) = max(abs(papa.ptemp(ii,boo) - (usstd(ii,4)-50)));
  else
    morethanmagnitude{ii} = [];
    maxmorethanmagnitude(ii) = 0.0;
  end
  if length(indexmore) == 0 & length(boo) > 0
    indexmore = boo;
  elseif length(indexmore) > 0 & length(boo) > 0
    indexmore = union(indexmore,boo);
  end
end

figure(2); 
plot(lessthancnt,usstd(1:97,2)*1013.25,morethancnt,usstd(1:97,2)*1013.25); set(gca,'ydir','reverse')
length(indexless)/length(papa.stemp)

addpath /home/sergio/MATLABCODE/PLOTMISC
figure(3); 
plot(maxlessthanmagnitude,usstd(1:97,2)*1013.25,maxmorethanmagnitude,usstd(1:97,2)*1013.25); set(gca,'ydir','reverse'); hold on
shadedErrorBarY(meanlessthanmagnitude,usstd(1:97,2)*1013.25,stdlessthanmagnitude,'b',0.2);
xlabel('How far below USStd-50 K'); ylabel('p(mb)')
fprintf(1,'%6i %6i %8.6f \n',[length(indexless) length(papa.stemp) length(indexless)/length(papa.stemp)])

figure(4); 
plot(papa.stemp(indexless),papa.spres(indexless),'.'); ylabel('spres'); xlabel('stemp')
set(gca,'ydir','reverse')

figure(5); 
semilogy(papa.ptemp(:,indexless),papa.plevs(:,indexless),'b',...
         mean(papa.ptemp(:,indexless)'),mean(papa.plevs(:,indexless)'),'c',...
         papa.ptemp(:,indexless(1)),papa.plevs(:,indexless(1)),'c',...
         usstd(:,4),1013.25*usstd(:,2),'r.-',...
         usstd(:,4)-50,1013.25*usstd(:,2),'r',usstd(:,4)+50,1013.25*usstd(:,2),'r') 
set(gca,'ydir','reverse'); axis([180 280 0 1000])
axis([180 260 100 1000])

figure(6); plot(papa.rlon(indexless),papa.rlat(indexless),'.')

allthelessdeltas = [];
for ii = 1 : 97
  allthelessdeltas = [allthelessdeltas lessthanmagnitude{ii}];
end
figure(7); delta = 0 : 0.25 : 20; plot(delta,histc(abs(allthelessdeltas),delta)); xlabel('Distance from TStd-50'); ylabel('hist')
figure(7); delta = 0 : 0.50 : 20; plot(delta,histc(abs(allthelessdeltas),delta)); xlabel('Distance from TStd-50'); ylabel('hist')
figure(7); delta = 0 : 1.00 : 20; plot(delta,histc(abs(allthelessdeltas),delta)); xlabel('Distance from TStd-50'); ylabel('hist')
grid on
