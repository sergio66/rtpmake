numdone = zeros(16,28);

nmax = 388; %% done till Aug 2019
nmax = 411; %% done till Aug 2020

iDorA = +1;
for ddd = 1 : 16     %% each timestep has 16 days
  for eee = 1 : 28    %% each tile has 28 ERA grid points
    if iDorA > 0
      dirx = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eee,'%02d') '/'];
    else
      dirx = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/ASC/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eee,'%02d') '/'];
    end
    thedir = dir([dirx '/*mat']);
    numdone(ddd,eee) = length(thedir); % ideally should be 365/16 * 17 ~ 388 since we go 2002/09 to 2019/08
  end
end

jett = jet; jett(1,:) = 1;
imagesc(numdone); title(['TimeSteps done out of ' num2str(nmax)]); 
ylabel('day(1-16)'); xlabel('ERA grid pt (1-28)'); colormap(jett); colorbar
cx = caxis; 
if cx(2)-nmax > 10
  cx(1) = cx(1)-10; cx(2) = cx(2)+10; caxis(cx);
else
  colormap(jet)
end

dstr = datestr(datetime('now'));
junk = [16*28*nmax nmax sum(numdone(:)) sum(numdone(:))*100/(16*28*nmax)];

fprintf(1,' %s done %6i of 16*28*%3i=%7i = %8.6f \n',dstr,junk);

