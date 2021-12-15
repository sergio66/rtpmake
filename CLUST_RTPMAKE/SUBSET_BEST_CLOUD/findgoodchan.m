%% see /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
function g = findgoodchan(h,ha,p,pa);

clrscenes = 1 : length(p.stemp);

addpath /home/sergio/MATLABCODE/CLOUD

g = dogoodchan;

lessthan0 = zeros(1,2378);

for ii = 1 : 2378
  lala = p.robs1(ii,:);
  oo = find(lala < 0);
  if length(oo) > 0
    lessthan0(ii) = 1;
  end
end
greaterthan0 = ones(1,2378);
oo = find(lessthan0 == 1);
fprintf(1,'found %4i channels where at least one rad < 0 \n',length(oo));
greaterthan0(oo) = 0;

fprintf(1,' >>> loop 0 : len(g) = %4i cldprofs = %4i \n',length(g),length(clrscenes))

figure(1); clf; plot(g,p.calflag(g,:)); title('calflag'); pause(0.1);
figure(2); clf; imagesc(g,1:length(p.stemp),double(p.calflag(g,:)')); 
  shading flat; colorbar; xlabel('g'); ylabel('profile num')
g0 = g;

calflag = p.calflag;
woo = find(p.calflag > 64);  calflag(woo) = -1;
woo = find(p.calflag <= 64); calflag(woo) = 1;
figure(3); clf; imagesc(1:2378,1:length(p.stemp),calflag);
  shading flat; colorbar; xlabel('chan num'); ylabel('profile num')

kk = 0;
iYes = +1;
while iYes > 0 & kk < 5
  kk = kk+1;
  %% want to remove all bad chans!!!!
  %% keep only channs which are good quality for > 95% of the profiles
  for ix = 1 : 3
    clear len_notbad
    for ii = 1 : length(g)
      notbad = find(p.calflag(g(ii),:) < 64);   %% see MATLABCODE/RATES_CLOUD/make_stats_clear_backwards.m
      len_notbad(ii) = length(notbad);
    end
    accept = find(len_notbad >= 0.90*length(p.stemp));
    g = g(accept);
  end
  figure(3); clf; plot(len_notbad(accept)); title('after checking chans ...')

  %% want to remove all bad profiles!!!! find profs where chans are not all good
  % for ii = 1 : length(p.stemp)
  %   notbad = find(p.calflag(g,ii) < 64);   %% see MATLABCODE/RATES_CLOUD/make_stats_clear_backwards.m
  %   len_notbad2(ii) = length(notbad);
  % end
  % accept = find(len_notbad2 == length(g));
  % clrscenes = intersect(clrscenes,accept); %%% <<<<<<<<<<<<<<<<<<<

  figure(4); clf; imagesc(g,clrscenes,double(p.calflag(g,clrscenes)')); 
  shading flat; colorbar; xlabel('g'); ylabel('profile num'); title('Qual FLag FINAL')

  ponk = double(p.calflag(g,clrscenes));
  ponk = ponk(:);
  if max(ponk) >= 64
    iYes = +1;
  else
    iYes = -1;
  end

  fprintf(1,' >>> loop %2i len(g) = %4i cldprofs = %4i \n',kk,length(g),length(clrscenes))
  %disp('cleaning up chans/profs ... ret')
  %pause
end

greaterthan0 = find(greaterthan0 == 1);
g = intersect(g,greaterthan0);
fprintf(1,'final tally of "good" channels = %4i \n',length(g));
plot(1:2378,mean(p.robs1'),g,mean(p.robs1(g,:)'))
