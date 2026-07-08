function yavg = average_Nprofiles(iN,x1,x2,x3,x4,x5,x6,x7,x8)

if nargin ~= (iN + 1)
  fprintf(1,'iN == %2i so expect those many profiles\n',iN)
end

names = fieldnames(x1);
yavg = x1;

for ii = 2 : iN
  junker = ['junk = x' num2str(ii) ';']; eval(junker);
  for ff = 1 : length(names)
    moo = names{ff};
    junker = ['yavg.' moo ' = yavg.' moo '+ junk.' moo ';']; 
    eval(junker); 
    %fprintf(1,'A %s \n',junker)
  end
end

for ff = 1 : length(names)
  moo = names{ff};
  junker = ['yavg.' moo ' = yavg.' moo '/' num2str(iN) ';']; 
  eval(junker); 
  %fprintf(1,'B %s \n',junker)
end
