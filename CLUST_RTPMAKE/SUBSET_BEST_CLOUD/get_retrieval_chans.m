function ind = get_retrieval_chans(h,g,iStemp_ColWV,iAIRSL2_or_Sergio);

if nargin == 3
  iAIRSL2_or_Sergio = +1;
end

if iAIRSL2_or_Sergio == -1
  disp('doing VERY SIMPLE choices for choosing channels')
  if iStemp_ColWV == 2
    ind = find((h.vchan >= 800 & h.vchan < 962) | (h.vchan >= 1100 & h.vchan < 1250));
  elseif iStemp_ColWV == 3
    ind = find((h.vchan >= 600 & h.vchan < 962) | (h.vchan >= 1100 & h.vchan < 1650));
  elseif iStemp_ColWV >= 4 & iStemp_ColWV <= 44
    ind = find(h.vchan >= 600 & h.vchan < 1650);
  elseif iStemp_ColWV >= 45
    ind = find(h.vchan < 2400);
  end

  ind = intersect(g,ind);

  Rx = instr_chans('airs',2,250); Rx = Rx(ind);
  lala = find(Rx < 1);
  ind = ind(lala);

  ind = ind(1:4:length(ind));
  fprintf(1,' sergio chanset is %4i chans long \n',length(ind));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iAIRSL2_or_Sergio == +1
  disp('doing AIRS L2 choices for choosing channels')
  
  fairs = instr_chans;
  [airsL2_chans_set,airsL2_IDs_set] = retrieval_channel_set();
  fieldnamesx = fieldnames(airsL2_IDs_set);
  fieldtypesx = [1 2 3 4 5 7 10 13];
  newchans = airsL2_IDs_set.neural_net;
  for nn = 2 : length(fieldtypesx)
    fprintf(1,'%2i %s \n',nn,fieldnamesx{fieldtypesx(nn)});
    str = ['x = airsL2_IDs_set.' fieldnamesx{fieldtypesx(nn)} ';'];
    eval(str);
    newchans = unique(union(newchans,x));
  end  
  newchans = newchans(fairs(newchans) <= 1700);
  ind = newchans;
  ind = intersect(g,ind);
  fprintf(1,' default AIRS L2 chanset is %4i chans long \n',length(ind));
  
elseif iAIRSL2_or_Sergio == +2
  disp('doing AIRS L2 choices for choosing channels, including more WV chans')
  
  fairs = instr_chans;
  [airsL2_chans_set,airsL2_IDs_set] = retrieval_channel_set();
  fieldnamesx = fieldnames(airsL2_IDs_set);
  fieldtypesx = [1 2 3 4 5 7 10 13];
  newchans = airsL2_IDs_set.neural_net;
  for nn = 2 : length(fieldtypesx)
    fprintf(1,'%2i %s \n',nn,fieldnamesx{fieldtypesx(nn)});
    str = ['x = airsL2_IDs_set.' fieldnamesx{fieldtypesx(nn)} ';'];
    eval(str);
    newchans = unique(union(newchans,x));
  end  
  newchans = newchans(fairs(newchans) <= 1700);
  ind = newchans;
  ind = intersect(g,ind);
  fprintf(1,' default AIRS L2 chanset is %4i chans long \n',length(ind));
  
  wvchans = find(fairs >= 1370 & fairs <= 1700);
  wvchans = intersect(g,wvchans);
  wvchans = wvchans(1:3:length(wvchans));
  ind = union(ind,wvchans);
  fprintf(1,' with more WV chans added on, AIRS L2 chanset is %4i chans long \n',length(ind));  
end
