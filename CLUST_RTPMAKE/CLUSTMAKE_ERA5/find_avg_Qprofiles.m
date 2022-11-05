clear  pav pQav globalavg globalQavg

bt1231 = rad2bt(1231,p2.rcalc(1520,:));
for ii = 1 : length(quants)-1
  boo = find(bt1231 >= booQ(ii));
  fprintf(1,'ii = %2i length(boo) = %5i \n',ii,length(boo))
  pav.stemp(ii) = nanmean(pnew_op.stemp(boo));
  pav.ptemp(:,ii) = nanmean(pnew_op.ptemp(:,boo),2);
  pav.gas_1(:,ii) = nanmean(pnew_op.gas_1(:,boo),2);
  pav.gas_3(:,ii) = nanmean(pnew_op.gas_3(:,boo),2);
  pav.plevs(:,ii) = nanmean(pnew_op.plevs(:,boo),2);
  pav.cldcalc(ii) = nanmean(pnew_op.rcalc(1520,boo),2);
  pav.clrcalc(ii) = nanmean(pnew_op.sarta_rclearcalc(1520,boo),2);

  if ii == 1
    globalavg = find_average_rtp(hnew_op,pnew_op,1,boo);
  else
    junk = find_average_rtp(hnew_op,pnew_op,1,boo);
    [hnew_op,globalavg] = cat_rtp(hnew_op,globalavg,hnew_op,junk);
  end

  %%%%%%%%%
  zoo = find(bt1231 >= booQ(ii) & bt1231 < booQ(ii+1));
  fprintf(1,'ii = %2i length(zoo) = %5i \n',ii,length(zoo))
  pQav.stemp(ii) = nanmean(pnew_op.stemp(zoo));
  pQav.ptemp(:,ii) = nanmean(pnew_op.ptemp(:,zoo),2);
  pQav.gas_1(:,ii) = nanmean(pnew_op.gas_1(:,zoo),2);
  pQav.gas_3(:,ii) = nanmean(pnew_op.gas_3(:,zoo),2);
  pQav.plevs(:,ii) = nanmean(pnew_op.plevs(:,zoo),2); 
  pQav.cldcalc(ii) = nanmean(pnew_op.rcalc(1520,zoo),2);
  pQav.clrcalc(ii) = nanmean(pnew_op.sarta_rclearcalc(1520,zoo),2);

  if ii == 1
    globalQavg = find_average_rtp(hnew_op,pnew_op,1,zoo);
  else
    junk = find_average_rtp(hnew_op,pnew_op,1,zoo);
    [hnew_op,globalQavg] = cat_rtp(hnew_op,globalQavg,hnew_op,junk);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = length(quants);
  boo = find(bt1231 == max(bt1231),1);
  pav.stemp(ii) =   pnew_op.stemp(boo);
  pav.ptemp(:,ii) = pnew_op.ptemp(:,boo);
  pav.gas_1(:,ii) = pnew_op.gas_1(:,boo);
  pav.gas_3(:,ii) = pnew_op.gas_3(:,boo);
  pav.plevs(:,ii) = pnew_op.plevs(:,boo);
  pav.cldcalc(ii) = nanmean(pnew_op.rcalc(1520,boo),2);
  pav.clrcalc(ii) = nanmean(pnew_op.sarta_rclearcalc(1520,boo),2);

[~,junk] = subset_rtp_allcloudfields(hnew_op,pnew_op,[],[],boo);
if isfield(junk,'calflag')
  junk = rmfield(junk,'calflag');
end
if isfield(junk,'pnote')
  junk = rmfield(junk,'pnote');
end
if isfield(junk,'iudef')
  junk = rmfield(junk,'iudef');
end
[hnew_op,globalavg] = cat_rtp(hnew_op,globalavg,hnew_op,junk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
