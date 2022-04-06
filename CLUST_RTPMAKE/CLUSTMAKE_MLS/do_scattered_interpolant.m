for iix = boo(1):boo(end)
  junk = squeeze(aT(:,:,iix)); [goodx,goody] = find(isfinite(junk)); goodx = unique(goodx); 
  if length(goodx) > 0
    fprintf(1,'scattered interp T %3i of %3i has %6i finite values of %6i \n',iix,55,length(goodx),72*45);
    x = mls_LAT(goodx,goody); x = double(x(:)); 
    y = mls_LON(goodx,goody); y = double(y(:)); 
    z = junk(goodx,goody);    z = double(z(:)); 
    FT = scatteredInterpolant(x,y,z,'linear');
    aTnew(iix,:)  = single(FT(double(pL3.rlat),double(pL3.rlon)));
  end

  junk = squeeze(aW(:,:,iix)); [goodx,goody] = find(isfinite(junk)); 
  if length(goodx) > 0
    fprintf(1,'scattered interp WV %3i of %3i has %6i finite values of %6i \n',iix,55,length(goodx),72*45);
    x = mls_LAT(goodx,goody); x = double(x(:)); 
    y = mls_LON(goodx,goody); y = double(y(:)); 
    z = junk(goodx,goody);    z = double(z(:)); 
    FW = scatteredInterpolant(x,y,z,'linear');
    aWnew(iix,:)  = single(FW(double(pL3.rlat),double(pL3.rlon)));
  end

  junk = squeeze(aO3(:,:,iix)); [goodx,goody] = find(isfinite(junk)); 
  if length(goodx) > 0
    fprintf(1,'scattered interp O3 %3i of %3i has %6i finite values of %6i \n',iix,55,length(goodx),72*45);
    x = mls_LAT(goodx,goody); x = double(x(:)); 
    y = mls_LON(goodx,goody); y = double(y(:)); 
    z = junk(goodx,goody);    z = double(z(:)); 
    FO3 = scatteredInterpolant(x,y,z,'linear');
    aO3new(iix,:)  = single(FT(double(pL3.rlat),double(pL3.rlon)));
  end
end
