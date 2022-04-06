for iix = 1 : 55
  junk = squeeze(aT(:,:,iix));  [goodx,goody] = find(isfinite(junk)); mls_tlev(iix) = length(goodx);
  junk = squeeze(aW(:,:,iix));  [goodx,goody] = find(isfinite(junk)); mls_Wlev(iix) = length(goodx);
  junk = squeeze(aO3(:,:,iix)); [goodx,goody] = find(isfinite(junk)); mls_O3lev(iix) = length(goodx);
end
boo  = find(mls_tlev > 0      & mls_Wlev > 0      & mls_O3lev > 0);          %% shows lots of bad NaNs
boo2 = find(mls_tlev == 45*72 & mls_Wlev == 45*72 & mls_O3lev == 45*72);     %% shows lots of bad NaNs

for iix = 1 : 55    
  FT = griddedInterpolant(mls_LAT,mls_LON,squeeze(aT(:,:,iix)),'linear');
  FW = griddedInterpolant(mls_LAT,mls_LON,squeeze(aW(:,:,iix)),'linear');
  FO3 = griddedInterpolant(mls_LAT,mls_LON,squeeze(aO3(:,:,iix)),'linear');
  aTnew(iix,:)  = FT(pL3.rlat,pL3.rlon);
  aWnew(iix,:)  = FW(pL3.rlat,pL3.rlon);
  aO3new(iix,:) = FO3(pL3.rlat,pL3.rlon);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% so do this
aTx = aT;
aWx = aW;
aO3x = aO3;
for iix = boo(1):boo(end)
  junk = squeeze(aT(:,:,iix));  [goodx,goody] = find(isfinite(junk)); goodx = unique(goodx);
    if sum(setdiff(goodx',3:43)) ~= 0
      fprintf(1,'WARNING aT(:,:,%3i) expecting NaNs at 1,2 44,45 but that is not happening \n',iix);
    end
    junknew = junk; junknew(1,:) = junk(3,:);  junknew(2,:) = junk(3,:); junknew(44,:) = junk(42,:);  junknew(45,:) = junk(42,:); 
    aTx(:,:,iix) = junknew;
  junk = squeeze(aW(:,:,iix));  [goodx,goody] = find(isfinite(junk)); goodx = unique(goodx);
    if sum(setdiff(goodx',3:43)) ~= 0
      fprintf(1,'WARNING aW(:,:,%3i) expecting NaNs at 1,2 44,45 but that is not happening \n',iix);
    end
    junknew = junk; junknew(1,:) = junk(3,:);  junknew(2,:) = junk(3,:); junknew(44,:) = junk(42,:);  junknew(45,:) = junk(42,:); 
    aWx(:,:,iix) = junknew;
  junk = squeeze(aO3(:,:,iix));  [goodx,goody] = find(isfinite(junk)); goodx = unique(goodx);
    if sum(setdiff(goodx',3:43)) ~= 0
      fprintf(1,'WARNING aO3(:,:,%3i) expecting NaNs at 1,2 44,45 but that is not happening \n',iix);
    end
    junknew = junk; junknew(1,:) = junk(3,:);  junknew(2,:) = junk(3,:); junknew(44,:) = junk(42,:);  junknew(45,:) = junk(42,:); 
    aO3x(:,:,iix) = junknew;
end

for iix = 1 : 55
  junk = squeeze(aTx(:,:,iix));  [goodx,goody] = find(isfinite(junk)); mls_tlevx(iix) = length(goodx);
  junk = squeeze(aWx(:,:,iix));  [goodx,goody] = find(isfinite(junk)); mls_Wlevx(iix) = length(goodx);
  junk = squeeze(aO3x(:,:,iix)); [goodx,goody] = find(isfinite(junk)); mls_O3levx(iix) = length(goodx);
end
boox  = find(mls_tlevx > 0      & mls_Wlevx > 0     & mls_O3levx > 0);        %% shows lots of no NaNs
boo2x = find(mls_tlevx == 45*72 & mls_Wlevx == 45*72 & mls_O3levx == 45*72);  %% shows lots of no NaNs

for iix = 1 : 55    
  FT = griddedInterpolant(mls_LAT,mls_LON,squeeze(aTx(:,:,iix)),'linear');
  FW = griddedInterpolant(mls_LAT,mls_LON,squeeze(aWx(:,:,iix)),'linear');
  FO3 = griddedInterpolant(mls_LAT,mls_LON,squeeze(aO3x(:,:,iix)),'linear');
  aTnew(iix,:)  = FT(pL3.rlat,pL3.rlon);
  aWnew(iix,:)  = FW(pL3.rlat,pL3.rlon);
  aO3new(iix,:) = FO3(pL3.rlat,pL3.rlon);
end

