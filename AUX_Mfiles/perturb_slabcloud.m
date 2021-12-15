function pnew = perturb_slab_cloud(h,ha,p,pa,pertSTR,klayers,sarta);

fip = mktemp('fip.rtp');
fop = mktemp('fop.rtp');
frp = mktemp('frp.rtp');
ugh1 = mktemp('ugh1');
ugh2 = mktemp('ugh2');

pPERT = p;

%% pertSTR = PTFAS;    % P = phase (1=water,2=ice) , T=cprtop,  F=fraction,  A=amount, S=size

if pertSTR(1) == '1'
  %% doing water
  ctype1 = find(p.ctype == 101);
  ctype2 = find(p.ctype2 == 101);
  if pertSTR(2) == '1'
    disp('water cprtop pert')
    dp = pPERT.cprbot(ctype1)-pPERT.cprtop(ctype1);
      pPERT.cprtop(ctype1)  = pPERT.cprtop(ctype1)  * 0.9;   %% move things UP by 10%
      pPERT.cprbot(ctype1)  = pPERT.cprtop(ctype1) + dp;
    dp = pPERT.cprbot2(ctype2)-pPERT.cprtop2(ctype2);
      pPERT.cprtop2(ctype2) = pPERT.cprtop2(ctype2) * 0.9;   %% move things UP by 10%
      pPERT.cprbot2(ctype2) = pPERT.cprtop2(ctype2) + dp;
  elseif pertSTR(2) == '2'
    disp('water cprtop pert by 20%')
    dp = pPERT.cprbot(ctype1)-pPERT.cprtop(ctype1);
      pPERT.cprtop(ctype1)  = pPERT.cprtop(ctype1)  * 0.8;   %% move things UP by 20%
      pPERT.cprbot(ctype1)  = pPERT.cprtop(ctype1) + dp;
    dp = pPERT.cprbot2(ctype2)-pPERT.cprtop2(ctype2);
      pPERT.cprtop2(ctype2) = pPERT.cprtop2(ctype2) * 0.8;   %% move things UP by 20%
      pPERT.cprbot2(ctype2) = pPERT.cprtop2(ctype2) + dp;
  elseif pertSTR(3) == '1'
    disp('water cfrac pert')
    pPERT.cprtop(ctype1)  = max(1,pPERT.cprtop(ctype1)  * 1.1);   %% increase frac by 10%
    pPERT.cprtop2(ctype2) = max(1,pPERT.cprtop2(ctype2) * 1.1);   %% increase frac by 10%
  elseif pertSTR(4) == '1'
    disp('water cngwat pert')
    pPERT.cngwat(ctype1)  = pPERT.cngwat(ctype1)  * 1.1;   %% increase amt by 10%
    pPERT.cngwat2(ctype2) = pPERT.cngwat2(ctype2) * 1.1;   %% increase amt by 10%
  elseif pertSTR(5) == '1'
    disp('water cpsize pert')
    pPERT.cpsize(ctype1)  = pPERT.cpsize(ctype1)  * 1.1;   %% increase size by 10%
    pPERT.cpsize2(ctype2) = pPERT.cpsize2(ctype2) * 1.1;   %% increase size  by 10%
  end
elseif pertSTR(1) == '2'
  %% doing ice
  ctype1 = find(p.ctype == 201);
  ctype2 = find(p.ctype2 == 201);
  if pertSTR(2) == '1'
    disp('ice cprtop pert')
    dp = pPERT.cprbot(ctype1)-pPERT.cprtop(ctype1);
      pPERT.cprtop(ctype1)  = pPERT.cprtop(ctype1)  * 0.9;   %% move things UP by 10%
      pPERT.cprbot(ctype1)  = pPERT.cprtop(ctype1) + dp;
    dp = pPERT.cprbot2(ctype2)-pPERT.cprtop2(ctype2);
      pPERT.cprtop2(ctype2) = pPERT.cprtop2(ctype2) * 0.9;   %% move things UP by 10%
      pPERT.cprbot2(ctype2) = pPERT.cprtop2(ctype2) + dp;
  elseif pertSTR(2) == '2'
    disp('ice cprtop pert by 20%')
    dp = pPERT.cprbot(ctype1)-pPERT.cprtop(ctype1);
      pPERT.cprtop(ctype1)  = pPERT.cprtop(ctype1)  * 0.8;   %% move things UP by 20%
      pPERT.cprbot(ctype1)  = pPERT.cprtop(ctype1) + dp;
    dp = pPERT.cprbot2(ctype2)-pPERT.cprtop2(ctype2);
      pPERT.cprtop2(ctype2) = pPERT.cprtop2(ctype2) * 0.8;   %% move things UP by 20%
      pPERT.cprbot2(ctype2) = pPERT.cprtop2(ctype2) + dp;
  elseif pertSTR(3) == '1'
    disp('ice cfrac pert')
    pPERT.cprtop(ctype1)  = max(1,pPERT.cprtop(ctype1)  * 1.1);   %% increase frac by 10%
    pPERT.cprtop2(ctype2) = max(1,pPERT.cprtop2(ctype2) * 1.1);   %% increase frac by 10%
  elseif pertSTR(4) == '1'
    disp('ice cngwat pert')
    pPERT.cngwat(ctype1)  = pPERT.cngwat(ctype1)  * 1.1;   %% increase amt by 10%
    pPERT.cngwat2(ctype2) = pPERT.cngwat2(ctype2) * 1.1;   %% increase amt by 10%
  elseif pertSTR(5) == '1'
    disp('ice cpsize pert')
    pPERT.cpsize(ctype1)  = pPERT.cpsize(ctype1)  * 1.1;   %% increase size by 10%
    pPERT.cpsize2(ctype2) = pPERT.cpsize2(ctype2) * 1.1;   %% increase size  by 10%
  end
else
  error('hmm .. do not understand cloud type')
end

rtpwrite(fip,h,ha,pPERT,pa);

klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh1];
fprintf(1,'running klayers : %s \n',klayerser)
eval(klayerser)

sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh1];
fprintf(1,'running sarta : %s \n',sartaer)
eval(sartaer)

[hx,hax,px,pax] = rtpread(frp);
pnew = pPERT;
pnew.rcalc = px.rcalc;


rmer = ['!/bin/rm '  fip ' ' fop ' ' frp ' ' ugh1 ' ' ugh2];
eval(rmer);