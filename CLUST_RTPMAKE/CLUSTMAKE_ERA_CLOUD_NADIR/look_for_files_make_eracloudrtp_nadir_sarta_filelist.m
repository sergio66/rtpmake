clear all

iRemove = -1;
iRemove = +1;

if iRemove > 0
  iYes = input('are you sure you want to remove zero size files (-1/+1) : ');
  if iYes < 0
    iRemove = -1
  end
end

figure(1); clf;
figure(2); clf;
figure(3); clf;
iCnt = 0;

yyS = 2002; yyE = 2002;
yyS = 2007; yyE = 2007;

iaY = input('Enter Start Stop years to check : ');
yyS = iaY(1);
yyE = iaY(2);

iTotal = 0;
iTotalExpect = 0;
for yy = yyS : yyE
  mmS = 1;
  mmE = 12;
  if mod(yy,4) == 0
    mos = [31 29 31 30 31 30 31 31 30 31 30 31];  %% leap year
  else
    mos = [31 28 31 30 31 30 31 31 30 31 30 31];  %% normal year
  end
  if yy == 2002
    mmS = 09;
  elseif yy == 2014
    mmE = 08;
  end
  for mm = mmS : mmE
    for dd = 1 : mos(mm)
      iTotalExpect = iTotalExpect + 24;
      iCnt = iCnt + 1;
      yyx(iCnt) = yy;
      mmx(iCnt) = mm;
      ddx(iCnt) = dd;
      ystr = num2str(yy,'%04d');
      mstr = num2str(mm,'%02d');
      dstr = num2str(dd,'%02d');
      slashstr = [ystr '/' mstr '/' dstr];
      dotstr   = [ystr '.' mstr '.' dstr];

      daysSOfar = dd + sum(mos(1:mm-1));
      dirAIRS = ['/asl/data/airs/AIRIBRAD/' num2str(yy,'%04d') '/' num2str(daysSOfar,'%03d') '/'];
      thedirAIRS = dir([dirAIRS '/AIRS.' dotstr '*.L1B.AIRS_Rad.v5.0.22*.hdf']);
      thedirAIRS = dir([dirAIRS '/AIRS.' dotstr '*.L1B.AIRS_Rad.v5.0.*.hdf']);      
      numAIRSfiles(iCnt) = length(thedirAIRS);

      dir0 = ['/asl/data/rtprod_airs/' slashstr '/'];
      filestr = [dir0 'sergio_nadir_cloudy_airs_l1b_era_sarta_baum_ice.' dotstr '.*.rtp'];
      filestr = [dir0 'rnd_nadir6track_cloudy_airs_l1b_era_sarta_baum_ice.' dotstr '.*.rtp'];
      thedir = dir(filestr);
      numfiles(iCnt) = length(thedir);
      sz = 0;
      zerosz = 0;
      for jj = 1 : length(thedir);
        sz  = sz + thedir(jj).bytes;
	if thedir(jj).bytes == 0
	  zerosz = zerosz + 1;
	  if iRemove > 0
       	    fname = [dir0 thedir(jj).name];
 	    rmer = ['!/bin/rm ' fname];
	    eval(rmer)
	  end
	end
      end
      iTotal = iTotal + length(thedir)-zerosz;
      avgsize(iCnt) = sz/(0.0000000000001+length(thedir)-zerosz);  %% dont count the size zero files!!!
      zerosize(iCnt) = zerosz;
    end %% dd loop
    
    figure(1); plot(1:iCnt,numfiles,'ko-',1:iCnt,zerosize,'r','linewidth',2); hold on
    figure(2); plot(1:iCnt,avgsize,'k', 'linewidth',2); hold on
    figure(3); plot(1:iCnt,numAIRSfiles,'k', 'linewidth',2); hold on    
    
    figure(1); line([iCnt iCnt],[0 24],'color','b');
    figure(2); line([iCnt iCnt],[0 10e7],'color','b');
    figure(3); line([iCnt iCnt],[0 240],'color','b');
    
    pause(0.1)
  end   %% mm loop
  figure(1); line([iCnt iCnt],[0 24],'color','r');
             text(iCnt-20,12,num2str(yy+1),'color','b','fontsize',10)  
  figure(2); line([iCnt iCnt],[0 10e7],'color','r');
             text(iCnt-20,5e7,num2str(yy+1),'color','b','fontsize',10)    
  figure(3); line([iCnt iCnt],[0 240],'color','r');
             text(iCnt-20,120,num2str(yy+1),'color','b','fontsize',10)    
  pause(0.1)
end     %% yy loop
fprintf(1,'made total of %5i files, expected %5i, AIRSlB expect %8.6f \n',iTotal,iTotalExpect,sum(numAIRSfiles)/10)
fprintf(1,'%5i zerosize files \n',sum(zerosize))

if (yyE == yyS)
  figure(1); text(1,12,num2str(yyS),'color','b','fontsize',10)
  figure(1); text(iCnt-20,12,num2str(yyS+1),'color','b','fontsize',10)

  figure(2); text(1,5e7,num2str(yyS),'color','b','fontsize',10)
  figure(2); text(iCnt-20,5e7,num2str(yyS+1),'color','b','fontsize',10)

  figure(3); text(1,120,num2str(yyS),'color','b','fontsize',10)
  figure(3); text(iCnt-20,120,num2str(yyS+1),'color','b','fontsize',10)
else
  figure(1); text(1,12,num2str(yyS),'color','b','fontsize',10)
  figure(1); text(iCnt-20,12,num2str(yyE+1),'color','b','fontsize',10)

  figure(2); text(1,5e7,num2str(yyS),'color','b','fontsize',10)
  figure(2); text(iCnt-20,5e7,num2str(yyE+1),'color','b','fontsize',10)

  figure(3); text(1,120,num2str(yyS),'color','b','fontsize',10)
  figure(3); text(iCnt-20,120,num2str(yyE+1),'color','b','fontsize',10)
end

figure(1); title('hourly files made'); hold off
figure(2); title('avg file size of each hour file'); hold off
figure(3); title('number of AIRS granule files found'); ax = axis; axis([ax(1) ax(2) 0 300]); hold off
