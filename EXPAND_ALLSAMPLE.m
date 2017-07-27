function EXPAND_ALLSAMPLE(dir,burnin)
% burnin: enter for percent scale
file=[dir,'/CHA_test.mat'];
load(file)
% 
% load('./Result/Test_06/CHA_test.mat'); % test
% 
NIT=length(cha.McCOMPRESS);
NCH=size(cha.McCOMPRESS(1).SMPMp,2);
NPOL=length(cha.MpCOMPRESS(1).NPOL);
NFLT=length(cha.McCOMPRESS(1).NFLT);
SUMPOL=zeros(NPOL,1);
SUMPOLPAIR=zeros(NPOL,NPOL);
SUMFLT=zeros(NFLT,1);
SUMFLTPAIR=zeros(NFLT,NFLT);
NDATAPOL=zeros(NPOL,1);
NDATAFLT=zeros(NFLT,1);
% 
binnum=[-127.5:1:127];
Mpbin=[-1E-7:1E-10:1E-7];
Mcbin=[-1:0.001:1];
SMPINT=50; % sampling interval
MpHIST=zeros(NPOL,size(Mpbin,2)-1);
McHIST=zeros(NFLT,size(Mcbin,2)-1);
BURNIN=floor(burnin*NIT/100)+1;
SMPID=[1:SMPINT:NCH];
SMPPOL=[];
SMPFLT=[];
for ii=BURNIN:NIT
  for jj=1:NPOL
    infid=cha.MpCOMPRESS(ii).NPOL(jj).Mpscale==Inf;
    estpol=[];
    if ~infid
      cbin=(binnum+128)./(2.55.*cha.MpCOMPRESS(ii).NPOL(jj).Mpscale)+cha.MpCOMPRESS(ii).NPOL(jj).MpMIN;
      for ll=1:length(cha.MpCOMPRESS(ii).NPOL(jj).MpHIST)
        estpol=[estpol ones(1,cha.MpCOMPRESS(ii).NPOL(jj).MpHIST(ll)).*cbin(ll)];
      end
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(estpol,Mpbin);
      ncha=sum(cha.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
      smppol=cha.MpCOMPRESS
    else
      ncha=sum(cha.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      estpol=ones(1,ncha).*cha.MpCOMPRESS(ii).NPOL(jj).MpMAX;
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(estpol,Mpbin);
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
    end
  end
  SUMPOL=SUMPOL+ncha.*cha.MpCOMPRESS(ii).MEANMp;
  SUMPOLPAIR=SUMPOLPAIR+(ncha-1).*cha.MpCOMPRESS(ii).COVMp+ncha.*cha.MpCOMPRESS(ii).MEANMp*cha.MpCOMPRESS(ii).MEANMp';
  SMPPOL=[SMPPOL (cha.MpCOMPRESS(ii).SMPMp(:,SMPID)+128)./(2.55.*cha.MpCOMPRESS.Mpscale)];
  for kk=1:NFLT
    infid=cha.McCOMPRESS(ii).NFLT(kk).Mcscale==Inf;
    estflt=[];
    if ~infid
      cbin=(binnum+128)./(2.55.*cha.McCOMPRESS(ii).NFLT(kk).Mcscale)+cha.McCOMPRESS(ii).NFLT(kk).McMIN;
      for mm=1:length(cha.McCOMPRESS(ii).NFLT(kk).McHIST)
        estflt=[estflt ones(1,cha.McCOMPRESS(ii).NFLT(kk).McHIST(mm)).*cbin(mm)];
      end
      McHIST(kk,:)=McHIST(kk,:)+histcounts(estflt,Mcbin);
      ncha=sum(cha.McCOMPRESS(ii).NFLT(kk).McHIST);
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    else
      ncha=sum(cha.McCOMPRESS(ii).NFLT(kk).McHIST);
      estflt=ones(1,ncha).*cha.McCOMPRESS(ii).NFLT(kk).McMAX;
      McHIST(kk,:)=McHIST(kk,:)+histcounts(estflt,Mcbin);
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    end
  end
  SUMFLT=SUMFLT+ncha.*cha.McCOMPRESS(ii).MEANMc;
  SUMFLTPAIR=SUMFLTPAIR+(ncha-1).*cha.McCOMPRESS(ii).COVMc+ncha.*cha.McCOMPRESS(ii).MEANMc*cha.McCOMPRESS(ii).MEANMc';
  SMPFLT=[SMPFLT cha.McCOMPRESS(ii).Mcscale.*cha.McCOMPRESS(ii).SMPMc(:,SMPID)];
end
% 
AVEPOL=SUMPOL./NDATAPOL;
AVEFLT=SUMFLT./NDATAFLT;
COVPOL=SUMPOLPAIR./NDATAPOL-(SUMPOL./NDATAPOL)*(SUMPOL./NDATAPOL)';
COVFLT=SUMFLTPAIR./NDATAFLT-(SUMFLT./NDATAFLT)*(SUMFLT./NDATAFLT)';
STDPOL=diag(COVPOL);
STDFLT=diag(COVFLT);
CORPOL=COVPOL./(sqrt(STDPOL)*sqrt(STDPOL'));
CORFLT=COVFLT./(sqrt(STDFLT)*sqrt(STDFLT'));
% 
TCHA.AVEPOL=AVEPOL;
TCHA.AVEFLT=AVEFLT;
TCHA.STDPOL=STDPOL;
TCHA.STDFLT=STDFLT;
TCHA.COVPOL=COVPOL;
TCHA.COVFLT=COVFLT;
TCHA.CORPOL=CORPOL;
TCHA.CORFLT=CORFLT;
TCHA.HISTPOL=MpHIST;
TCHA.HISTFLT=McHIST;
TCHA.NDATPOL=NDATAPOL;
TCHA.NDATFLT=NDATAFLT;
% 
outfile=[dir,'/TCHA.mat'];
save(outfile,'TCHA');
% save('./Result/Test_06/TCHA_test.mat','TCHA') % test

end
