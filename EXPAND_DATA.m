function EXPAND_DATA(dir)
% 
file=[dir,'/CHA_test.mat'];
load(file)
% 
% load('./Result/Test_06/CHA_test.mat'); % test
% 
NIT=length(CHA.McCOMPRESS);
NPOL=length(CHA.MpCOMPRESS(1).NPOL);
NFLT=length(CHA.McCOMPRESS(1).NFLT);
sumpol=zeros(NPOL,1);
SUMPOL=zeros(NPOL,1);
SUMPOL2=zeros(NPOL,1);
SUMPOLPAIR=zeros(NPOL,NPOL);
sumflt=zeros(NFLT,1);
SUMFLT=zeros(NFLT,1);
SUMFLT2=zeros(NFLT,1);
SUMFLTPAIR=zeros(NFLT,NFLT);
NDATAPOL=zeros(NPOL,1);
NDATAFLT=zeros(NFLT,1);
% 
binnum=[0.5:1:127];
Mpbin=[-10^9:10^7:10^9];
Mcbin=[-1:0.02:1];
MpHIST=zeros(NPOL,size(Mpbin,2));
McHIST=zeros(NPOL,size(Mcbin,2));

for ii=1:NIT
  for jj=1:NPOL
    infid=CHA.MpCOMPRESS(ii).NPOL(jj).Mpscale==Inf;
    estpol=[];
    if ~infid
      cbin=binnum./CHA.MpCOMPRESS(ii).NPOL(jj).Mpscale+CHA.MpCOMPRESS(ii).NPOL(jj).MpMIN;
      for ll=1:length(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST)
        estpol=[estpol ones(1,CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST(ll)).*cbin(ll)];
      end
      MpHIST(jj)=MpHIST(jj)+histcounts(estpol,Mpbin);
      ncha=sum(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
    else
      ncha=sum(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      estpol=ones(1,ncha).*CHA.MpCOMPRESS(ii).POL(jj).MpMAX;
      MpHIST(jj)=MpHIST(jj)+histcounts(estpol,Mpbin);
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
    end
  end
  SUMPOLPAIR=SUMPOLPAIR+(ncha-1).*CHA.MpCOMPRESS(ii).COVMp+ncha.*CHA.MpCOMPRESS(ii).MEANMp*CHA.MpCOMPRESS(ii).MEANMp';
  for kk=1:NFLT
    infid=CHA.McCOMPRESS(ii).NFLT(kk).Mcscale==Inf;
    if ~infid
      cbin=binnum./CHA.McCOMPRESS(ii).NFLT(kk).Mcscale+CHA.McCOMPRESS(ii).NFLT(kk).McMIN;
      datasum=sum(cbin.*CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      datasum2=sum(cbin.^2.*CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      ncha=sum(CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      sumflt(kk)=datasum;
      SUMFLT(kk)=SUMFLT(kk)+datasum;
      SUMFLT2(kk)=SUMFLT2(kk)+datasum2;
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    else
      ncha=sum(CHA.McCOMPRESS(ii).NFLT(kk).MpHIST);
      datasum=ncha*CHA.McCOMPRESS(ii).NFLT(kk).McMAX;
      datasum2=ncha*CHA.McCOMPRESS(ii).NFLT(kk).McMAX^2;
      sumflt(kk)=datasum;
      SUMFLT(kk)=SUMFLT(kk)+datasum;
      SUMFLT2(kk)=SUMFLT2(kk)+datasum2;
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    end
  end
  SUMFLTPAIR=SUMFLTPAIR+(ncha-1).*CHA.McCOMPRESS(ii).COVMc+ncha.*(sumflt./ncha)*(sumflt./ncha)';
end
% 
AVEPOL=SUMPOL./NDATAPOL;
AVEFLT=SUMFLT./NDATAFLT;
STDPOL=SUMPOL2./NDATAPOL-AVEPOL.^2;
STDFLT=SUMFLT2./NDATAFLT-AVEFLT.^2;
COVPOL=SUMPOLPAIR./NDATAPOL-AVEPOL*AVEPOL';
COVFLT=SUMFLTPAIR./NDATAFLT-AVEFLT*AVEFLT';
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
% 
outfile=[dir,'/TCHA.mat'];
save(outfile,'TCHA');
% save('./Result/Test_06/TCHA_test.mat','TCHA') % test

end