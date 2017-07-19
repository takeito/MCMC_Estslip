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
SUMPOL=zeros(NPOL,1);
SUMPOLPAIR=zeros(NPOL,NPOL);
SUMFLT=zeros(NFLT,1);
SUMFLTPAIR=zeros(NFLT,NFLT);
NDATAPOL=zeros(NPOL,1);
NDATAFLT=zeros(NFLT,1);
% 
binnum=[-127.5:1:127];
Mpbin=[-10^10:10^8:10^10];
Mcbin=[-1:0.01:1];
MpHIST=zeros(NPOL,size(Mpbin,2)-1);
McHIST=zeros(NFLT,size(Mcbin,2)-1);

for ii=1:NIT
  for jj=1:NPOL
    infid=CHA.MpCOMPRESS(ii).NPOL(jj).Mpscale==Inf;
    estpol=[];
    if ~infid
      cbin=(binnum+128)./(2.55.*CHA.MpCOMPRESS(ii).NPOL(jj).Mpscale)+CHA.MpCOMPRESS(ii).NPOL(jj).MpMIN;
      for ll=1:length(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST)
        estpol=[estpol ones(1,CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST(ll)).*cbin(ll)];
      end
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(estpol,Mpbin);
      ncha=sum(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
    else
      ncha=sum(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      estpol=ones(1,ncha).*CHA.MpCOMPRESS(ii).NPOL(jj).MpMAX;
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(estpol,Mpbin);
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
    end
  end
  SUMPOL=SUMPOL+ncha.*CHA.MpCOMPRESS(ii).MEANMp;
  SUMPOLPAIR=SUMPOLPAIR+(ncha-1).*CHA.MpCOMPRESS(ii).COVMp+ncha.*CHA.MpCOMPRESS(ii).MEANMp*CHA.MpCOMPRESS(ii).MEANMp';
  for kk=1:NFLT
    infid=CHA.McCOMPRESS(ii).NFLT(kk).Mcscale==Inf;
    estflt=[];
    if ~infid
      cbin=(binnum+128)./(2.55.*CHA.McCOMPRESS(ii).NFLT(kk).Mcscale)+CHA.McCOMPRESS(ii).NFLT(kk).McMIN;
      for mm=1:length(CHA.McCOMPRESS(ii).NFLT(kk).McHIST)
        estflt=[estflt ones(1,CHA.McCOMPRESS(ii).NFLT(kk).McHIST(mm)).*cbin(mm)];
      end
      McHIST(kk,:)=McHIST(kk,:)+histcounts(estflt,Mcbin);
      ncha=sum(CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    else
      ncha=sum(CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      estflt=ones(1,ncha).*CHA.McCOMPRESS(ii).NFLT(kk).McMAX;
      McHIST(kk,:)=McHIST(kk,:)+histcounts(estflt,Mcbin);
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    end
  end
  SUMFLT=SUMFLT+ncha.*CHA.McCOMPRESS(ii).MEANMc;
  SUMFLTPAIR=SUMFLTPAIR+(ncha-1).*CHA.McCOMPRESS(ii).COVMc+ncha.*CHA.McCOMPRESS(ii).MEANMc*CHA.McCOMPRESS(ii).MEANMc';
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
% 
outfile=[dir,'/TCHA.mat'];
save(outfile,'TCHA');
% save('./Result/Test_06/TCHA_test.mat','TCHA') % test

end