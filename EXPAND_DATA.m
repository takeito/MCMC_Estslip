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
SUMPOL2=zeros(NPOL,1);
SUMFLT=zeros(NFLT,1);
SUMFLT2=zeros(NFLT,1);
NDATAPOL=zeros(NPOL,1);
NDATAFLT=zeros(NFLT,1);
% 
binnum=[0.5:1:127];

for ii=1:NIT
  for jj=1:NPOL
    infid=CHA.MpCOMPRESS(ii).NPOL(jj).Mpscale==Inf;
    if ~infid
      cbin=binnum./CHA.MpCOMPRESS(ii).NPOL(jj).Mpscale+CHA.MpCOMPRESS(ii).NPOL(jj).MpMIN;
      datasum=sum(cbin.*CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      datasum2=sum(cbin.^2.*CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      ncha=sum(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      SUMPOL(jj)=SUMPOL(jj)+datasum;
      SUMPOL2(jj)=SUMPOL2(jj)+datasum2;
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
    else
      ncha=sum(CHA.MpCOMPRESS(ii).NPOL(jj).MpHIST);
      datasum=ncha*CHA.MpCOMPRESS(ii).NPOL(jj).MpMAX;
      datasum2=ncha*CHA.MpCOMPRESS(ii).NPOL(jj).MpMAX^2;
      SUMPOL(jj)=SUMPOL(jj)+datasum;
      SUMPOL2(jj)=SUMPOL2(jj)+datasum2;
      NDATAPOL(jj)=NDATAPOL(jj)+ncha;
    end
  end
  for kk=1:NFLT
    infid=CHA.McCOMPRESS(ii).NFLT(kk).Mcscale==Inf;
    if ~infid
      cbin=binnum./CHA.McCOMPRESS(ii).NFLT(kk).Mcscale+CHA.McCOMPRESS(ii).NFLT(kk).McMIN;
      datasum=sum(cbin.*CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      datasum2=sum(cbin.^2.*CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      ncha=sum(CHA.McCOMPRESS(ii).NFLT(kk).McHIST);
      SUMFLT(kk)=SUMFLT(kk)+datasum;
      SUMFLT2(kk)=SUMFLT2(kk)+datasum2;
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    else
      ncha=sum(CHA.McCOMPRESS(ii).NFLT(kk).MpHIST);
      datasum=ncha*CHA.McCOMPRESS(ii).NFLT(kk).McMAX;
      datasum2=ncha*CHA.McCOMPRESS(ii).NFLT(kk).McMAX^2;
      SUMFLT(kk)=SUMFLT(kk)+datasum;
      SUMFLT2(kk)=SUMFLT2(kk)+datasum2;
      NDATAFLT(kk)=NDATAFLT(kk)+ncha;
    end
  end
end
% 
AVEPOL=SUMPOL./NDATAPOL;
AVEFLT=SUMFLT./NDATAFLT;
STDPOL=SUMPOL2./NDATAPOL-AVEPOL.^2;
STDFLT=SUMFLT2./NDATAFLT-AVEFLT.^2;
% 
TCHA.AVEPOL=AVEPOL;
TCHA.AVEFLT=AVEFLT;
TCHA.STDPOL=STDPOL;
TCHA.STDFLT=STDFLT;
% 
outfile=[dir,'/TCHA.mat'];
save(outfile,'TCHA');
% save('./Result/Test_06/TCHA_test.mat','TCHA') % test

end