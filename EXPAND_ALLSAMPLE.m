% function EXPAND_ALLSAMPLE(dir,burnin)
function EXPAND_ALLSAMPLE(cha,burnin)
% burnin: enter for percent scale
% file=[dir,'/CHA_test.mat'];
% load(file)
% 
% load('./Result/Test_06/CHA_test.mat'); % test
% 
NIT=length(cha.McCOMPRESS);
NCH=size(cha.MpCOMPRESS(1).SMPMp,2);
NPOL=length(cha.MpCOMPRESS(1).NPOL);
NFLT=length(cha.McCOMPRESS(1).NFLT);
SUMPOL=zeros(NPOL,1);
SUMPOLPAIR=zeros(NPOL,NPOL);
SUMFLT=zeros(NFLT,1);
SUMFLTPAIR=zeros(NFLT,NFLT);
NDATAPOL=zeros(NPOL,1);
NDATAFLT=zeros(NFLT,1);
% 
Mpbin=[-1E-7:1E-10:1E-7];
Mcbin=[-1:0.001:1];
SMPINT=50; % sampling interval
MpHIST=zeros(NPOL,size(Mpbin,2)-1);
McHIST=zeros(NFLT,size(Mcbin,2)-1);
BURNIN=floor(burnin*NIT/100)+1;
SMPID=[1:SMPINT:NCH];
smppol=zeros(NPOL,NCH);
smpflt=zeros(NFLT,NCH);
SMPPOL=[];
SMPFLT=[];
% infid=zeros(NPOL,1);
% mpscale=zeros(NPOL,1);
% mpmin=zeros(NPOL,1);
for ii=BURNIN:NIT
  for jj=1:NPOL
    infid=cha.MpCOMPRESS(ii).NPOL(jj).Mpscale==Inf;
    if ~infid
      smppol(jj,:)=double((cha.MpCOMPRESS(ii).SMPMp+128))./(2.55.*cha.MpCOMPRESS(ii).NPOL(jj).Mpscale)+cha.MpCOMPRESS(ii).NPOL(jj).MpMIN;
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(smppol,Mpbin);
      NDATAPOL(jj)=NDATAPOL(jj)+NCH;
    else
      smppol=ones(1,NCH).*cha.MpCOMPRESS(ii).NPOL(jj).MpMAX;
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(smppol,Mpbin);
      NDATAPOL(jj)=NDATAPOL(jj)+NCH;
    end
  end
  SUMPOL=SUMPOL+NCH.*cha.MpCOMPRESS(ii).MEANMp;
  SUMPOLPAIR=SUMPOLPAIR+(NCH-1).*cha.MpCOMPRESS(ii).COVMp+NCH.*cha.MpCOMPRESS(ii).MEANMp*cha.MpCOMPRESS(ii).MEANMp';
  SMPPOL=[SMPPOL smppol(:,SMPID)];
  for kk=1:NFLT
    infid=cha.McCOMPRESS(ii).NFLT(kk).Mcscale==Inf;
    if ~infid
      smpflt(kk,:)=double((cha.McCOMPRESS(ii).SMPMc+128))./(2.55.*cha.McCOMPRESS(ii).NFLT(kk).Mcscale)+cha.McCOMPRESS(ii).NFLT(kk).McMIN;
      MpHIST(kk,:)=McHIST(kk,:)+histcounts(smpflt,Mcbin);
      NDATAFLT(kk)=NDATAFLT(kk)+NCH;
    else
      smpflt=ones(1,NCH).*cha.McCOMPRESS(ii).NFLT(kk).McMAX;
      McHIST(kk,:)=McHIST(kk,:)+histcounts(smpflt,Mcbin);
      NDATAFLT(kk)=NDATAFLT(kk)+NCH;
    end
  end
  SUMFLT=SUMFLT+NCH.*cha.McCOMPRESS(ii).MEANMc;
  SUMFLTPAIR=SUMFLTPAIR+(NCH-1).*cha.McCOMPRESS(ii).COVMc+NCH.*cha.McCOMPRESS(ii).MEANMc*cha.McCOMPRESS(ii).MEANMc';
  SMPFLT=[SMPFLT smpflt(:,SMPID)];
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
TCHA.SMPPOL=SMPPPOL;
TCHA.SMPFLT=SMPPFLT;
% 
outfile=[dir,'/TCHA.mat'];
save(outfile,'TCHA');
% save('./Result/Test_06/TCHA_test.mat','TCHA') % test

end
