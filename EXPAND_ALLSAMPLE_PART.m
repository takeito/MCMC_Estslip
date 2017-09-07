function EXPAND_ALLSAMPLE_PART(DIR,burnin)
% function EXPAND_ALLSAMPLE(cha,burnin)
% burnin: enter for percent scale
[INPUT]=READ_CHAFILE(DIR);
[TCHA]=cal_avestdbin(INPUT,burnin);
WRITE_TCHA(TCHA,DIR)
end
%% Read all cha mat file
function [INPUT]=READ_CHAFILE(DIR)
EXT='CHA_test*.mat';
file=dir([DIR,'/',EXT]);
[Nit,~]=size(file);
for ii=1:Nit
  INPUT(ii).fname=fullfile(file(ii).folder,file(ii).name);
end
end
%% Save routine.
function WRITE_TCHA(TCHA,DIR)
outfile=[DIR,'/TCHA.mat'];
save(outfile,'TCHA','-v7.3');
end
%% Calculate average, std, bin
function [TCHA]=cal_avestdbin(INPUT,burnin)
% load('./Result/Test_06/CHA_test.mat'); % test
NIT=size(INPUT,2);
Mpbin=[-1E-7:1E-9:1E-7];
Mcbin=[-1:0.01:1];
SMPINT=50; % sampling interval
BURNIN=floor(burnin*NIT/100)+1;
ACCTOTAL=0;
for ii=1:NIT
  load(INPUT(ii).fname);
  ACCFLAG=isfield(cha,'AJR');
  if ii==1
    NCH=size(cha.MpCOMPRESS.SMPMp,2);
    NPOL=length(cha.MpCOMPRESS.NPOL);
    NFLT=length(cha.McCOMPRESS.NFLT);
    SUMPOL=zeros(NPOL,1);
    SUMFLT=zeros(NFLT,1);    
    SUMPOLPAIR=zeros(NPOL,NPOL);
    SUMFLTPAIR=zeros(NFLT,NFLT);
    NDATAPOL=zeros(NPOL,1);
    NDATAFLT=zeros(NFLT,1);
    MpHIST=zeros(NPOL,size(Mpbin,2)-1);
    McHIST=zeros(NFLT,size(Mcbin,2)-1);
    SMPID=[1:SMPINT:NCH];
    smppol=zeros(NPOL,NCH);
    smpflt=zeros(NFLT,NCH);
    MEDPOL=zeros(NPOL,1);
    MEDFLT=zeros(NFLT,1);
    SMPPOL=[];
    SMPFLT=[];
  end
  if ii>BURNIN
    for jj=1:NPOL
      infid=cha.MpCOMPRESS.NPOL(jj).Mpscale==Inf;
      if ~infid
        smppol(jj,:)=(double(cha.MpCOMPRESS.SMPMp(jj,:))+128)./(2.55.*cha.MpCOMPRESS.NPOL(jj).Mpscale)+cha.MpCOMPRESS.NPOL(jj).MpMIN;
      else
        smppol(jj,:)=ones(1,NCH).*cha.MpCOMPRESS.NPOL(jj).MpMAX;
      end
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(smppol(jj,:),Mpbin);
      NDATAPOL(jj)=NDATAPOL(jj)+NCH;
    end
    for kk=1:NFLT
      infid=cha.McCOMPRESS.NFLT(kk).Mcscale==Inf;
      if ~infid
        smpflt(kk,:)=(double(cha.McCOMPRESS.SMPMc(kk,:))+128)./(2.55.*cha.McCOMPRESS.NFLT(kk).Mcscale)+cha.McCOMPRESS.NFLT(kk).McMIN;
      else
        smpflt(kk,:)=ones(1,NCH).*cha.McCOMPRESS.NFLT(kk).McMAX;
      end
      McHIST(kk,:)=McHIST(kk,:)+histcounts(smpflt(kk,:),Mcbin);
      NDATAFLT(kk)=NDATAFLT(kk)+NCH;
    end
    SUMPOL=SUMPOL+NCH.*cha.MpCOMPRESS.MEANMp;
    SUMFLT=SUMFLT+NCH.*cha.McCOMPRESS.MEANMc;
    SUMPOLPAIR=SUMPOLPAIR+(NCH-1).*cha.MpCOMPRESS.COVMp+NCH.*cha.MpCOMPRESS.MEANMp*cha.MpCOMPRESS.MEANMp';
    SUMFLTPAIR=SUMFLTPAIR+(NCH-1).*cha.McCOMPRESS.COVMc+NCH.*cha.McCOMPRESS.MEANMc*cha.McCOMPRESS.MEANMc';    
  else
    for jj=1:NPOL
      infid=cha.MpCOMPRESS.NPOL(jj).Mpscale==Inf;
      if ~infid
        smppol(jj,:)=(double(cha.MpCOMPRESS.SMPMp(jj,:))+128)./(2.55.*cha.MpCOMPRESS.NPOL(jj).Mpscale)+cha.MpCOMPRESS.NPOL(jj).MpMIN;
      else
        smppol(jj,:)=ones(1,NCH).*cha.MpCOMPRESS.NPOL(jj).MpMAX;
      end
    end
    for kk=1:NFLT
      infid=cha.McCOMPRESS.NFLT(kk).Mcscale==Inf;
      if ~infid
        smpflt(kk,:)=(double(cha.McCOMPRESS.SMPMc(kk,:))+128)./(2.55.*cha.McCOMPRESS.NFLT(kk).Mcscale)+cha.McCOMPRESS.NFLT(kk).McMIN;
      else
        smpflt(kk,:)=ones(1,NCH).*cha.McCOMPRESS.NFLT(kk).McMAX;
      end
    end
  end
  if ACCFLAG
    ACCTOTAL=ACCTOTAL+cha.AJR;
  end
  SMPPOL=[SMPPOL smppol(:,SMPID)];
  SMPFLT=[SMPFLT smpflt(:,SMPID)];
  clear cha
  fprintf('Now finised at %i/%i\n',ii,NIT)
end
AVEPOL=SUMPOL./NDATAPOL;
AVEFLT=SUMFLT./NDATAFLT;
[~,MEDPOLID]=max(MpHIST,[],2);
[~,MEDFLTID]=max(McHIST,[],2);
for NP=1:NPOL
  MEDPOL(NP)=0.5*(Mpbin(MEDPOLID(NP))+Mpbin(MEDPOLID(NP)+1));
end
for NF=1:NFLT
  MEDFLT(NF)=0.5*(Mcbin(MEDFLTID(NF))+Mcbin(MEDFLTID(NF)+1));
end
COVPOL=SUMPOLPAIR./NDATAPOL-(SUMPOL./NDATAPOL)*(SUMPOL./NDATAPOL)';
COVFLT=SUMFLTPAIR./NDATAFLT-(SUMFLT./NDATAFLT)*(SUMFLT./NDATAFLT)';
STDPOL=diag(COVPOL);
STDFLT=diag(COVFLT);
CORPOL=COVPOL./(sqrt(STDPOL)*sqrt(STDPOL'));
CORFLT=COVFLT./(sqrt(STDFLT)*sqrt(STDFLT'));
% Output
if ACCFLAG
  TCHA.ACCTOTAL=ACCTOTAL;
end
TCHA.Burnin=burnin;
TCHA.Smpint=SMPINT;
TCHA.Mpbin=Mpbin;
TCHA.Mcbin=Mcbin;
TCHA.AVEPOL=single(AVEPOL);
TCHA.AVEFLT=single(AVEFLT);
TCHA.MEDPOL=single(MEDPOL);
TCHA.MEDFLT=single(MEDFLT);
TCHA.STDPOL=single(STDPOL);
TCHA.STDFLT=single(STDFLT);
TCHA.COVPOL=single(COVPOL);
TCHA.COVFLT=single(COVFLT);
TCHA.CORPOL=single(CORPOL);
TCHA.CORFLT=single(CORFLT);
TCHA.HISTPOL=single(MpHIST);
TCHA.HISTFLT=single(McHIST);
TCHA.NDATPOL=single(NDATAPOL);
TCHA.NDATFLT=single(NDATAFLT);
TCHA.SMPPOL=single(SMPPOL);
TCHA.SMPFLT=single(SMPFLT);
end
