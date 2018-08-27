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
% sfactor=2^8;
sfactor=2^16;
Mpbin=[-1E-7:1E-10:1E-7];
Mcbin=[-1:1E-3:1];
Mibin=[-1E-7:1E-10:1E-7];
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
    NINE=length(cha.MiCOMPRESS.NINE);
    SUMPOL=zeros(NPOL,1);
    SUMFLT=zeros(NFLT,1);    
    SUMINE=zeros(NINE,1);    
    SUMPOLPAIR=zeros(NPOL,NPOL);
    SUMFLTPAIR=zeros(NFLT,NFLT);
    SUMINEPAIR=zeros(NINE,NINE);
    NDATAPOL=0;
    NDATAFLT=0;
    NDATAINE=0;
    MpHIST=zeros(NPOL,size(Mpbin,2)-1);
    McHIST=zeros(NFLT,size(Mcbin,2)-1);
    MiHIST=zeros(NINE,size(Mibin,2)-1);
    SMPID=[1:SMPINT:NCH];
    smppol=zeros(NPOL,NCH);
    smpflt=zeros(NFLT,NCH);
    smpine=zeros(NINE,NCH);
    MEDPOL=zeros(NPOL,1);
    MEDFLT=zeros(NFLT,1);
    MEDINE=zeros(NINE,1);
    SMPPOL=[];
    SMPFLT=[];
    SMPINE=[];
  end
  if ii>BURNIN
    for jj=1:NPOL
      infid=cha.MpCOMPRESS.NPOL(jj).Mpscale==Inf;
      if ~infid
        smppol(jj,:)=(double(cha.MpCOMPRESS.SMPMp(jj,:))+(sfactor/2))./((sfactor-1).*cha.MpCOMPRESS.NPOL(jj).Mpscale)+cha.MpCOMPRESS.NPOL(jj).MpMIN;
      else
        smppol(jj,:)=ones(1,NCH).*cha.MpCOMPRESS.NPOL(jj).MpMAX;
      end
      MpHIST(jj,:)=MpHIST(jj,:)+histcounts(smppol(jj,:),Mpbin);
    end
    for kk=1:NFLT
      infid=cha.McCOMPRESS.NFLT(kk).Mcscale==Inf;
      if ~infid
        smpflt(kk,:)=(double(cha.McCOMPRESS.SMPMc(kk,:))+(sfactor/2))./((sfactor-1).*cha.McCOMPRESS.NFLT(kk).Mcscale)+cha.McCOMPRESS.NFLT(kk).McMIN;
      else
        smpflt(kk,:)=ones(1,NCH).*cha.McCOMPRESS.NFLT(kk).McMAX;
      end
      McHIST(kk,:)=McHIST(kk,:)+histcounts(smpflt(kk,:),Mcbin);
    end
    for ll=1:NINE
      infid=cha.MiCOMPRESS.NINE(ll).Miscale==Inf;
      if ~infid
        smpine(ll,:)=(double(cha.MiCOMPRESS.SMPMi(ll,:))+(sfactor/2))./((sfactor-1).*cha.MiCOMPRESS.NINE(ll).Miscale)+cha.MiCOMPRESS.NINE(ll).MiMIN;
      else
        smpine(ll,:)=ones(1,NCH).*cha.MiCOMPRESS.NINE(ll).MiMAX;
      end
      MiHIST(ll,:)=MiHIST(ll,:)+histcounts(smpine(ll,:),Mibin);
    end
    NDATAPOL=NDATAPOL+NCH;
    NDATAFLT=NDATAFLT+NCH;
    NDATAINE=NDATAINE+NCH;
    SUMPOL=SUMPOL+sum(smppol,2);
    SUMFLT=SUMFLT+sum(smpflt,2);
    SUMINE=SUMINE+sum(smpine,2);
    SUMPOLPAIR=SUMPOLPAIR+smppol*smppol';
    SUMFLTPAIR=SUMFLTPAIR+smpflt*smpflt';
    SUMINEPAIR=SUMINEPAIR+smpine*smpine';
  else
    for jj=1:NPOL
      infid=cha.MpCOMPRESS.NPOL(jj).Mpscale==Inf;
      if ~infid
        smppol(jj,:)=(double(cha.MpCOMPRESS.SMPMp(jj,:))+(sfactor/2))./((sfactor-1).*cha.MpCOMPRESS.NPOL(jj).Mpscale)+cha.MpCOMPRESS.NPOL(jj).MpMIN;
      else
        smppol(jj,:)=ones(1,NCH).*cha.MpCOMPRESS.NPOL(jj).MpMAX;
      end
    end
    for kk=1:NFLT
      infid=cha.McCOMPRESS.NFLT(kk).Mcscale==Inf;
      if ~infid
        smpflt(kk,:)=(double(cha.McCOMPRESS.SMPMc(kk,:))+(sfactor/2))./((sfactor-1).*cha.McCOMPRESS.NFLT(kk).Mcscale)+cha.McCOMPRESS.NFLT(kk).McMIN;
      else
        smpflt(kk,:)=ones(1,NCH).*cha.McCOMPRESS.NFLT(kk).McMAX;
      end
    end
    for ll=1:NINE
      infid=cha.MiCOMPRESS.NINE(ll).Miscale==Inf;
      if ~infid
        smpine(ll,:)=(double(cha.MiCOMPRESS.SMPMi(ll,:))+(sfactor/2))./((sfactor-1).*cha.MiCOMPRESS.NINE(ll).Miscale)+cha.MiCOMPRESS.NINE(ll).MiMIN;
      else
        smpine(ll,:)=ones(1,NCH).*cha.MiCOMPRESS.NINE(ll).MiMAX;
      end
    end
  end
  if ACCFLAG
    ACCTOTAL=ACCTOTAL+cha.AJR;
  end
  SMPPOL=[SMPPOL smppol(:,SMPID)];
  SMPFLT=[SMPFLT smpflt(:,SMPID)];
  SMPINE=[SMPINE smpine(:,SMPID)];
  clear cha
  fprintf('Now finised at %i/%i\n',ii,NIT)
end
AVEPOL=SUMPOL./NDATAPOL;
AVEFLT=SUMFLT./NDATAFLT;
AVEINE=SUMINE./NDATAINE;
[~,MEDPOLID]=max(MpHIST,[],2);
[~,MEDFLTID]=max(McHIST,[],2);
[~,MEDINEID]=max(MiHIST,[],2);
for NP=1:NPOL
  MEDPOL(NP)=0.5*(Mpbin(MEDPOLID(NP))+Mpbin(MEDPOLID(NP)+1));
end
for NF=1:NFLT
  MEDFLT(NF)=0.5*(Mcbin(MEDFLTID(NF))+Mcbin(MEDFLTID(NF)+1));
end
for NI=1:NINE
  MEDINE(NI)=0.5*(Mibin(MEDINEID(NI))+Mibin(MEDINEID(NI)+1));
end
COVPOL=SUMPOLPAIR./NDATAPOL-AVEPOL*AVEPOL';
COVFLT=SUMFLTPAIR./NDATAFLT-AVEFLT*AVEFLT';
COVINE=SUMINEPAIR./NDATAINE-AVEINE*AVEINE';
STDPOL=sqrt(diag(COVPOL));
STDFLT=sqrt(diag(COVFLT));
STDINE=sqrt(diag(COVINE));
CORPOL=COVPOL./(STDPOL*STDPOL');
CORFLT=COVFLT./(STDFLT*STDFLT');
CORINE=COVINE./(STDINE*STDINE');
% Output
if ACCFLAG
  TCHA.ACCTOTAL=ACCTOTAL;
end
TCHA.Burnin=burnin;
TCHA.Smpint=SMPINT;
TCHA.Mpbin=Mpbin;
TCHA.Mcbin=Mcbin;
TCHA.Mibin=Mibin;
TCHA.AVEPOL=single(AVEPOL);
TCHA.AVEFLT=single(AVEFLT);
TCHA.AVEINE=single(AVEINE);
TCHA.MEDPOL=single(MEDPOL);
TCHA.MEDFLT=single(MEDFLT);
TCHA.MEDINE=single(MEDINE);
TCHA.STDPOL=single(STDPOL);
TCHA.STDFLT=single(STDFLT);
TCHA.STDINE=single(STDINE);
TCHA.COVPOL=single(COVPOL);
TCHA.COVFLT=single(COVFLT);
TCHA.COVINE=single(COVINE);
TCHA.CORPOL=single(CORPOL);
TCHA.CORFLT=single(CORFLT);
TCHA.CORINE=single(CORINE);
TCHA.HISTPOL=single(MpHIST);
TCHA.HISTFLT=single(McHIST);
TCHA.HISTINE=single(MiHIST);
TCHA.NDATPOL=single(NDATAPOL);
TCHA.NDATFLT=single(NDATAFLT);
TCHA.NDATINE=single(NDATAINE);
TCHA.SMPPOL=single(SMPPOL);
TCHA.SMPFLT=single(SMPFLT);
TCHA.SMPINE=single(SMPINE);
end
