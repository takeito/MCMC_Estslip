function EXPAND_DATA(dir)
% 
[PRM]=READ_PARAMETERS(INPUT);
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
TCHA.NDATPOL=NDATAPOL;
TCHA.NDATFLT=NDATAFLT;
% 
outfile=[dir,'/TCHA.mat'];
save(outfile,'TCHA');
% save('./Result/Test_06/TCHA_test.mat','TCHA') % test

end
%% READ PARAMETER FOR MCMC Inversion 
function [PRM]=READ_PARAMETERS(INPUT)
% MCMC Inversion for Geodetic 
% Coded    by Takeo Ito 2011/11/08 (ver 1.0)
% Modified by Takeo Ito 2012/10/26 (ver 1.1)
% Modified by Takeo Ito 2015/11/11 (ver 1.2)
% Modified by Takeo Ito 2016/07/06 (ver 1.3)
%
Fid=fopen(INPUT.Parfile,'r');
PRM.HOME_D=pwd;
FileOBS=fscanf(Fid,'%s \n',[1,1]);
PRM.FileOBS=fullfile(PRM.HOME_D,FileOBS);
[~]=fgetl(Fid);
DIRBlock=fscanf(Fid,'%s \n',[1,1]);
PRM.DIRBlock=fullfile(PRM.HOME_D,DIRBlock);
[~]=fgetl(Fid);
DIRBlock_Interface=fscanf(Fid,'%s \n',[1,1]);
PRM.DIRBlock_Interface=fullfile(PRM.HOME_D,DIRBlock_Interface);
[~]=fgetl(Fid);
%
PRM.GPU=fscanf(Fid,'%d \n',[1,1]);
[~]=fgetl(Fid);
PRM.ITR=fscanf(Fid,'%d \n',[1,1]);
[~]=fgetl(Fid);
PRM.CHA=fscanf(Fid,'%d \n',[1,1]);
[~]=fgetl(Fid);
PRM.KEP=fscanf(Fid,'%d \n',[1,1]);
[~]=fgetl(Fid);
PRM.RWD=fscanf(Fid,'%f \n',[1,1]);
fclose(Fid);
%====================================================
tmp=load(INPUT.Optfile);
PRM.num=size(tmp,1);
PRM.OptB1=tmp(:,1);
PRM.OptB2=tmp(:,2);
PRM.OptINT=tmp(:,3);
%====================================================
fprintf('==================\nINPUT PARAMETERS\n==================\n') 
fprintf('HOME_D             : %s \n',PRM.HOME_D) 
fprintf('FileOBS            : %s \n',PRM.FileOBS) 
fprintf('DIRBlock           : %s \n',PRM.DIRBlock)
fprintf('DIRBlock_Interface : %s \n',PRM.DIRBlock_Interface) 
fprintf('GPUdev (CPU:99)    : %i \n',PRM.GPU) 
fprintf('ITR(Max_Nitr)      : %i \n',PRM.ITR) 
fprintf('CHA(Chain)         : %i \n',PRM.CHA) 
fprintf('KEP(KEEP)          : %i \n',PRM.KEP) 
fprintf('RWD(Walk_dis)      : %4.2f \n',PRM.RWD) 
fprintf('==================\n') 
%====================================================
disp('PASS READ_PARAMETERS')
end
