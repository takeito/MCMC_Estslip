function out_allanalysis_v3r(DIR)
% This is expander script for any TCHA.mat files which is computed with
% REMC method (using EST_BLOCK_slip_REMC)
Parfile='PARAMETER/parameter_export.txt';
[Par]=ReadPara(Parfile);
fprintf('Now loading %s ...',[DIR,'/PRM.mat'])
load([DIR,'/PRM.mat']);fprintf('load\n')
fprintf('Now loading %s ...',[DIR,'/OBS.mat'])
load([DIR,'/OBS.mat']);fprintf('load\n')
fprintf('Now loading %s ...',[DIR,'/BLK.mat'])
load([DIR,'/BLK.mat']);fprintf('load\n')
fprintf('Now loading %s ...',[DIR,'/TRI.mat'])
load([DIR,'/TRI.mat']);fprintf('load\n')
fprintf('Now loading %s ...',[DIR,'/GRN.mat'])
load([DIR,'/GRN.mat']);fprintf('load\n')
fprintf('Now loading %s ...',[DIR,'/TCHA.mat'])
load([DIR,'/TCHA.mat']);fprintf('load\n')
[G,DPT]=expand_parallelization(D,G);
% 
[SDR]=coupling2sdr(TCHA,D,DPT,G);
ExportCoupling(DIR,TCHA,BLK,SDR);
ExportInternalDeformation(DIR,TCHA,BLK);
for CP=1:size(Par.Coupling_Pair,2)
  ExportCouplingPair(DIR,BLK,TCHA,SDR,Par.Coupling_Pair(CP));
end
out_epole_allchain(DIR,TCHA,BLK,Par.BLKNAME);
[TRIg,~,GRD]=MAKE_PART_GREEN(BLK,Par.Grid_Setting);
out_vector_allchain_v2(DIR,BLK,TCHA,G,D,DPT,GRD,TRIg,OBS);
RelativeMotion_allchain(DIR,BLK,TCHA);
% % 
for EL=1:size(Par.Elastic_Pair,2)
  out_elastic_pair_allchain_v2(DIR,BLK,TCHA,G,D,OBS,GRD,TRI,TRIg,Par.Elastic_Pair(EL))  %OG(revised),OGnew
end
%
end
%% Expand and Parallelization 
function [G,DPT]=expand_parallelization(D,G)
G(1).TB=full(G(1).TB);
DPT.OBS=repmat(D(1).OBS,1,NReplica);
DPT.ERR=repmat(D(1).ERR,1,NReplica);
DPT.MID   = repmat(D(1).MID,1,NReplica);
DPT.CFINV = repmat(D(1).CFINV,1,NReplica);
end
%% Read export parameter file
function [PAR]=ReadPara(Parfile)
% Note:
% Prepare the export parameter file in the 'PARAMETER' folder as bellows,
%--example from here--
% # BLKNAME
% AM,PAC,OK,PHS,IMP
% # Coupling_Pair
% Inland
% 1 3
% 2 3
% 3 5
% 4 6
% # Elastic_Pair
% MTL
% 7 9
% 9 10
% PAC
% 11 12
% 10 11
% 9 11
% # GRID_SETTING
% 120 150 
% 20 50 
% 0.4
% --- END HERE ---
%--end of example--
PAR=[];
PDIR=pwd;
Parfile=fullfile(PDIR,Parfile);
Fid=fopen(Parfile,'r');
if Fid~=0
  PAR.BLKNAME=[];
  PAR.Coupling_Pair=[];
  PAR.Elastic_Pair=[];
  PAR.Grid_Setting=[];
  tline=fgetl(Fid);
  while 1
    switch tline
      case '# BLKNAME'
        while 1
          tline=fgetl(Fid);
          Tline=strtrim(strsplit(tline));
          if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
            PAR.BLKNAME=[PAR.BLKNAME, Tline];
          else
            break
          end
        end
      case '# Coupling_Pair'
        NCo=0;
        tline=fgetl(Fid);
        while 1
          Tline=strtrim(strsplit(tline));
          if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
            NCo=NCo+1;
%             tline=fscanf(Fid,'%s \n',[1,1]);
            PAR.Coupling_Pair(NCo).NAME=Tline;
            PAR.Coupling_Pair(NCo).pair=[];
            while 1
              tline=fgetl(Fid);
              Tline=strtrim(strsplit(tline));
              Tline=str2num(char(Tline));
              if ~isempty(Tline)
                PAR.Coupling_Pair(NCo).pair=[PAR.Coupling_Pair(NCo).pair; Tline'];
              else
                break
              end
            end
          else
            break
          end
        end
      case '# Elastic_Pair'
      NEl=0;
      tline=fgetl(Fid);
      while 1
        Tline=strtrim(strsplit(tline));
        if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
          NEl=NEl+1;
%           Tline=fscanf(Fid,'%s \n',[1,1]);
          PAR.Elastic_Pair(NEl).NAME=Tline;
          PAR.Elastic_Pair(NEl).pair=[];
          while 1
            tline=fgetl(Fid);
            Tline=strtrim(strsplit(tline));
            Tline=str2num(char(Tline));
            if ~isempty(Tline)
              PAR.Elastic_Pair(NEl).pair=[PAR.Elastic_Pair(NEl).pair; Tline'];
            else
              break
            end
          end
        else
          break
        end
      end
      case '# GRID_SETTING'
        while 1
          tline=fgetl(Fid);
          Tline=strtrim(strsplit(tline));
          if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
            Tline=str2num(char(Tline));
            if ~isempty(Tline)
                PAR.Grid_Setting=[PAR.Grid_Setting, Tline'];
            end
          else
            break
          end
        end
        otherwise
        tline=fgetl(Fid);
    end
    if strcmpi(tline,'--- END HERE ---'); break; end
  end
end

end
%% Export Euler pole for each block
function out_epole_allchain(DIR,TCHA,BLK,NAMEMAT)
% 
folder=[DIR,'/euler_pole'];
exid=exist(folder);
if exid~=7; mkdir(folder); end
NPOL=size(TCHA.AVEPOL,1)/TCHA.NReplica;
TCHA.AVEPOL=reshape(TCHA.AVEPOL,NPOL,TCHA.NReplica);
NR=0;
for REP=1:TCHA.NReplica
  subfolder=[folder,'/replica',num2str(REP)];
  exid=exist(subfolder);
  if exid~=7; mkdir(subfolder); end
  FID=fopen([subfolder,'/est_euler_pole.txt'],'w');
  fprintf(FID,'BLOCK_No. BLOCK_Name lat(deg) lon(deg) ang(deg/my) sigxx sigxy sigxz sigyy sigyz sigzz (1e-8 (rad/Myr)^2) \n');
  fprintf('BLOCK_No. BLOCK_Name lat(deg) lon(deg) ang(deg/my) sigxx sigxy sigxz sigyy sigyz sigzz (1e-8 (rad/Myr)^2) \n');
  if isempty(NAMEMAT)
    for kk=1:BLK(1).NBlock
      NAMEMAT{ii}=num2str(kk,'%02i');
    end
  end
  for BK=1:BLK(1).NBlock
      [latp,lonp,ang]=xyzp2lla(TCHA.AVEPOL(3.*BK-2,REP),TCHA.AVEPOL(3.*BK-1,REP),TCHA.AVEPOL(3.*BK,REP));
      [a,b,c,d,e,f]=out_cov(TCHA.COVPOL(NR+3.*BK-2:NR+3.*BK,NR+3.*BK-2:NR+3.*BK));
      fprintf('%2i %s %7.2f %8.2f %9.2e %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',...
          BK,NAMEMAT{BK},mean(latp),mean(lonp),mean(ang),a,b,c,d,e,f);
      fprintf(FID,'%2i %s %7.2f %8.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',...
          BK,NAMEMAT{BK},mean(latp),mean(lonp),mean(ang),a,b,c,d,e,f);
  end
  NR=NR+NPOL;
  fclose(FID);
end
% 
end
%% Export vectors of calculation, observation, residual vector
function out_vector_allchain_v2(DIR,BLK,TCHA,G,D,DPT,GRD,TRIg,OBS)
% 
calvec=calc_sampling_vector(OBS,BLK,TCHA,D,DPT,G);
resvec=calc_residual_vector(OBS,BLK,TCHA,calvec);
[grdvec,GRD]=calc_vector_atmesh(BLK,TCHA,D,G,GRD,TRIg);
% 
WRITE_VECTOR(TCHA,DIR,OBS,BLK,calvec,resvec,grdvec,GRD);
end
%% Export rigid and relative motion
function RelativeMotion_allchain(DIR,BLK,TCHA)
% 
Est_Motion_BLOCKS(DIR,TCHA,BLK)
% 
end
%% Export elastic vectors resulting from coupling at block boundaries
function out_elastic_pair_allchain_v2(DIR,BLK,TCHA,G,D,OBS,GRD,TRI,TRIg,Elastic_Pair)
NAME=Elastic_Pair.NAME{1};
PAIR=Elastic_Pair.pair;
% 
% calvec=calc_sampling_vector(OBS,BLK,TCHA,D,G);
[grdvec,GRD]=calc_vector_atmesh_pair(BLK,TCHA,D,G,GRD,TRIg,PAIR);
CALVECgrd.ela=[GRD(1).ALON;GRD(1).ALAT;grdvec.ELA(1:3:end)';grdvec.ELA(2:3:end)';grdvec.ELA(3:3:end)'];
[sitevec,OBS]=calc_vector_atmesh_pair(BLK,TCHA,D,DPT,G,OBS,TRI,PAIR);
CALVECsite.ela=[OBS(1).ALON;OBS(1).ALAT;sitevec.ELA(1:3:end)';sitevec.ELA(2:3:end)';sitevec.ELA(3:3:end)'];

oDIR=[DIR,'/vector'];
exid=exist(oDIR);
if exid~=7; mkdir(oDIR); end
FIDgrd=fopen([oDIR,'/CAL_vector_ela_',NAME,'_grid.txt'],'w');
FIDsite=fopen([oDIR,'/CAL_vector_ela_',NAME,'_site.txt'],'w');
fprintf(FIDgrd,'%f %f %f %f %f\n',CALVECgrd.ela);
fprintf(FIDsite,'%f %f %f %f %f\n',CALVECsite.ela);
fclose(FIDgrd);
fclose(FIDsite);
% 
end
%% Show results for makeing FIGURES
function ExportCoupling(DIR,TCHA,BLK,SDR)
NN=1;
folder=[DIR,'/coupling'];
exid=exist(folder);
if exid~=7; mkdir(folder); end
NFLT=size(TCHA.AVEFLT,1)/TCHA.NReplica;
TCHA.AVEFLT=reshape(TCHA.AVEFLT,NFLT,TCHA.NReplica);
TCHA.MEDFLT=reshape(TCHA.MEDFLT,NFLT,TCHA.NReplica);
TCHA.STDFLT=reshape(TCHA.STDFLT,NFLT,TCHA.NReplica);
for REP=1:TCHA.NReplica
  subfolder=[folder,'/replica',num2str(REP)];
  exid=exist(subfolder);
  if exid~=7; mkdir(subfolder); end
  FIDstdinfo=fopen([folder,'/Std_info.txt'],'w');
  fprintf(FIDstdinfo,'NB1 NB2 STDmax STDmin\n');
  for NB1=1:BLK(1).NBlock
    for NB2=NB1+1:BLK(1).NBlock
      NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
      if NF~=0
        FIDmain = fopen([subfolder,'/C_',num2str(NB1),'_',num2str(NB2),'.txt'],'w');
        FLTNUM = NN:NN+NF-1;
        AVECP = TCHA.AVEFLT(FLTNUM,REP);
        MEDCP = TCHA.MEDFLT(FLTNUM,REP);
        SDRs  =    SDR.flax(FLTNUM,REP);
        STD   = TCHA.STDFLT(FLTNUM,REP);
        clon = mean(BLK(1).BOUND(NB1,NB2).blon,2);
        clat = mean(BLK(1).BOUND(NB1,NB2).blat,2);
        cdep = mean(BLK(1).BOUND(NB1,NB2).bdep,2);
        outdata = [FLTNUM' ...
                   BLK(1).BOUND(NB1,NB2).blon ...
                   BLK(1).BOUND(NB1,NB2).blat ...
                   BLK(1).BOUND(NB1,NB2).bdep ...
                   clon clat cdep ...
                   AVECP MEDCP SDRs STD];
        fprintf(FIDmain,'# Contents\n');
        fprintf(FIDmain,'# TRI_No. Lon1 Lon2 Lon3 Lat1 Lat2 Lat3 C_Lon C_Lat C_Dep Mean_Coupling Median_Coupling SDR sigma\n');
        fprintf(FIDmain,'%5d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f\n',outdata');
        fprintf(FIDstdinfo,'%d %d %f %f\n',NB1,NB2,min(STD),max(STD));
        fclose(FIDmain);
        NN=NN+NF;
      end
    end
  end
end
fclose(FIDstdinfo);
end
%% Export coupling and SDR for each boundary
function ExportCouplingPair(DIR,BLK,TCHA,sdr,Coupling_Pair)
PAIR=Coupling_Pair.pair;
name=Coupling_Pair.NAME{1};
NN=1;
folder=[DIR,'/coupling'];
exid=exist(folder);
if exid~=7; mkdir(folder); end
for REP=1:TCHA.NReplica
  subfolder=[folder,'/replica',num2str(REP)];
  exid=exist(subfolder);
  if exid~=7; mkdir(subfolder); end
  FID = fopen([subfolder,'/CouplingTrace_',name,'.txt'],'w');
  fprintf(FID,'# Contents');
  fprintf(FID,'# FLT_No. Lon1 Lon2 Lat1 Lat2 C_Lon C_Lat Mean_Coupling Median_Coupling SDR');
  for NB1=1:BLK(1).NBlock
    for NB2=NB1+1:BLK(1).NBlock
      NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
      if NF~=0
        FLAG=0;
        for PN=1:size(PAIR,1)
          PAIRID=ismember([NB1 NB2],PAIR(PN,:));
          ISPAIR=sum(PAIRID);
          if ISPAIR==2;FLAG=1;break;end
        end
        FLTNUM=NN:NN+NF-1;
        if FLAG==1
          ID=BLK(1).BOUND(NB1,NB2).bdep==0;
          acID=sum(BLK(1).BOUND(NB1,NB2).bdep==0,2)==2;
          rmID=find(sum(BLK(1).BOUND(NB1,NB2).bdep~=0,2)==2);
          ID(rmID,:)=false;ID=ID';
          tmplon=BLK(1).BOUND(NB1,NB2).blon';
          tmplat=BLK(1).BOUND(NB1,NB2).blat';
          expLON=tmplon(ID);expLON=reshape(expLON',2,length(expLON)/2);
          expLAT=tmplat(ID);expLAT=reshape(expLAT',2,length(expLAT)/2);
          meanLON=mean(expLON,1);
          meanLAT=mean(expLAT,1);
          AVECP=TCHA.AVEFLT(FLTNUM,:);
          AVECP=AVECP(acID);
          MEDCP=TCHA.MEDFLT(FLTNUM,:);
          MEDCP=MEDCP(acID);
          SDR=sdr.flax(FLTNUM,:);
          SDR=SDR(acID);
          FLTNUM=FLTNUM(acID);
          outdata=[FLTNUM ...
                   expLON expLAT ... 
                   meanLON meanLAT ...
                   AVECP MEDCP SDR];
          fprintf(FID,'%5d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f \n',outdata');
        end
        NN=NN+NF;
      end
    end
  end
  fclose(FID);
end
end
%% Export strain rates of internal deformation.
function ExportInternalDeformation(DIR,TCHA,BLK)
NN=1;
folder=[DIR,'/rigid'];
exid=exist(folder);
if exid~=7; mkdir(folder); end
NINE=size(TCHA.AVEINE,1)/TCHA.NReplica;
TCHA.AVEINE=reshape(TCHA.AVEINE,NINE,TCHA.NReplica);
for REP=1:TCHA.NReplica
  subfolder=[folder,'/replica',num2str(REP)'];
  exid=exist(subfolder);
  if exid~=7; mkdir(subfolder); end
  FIDinternal=fopen([folder,'/Internal_Deformation.txt'],'w');
  fprintf(FIDinternal,'Block Latitude Longitude exx exy eyy emax emin thetaP shearMAX sig_exx sig_exy sig_eyy sig_emax sig_emin sig_shearMAX [nanostrain/yr] \n');
  for NB=1:BLK(1).NBlock
    exx=TCHA.AVEINE(3*NB-2,REP);
    exy=TCHA.AVEINE(3*NB-1,REP);
    eyy=TCHA.AVEINE(3*NB  ,REP);
    sigexx=TCHA.STDINE(3*NB-2,REP);
    sigexy=TCHA.STDINE(3*NB-1,REP);
    sigeyy=TCHA.STDINE(3*NB  ,REP);
    E=[exx exy;...
        exy eyy];
    [eigV,eigD]=eig(E);
    e1=eigD(1,1); e2=eigD(2,2);
    v1=eigV(:,1); v2=eigV(:,2);
    if e1>=e2
      emax=e1; axmax=v1;
      emin=e2; axmin=v2;
    else
      emax=e2; axmax=v2;
      emin=e1; axmin=v1;
    end
    thetaP=rad2deg(atan2(axmax(2),axmax(1)));
    if thetaP<0; thetaP=thetaP+360; end
    shearMAX=sqrt((1/4)*(exx-eyy)^2+exy^2);
    sigemax    =sqrt( ( 0.5 + 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigexx^2 ...
                     +(          1*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *  exy       )^2 *sigexy^2 ...
                     +( 0.5 - 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigeyy^2 );
    sigemin    =sqrt( ( 0.5 - 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigexx^2 ...
                     +(         -1*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *  exy       )^2 *sigexy^2 ...
                     +( 0.5 + 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigeyy^2 );
    sigshearMAX=sqrt( (       0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigexx^2 ...
                     +(          1*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *  exy       )^2 *sigexy^2 ...
                     +(      -0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigeyy^2 );
    fprintf(FIDinternal,'%2d %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',...
            NB,BLK(NB).LATinter,BLK(NB).LONinter,exx*1e9,exy*1e9,eyy*1e9,emax*1e9,emin*1e9,thetaP,...
            shearMAX*1e9,sigexx*1e9,sigexy*1e9,sigeyy*1e9,sigemax*1e9,sigemin*1e9,sigshearMAX*1e9);
  end
end
fclose(FIDinternal);
end
%% Translate coupling to slip deficit rate.
function [SDR]=coupling2sdr(TCHA,D,DPT,G)

Mp.SMP=TCHA.AVEPOL;
Mc.N=size(TCHA.MEDFLT,1)/3;
Mc.SMPMAT=reshape(...
    repmat(TCHA.MEDFLT,3*D.CNT,1)...
    ,3*Mc.N,TCHA.NReplica*D.CNT);
Mc.SMPMAT=reshape(Mc.SMPMAT(DPT.MID),3*Mc.N,TCHA.NReplica);
SDR.tmp=(G.TB*Mp.SMP).*DPT.CFINV.*Mc.SMPMAT;
SFLTNUM=sum(D.MID)./3;
HH=0;
H1=0;
SDR.flax=[];
TCHA.MEDFLT=reshape(TCHA.MEDFLT,Mc.N,TCHA.NReplica);
for ii=1:length(SFLTNUM)
  sdr.all=SDR.tmp(HH+1:HH+3*SFLTNUM(ii),:);
  sdr.str=sdr.all(1:SFLTNUM(ii),:);
  sdr.dip=sdr.all(SFLTNUM(ii)+1:SFLTNUM(ii)*2,:);
  sdr.tns=sdr.all(2*SFLTNUM(ii)+1:SFLTNUM(ii)*3,:);
  sdr.scaler=sqrt(sdr.str.^2+sdr.dip.^2+sdr.tns.^2);
  sdr.cp=TCHA.MEDFLT(H1+1:H1+SFLTNUM(ii),:);
  minus=sign(sdr.cp);
  sdr.scaler=minus.*sdr.scaler;
  SDR.flax=[SDR.flax; sdr.scaler];
  HH=HH+3*SFLTNUM(ii);
  H1=H1+SFLTNUM(ii);
end
clear SDR.tmp

end
%====================================================
function [lat,lon,ang]=xyzp2lla(X,Y,Z)
% XYZP2LLA  Converts Shpear coordinates from cartesian. Vectorized.
% GRS80
% CODE BY T.ITO 2017/03/11     ver0.1
% lat: deg, lon: deg, ang: deg/m.y.
lat=atan2(Z,sqrt(X.*X+Y.*Y)).*180/pi;
lon=atan2(Y,X).*180/pi;
ang=sqrt(X.*X+Y.*Y+Z.*Z).*(1e6.*(180./pi));
end
%====================================================
function [a,b,c,d,e,f]=out_cov(COVMAT)
COVPOL=COVMAT.*1e12.*1e8;
a=COVPOL(1,1);
b=COVPOL(1,2);
c=COVPOL(1,3);
d=COVPOL(2,2);
e=COVPOL(2,3);
f=COVPOL(3,3);
end
%%
function WRITE_VECTOR(TCHA,DIR,OBS,BLK,calvec,resvec,grdvec,GRD)
folder=[DIR,'/vector'];
exid=exist(folder);
if exid~=7; mkdir(folder); end
% 
OBSVEC=[OBS(1).ALON;OBS(1).ALAT;OBS(1).EVEC;OBS(1).NVEC;OBS(1).HVEC];
NR=0;
for REP=1:TCHA.NReplica
  subfolder=[folder,'/replica',num2str(REP)];
  exid=exist(subfolder);
  if exid~=7; mkdir(subfolder); end
  RESID=[REP,REP+TCHA.NReplica,REP+2*TCHA.NReplica];
  CALVEC.RIG=[OBS(1).ALON;OBS(1).ALAT;calvec.RIG(1:3:end,REP)';calvec.RIG(2:3:end,REP)';calvec.RIG(3:3:end,REP)'];
  CALVEC.ELA=[OBS(1).ALON;OBS(1).ALAT;calvec.ELA(1:3:end,REP)';calvec.ELA(2:3:end,REP)';calvec.ELA(3:3:end,REP)'];
  CALVEC.SUM=[OBS(1).ALON;OBS(1).ALAT;calvec.SUM(1:3:end,REP)';calvec.SUM(2:3:end,REP)';calvec.SUM(3:3:end,REP)'];
  RESVEC.SUM=[OBS(1).ALON;OBS(1).ALAT;resvec.SUM(RESID,:)];
  CALVECgrid.rig=[GRD(1).ALON;GRD(1).ALAT;grdvec.RIG(1:3:end,REP)';grdvec.RIG(2:3:end,REP)';grdvec.RIG(3:3:end,REP)'];
  CALVECgrid.ela=[GRD(1).ALON;GRD(1).ALAT;grdvec.ELA(1:3:end,REP)';grdvec.ELA(2:3:end,REP)';grdvec.ELA(3:3:end,REP)'];
% 
  FID=fopen([subfolder,'/OBS_vector.txt'],'w');
  fprintf(FID,'%f %f %f %f %f\n',OBSVEC);
  fclose(FID);
  FID=fopen([subfolder,'/CAL_vector.txt'],'w');
  fprintf(FID,'%f %f %f %f %f\n',CALVEC.SUM);
  fclose(FID);
  FID=fopen([subfolder,'/CAL_vector_rig_site.txt'],'w');
  fprintf(FID,'%f %f %f %f %f\n',CALVEC.RIG);
  fclose(FID);
  FID=fopen([subfolder,'/CAL_vector_ela_site.txt'],'w');
  fprintf(FID,'%f %f %f %f %f\n',CALVEC.ELA);
  fclose(FID);
  FID=fopen([subfolder,'/CAL_vector_rig_grid.txt'],'w');
  fprintf(FID,'%f %f %f %f %f\n',CALVECgrid.rig);
  fclose(FID);
  FID=fopen([subfolder,'/CAL_vector_ela_grid.txt'],'w');
  fprintf(FID,'%f %f %f %f %f\n',CALVECgrid.ela);
  fclose(FID);
  FID=fopen([subfolder,'/RES_vector.txt'],'w');
  fprintf(FID,'%f %f %f %f %f\n',RESVEC.SUM);
  fclose(FID);
  FID=fopen([subfolder,'/RMS_vector_block.txt'],'w');
  fprintf(FID,'BLOCK  RMS(mm/yr)\n');
  for ii=1:BLK(1).NBlock
    fprintf(FID,'%d %f\n',ii,resvec.RMS(ii));
  end
fclose(FID);
end
% 
end
%%
function [BO]=READ_local_bound(INPUT)
FID=fopen(INPUT.BOlocal1,'r');
data=fscanf(FID,'%f %f\n',[2 Inf]);
BO(1).lon=data(1,:);
BO(1).lat=data(2,:);
fclose(FID);
% 
FID=fopen(INPUT.BOlocal2,'r');
data=fscanf(FID,'%f %f\n',[2 Inf]);
BO(2).lon=data(1,:);
BO(2).lat=data(2,:);
fclose(FID);
% 
FID=fopen(INPUT.BOlocal3,'r');
data=fscanf(FID,'%f %f\n',[2 Inf]);
BO(3).lon=data(1,:);
BO(3).lat=data(2,:);
fclose(FID);
end
%%
function reslocal=residual_vector_local(BO,RESVEC)
BOID=inpolygon(RESVEC(1,:),RESVEC(2,:),BO.lon,BO.lat);
reslocal=RESVEC(:,BOID);
end
%%
function RESVEC=calc_residual_vector(OBS,BLK,TCHA,calvec)
obsV=[repmat(OBS(1).EVEC,TCHA.NReplica,1);...
      repmat(OBS(1).NVEC,TCHA.NReplica,1);...
      repmat(OBS(1).HVEC,TCHA.NReplica,1)];
calV=[calvec.SUM(1:3:end,:)';...
      calvec.SUM(2:3:end,:)';...
      calvec.SUM(3:3:end,:)'];
RESVEC.SUM=obsV-calV;

RESVEC.RMS=zeros(BLK(1).NBlock,TCHA.NReplica);
for REP=1:TCHA.NReplica
  for ii=1:BLK(1).NBlock
    ID=inpolygon(OBS(1).ALON,OBS(1).ALAT,BLK(ii).LON,BLK(ii).LAT);
    REPID=[REP,REP+TCHA.NReplica,REP+2*TCHA.NReplica];
    RESTMP=RESVEC.SUM(REPID,ID);
    RESVEC.RMS(ii,REP)=sqrt(sum(sum(RESTMP.^2))/sum(ID));
  end
end
end
%% 
function [GRDvec,GRD]=calc_vector_atmesh(BLK,TCHA,D,DPT,G,GRD,TRIg)
% minlon=120; maxlon=150;
% minlat=20 ; maxlat=50;
% interval= 0.4;
% [Xm,Ym]=meshgrid(minlon:interval:maxlon,minlat:interval:maxlat);
% XM=Xm(:);
% YM=Ym(:);
% GRD(1).ALON=XM';
% GRD(1).ALAT=YM';
% GRD(1).AHIG=zeros(size(GRD(1).ALON));
[Dg,Gg,GRD]=ReshapeGreen(BLK,GRD,TRIg);
% 
Mp.SMP=TCHA.AVEPOL;
Mc.N=size(TCHA.MEDFLT,1)/3;
Mc.SMPMAT=reshape(...
    repmat(TCHA.MEDFLT,3*D.CNT,1)...
    ,3*Mc.N,TCHA.NReplica*D.CNT);
Mc.SMPMAT=reshape(Mc.SMPMAT(DPT.MID),3*Mc.N,TCHA.NReplica);
% CALC APRIORI AND RESIDUAL COUPLING RATE SECTION
GRDvec.RIG=Gg.P*Mp.SMP;
GRDvec.ELA=Gg.C*((G.TB*Mp.SMP).*DPT.CFINV.*Mc.SMPMAT);
GRDvec.SUM=GRDvec.RIG+GRDvec.ELA;
% 
end
%% CALCULATE VECTOR BASED ON SAMPLED PARAMETER AND GREEN FUNCTION
function CALvec=calc_sampling_vector(OBS,BLK,TCHA,D,DPT,G)
% 
Mp.SMP=TCHA.AVEPOL;
Mc.N=size(TCHA.MEDFLT,1)/3;
Mc.SMPMAT=reshape(...
    repmat(TCHA.MEDFLT,3*D.CNT,1)...
    ,3*Mc.N,TCHA.NReplica*D.CNT);
Mc.SMPMAT=reshape(Mc.SMPMAT(DPT.MID),3*Mc.N,TCHA.NReplica);
% CALC APRIORI AND RESIDUAL COUPLING RATE SECTION
CALvec.RIG=G.P*Mp.SMP;
CALvec.ELA=G.C*((G.TB*Mp.SMP).*D(1).CFINV.*Mc.SMPMAT);
CALvec.SUM=CALvec.RIG+CALvec.ELA;

end
%% MAKE PART GREEN FUNCTION v2
function [TRIg,Dg,GRD]=MAKE_PART_GREEN(BLK,Grid_Setting)
% 
minlon=Grid_Setting(1);
maxlon=Grid_Setting(2);
minlat=Grid_Setting(3);
maxlat=Grid_Setting(4);
interval=Grid_Setting(5);
[Xm,Ym]=meshgrid(minlon:interval:maxlon,minlat:interval:maxlat);
XM=Xm(:);
YM=Ym(:);
GRD(1).ALON=XM';
GRD(1).ALAT=YM';
GRD(1).AHIG=zeros(size(GRD(1).ALON));
% 
N=0;
GRD(1).NGRD=length(GRD(1).ALON);
for N=1:GRD(1).NGRD
  GRD(1).AXYZ(N,:)=conv2ell_hig(GRD(1).ALAT(N),GRD(1).ALON(N),GRD(1).AHIG(N));
end
GRD(1).ABLK=zeros(GRD(1).NGRD,1);
% 
for N=1:BLK(1).NBlock
  IND=inpolygon(GRD(1).ALON,GRD(1).ALAT,BLK(N).LON,BLK(N).LAT);
  GRD(1).ABLK(IND)=N;
  GRD(N).NBLK=sum(IND);
  GRD(N).LAT=GRD(1).ALAT(IND);
  GRD(N).LON=GRD(1).ALON(IND);
  GRD(N).HIG=GRD(1).AHIG(IND);
  GRD(N).OXYZ=conv2ell_hig(GRD(N).LAT,GRD(N).LON,GRD(N).HIG);
end
Dg(1).IND=find(GRD(1).ABLK~=0)';
GRD(1).ALON=GRD(1).ALON(:,Dg(1).IND);
GRD(1).ALAT=GRD(1).ALAT(:,Dg(1).IND);
GRD(1).AHIG=GRD(1).AHIG(:,Dg(1).IND);
GRD(1).ABLK=GRD(1).ABLK(Dg(1).IND,:);
GRD(1).AXYZ=GRD(1).AXYZ(Dg(1).IND,:);
% NGRD=size(GRD(1).ALON,2);
% 
[TRIg]=GREEN_TRI(BLK,GRD);
% 
end
%% MAKE PARTS OF GREEN FUNCTION AT MESH GRID
function [Dg,Gg,GRD]=ReshapeGreen(BLK,GRD,TRIg)
% 
NGRD=size(GRD(1).ALON,2);
TMP.P=zeros(3*NGRD,3.*BLK(1).NBlock);
MC=1;
MT=1;
MR=1;
for NB1=1:BLK(1).NBlock
  Dg(1).grdid=zeros(1,NGRD);
  Dg(1).grdid(1,GRD(1).ABLK==NB1)=true;
  Dg(1).GRDID(:,NB1)=reshape(repmat(Dg(1).grdid,3,1),3*NGRD,1);
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(TRIg(1).BOUND(NB1,NB2).clon,2);
    if NF~=0
      TMP.C(1:3*NGRD,MC     :MC+  NF-1)=TRIg(1).BOUND(NB1,NB2).GSTR;
      TMP.C(1:3*NGRD,MC+  NF:MC+2*NF-1)=TRIg(1).BOUND(NB1,NB2).GDIP;
      TMP.C(1:3*NGRD,MC+2*NF:MC+3*NF-1)=TRIg(1).BOUND(NB1,NB2).GTNS;
      MC=MC+3*NF;
      MT=MT+2*NF;
      MR=MR+  NF;
    end
  end
%   
  IND=GRD(1).ABLK==NB1;  
  NIND=[zeros(size(IND)),IND,zeros(size(IND))]; NIND=logical(reshape(NIND',3*NGRD,1));
  EIND=[IND,zeros(size(IND)),zeros(size(IND))]; EIND=logical(reshape(EIND',3*NGRD,1));
  TMP.P(EIND,3*NB1-2)=-GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,3);
  TMP.P(EIND,3*NB1-1)=-GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,3);
  TMP.P(EIND,3*NB1  )= GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,2)                    +GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,1);
  TMP.P(NIND,3*NB1-2)= GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,3)+GRD(1).AXYZ(IND,6).*GRD(1).AXYZ(IND,2);
  TMP.P(NIND,3*NB1-1)=-GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,3)-GRD(1).AXYZ(IND,6).*GRD(1).AXYZ(IND,1);
  TMP.P(NIND,3*NB1  )= GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,2)-GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,1);
end
% 
Gg(1).C  =TMP.C;
Gg(1).P  =TMP.P;

end
%% MAKE GREEN FUNCTION
function [TRI]=GREEN_TRI(BLK,GRD)
% Coded by Takeo Ito 2017/01/02 (ver 1.1)
PR=0.25;
ND=size(GRD(1).ALAT,2);
%
ALAT=mean(GRD(1).ALAT(:));
ALON=mean(GRD(1).ALON(:));
[GRDx,GRDy]=PLTXY(GRD(1).ALAT,GRD(1).ALON,ALAT,ALON);
GRDz=-1e-3.*GRD(1).AHIG;
%
TRI(1).TNF=0;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
    TRI(1).BOUND(NB1,NB2).clat=[];
    TRI(1).BOUND(NB1,NB2).clon=[];
    TRI(1).BOUND(NB1,NB2).cdep=[];
    if NF~=0
      TRI(1).BOUND(NB1,NB2).GSTR=zeros(3*ND,NF);
      TRI(1).BOUND(NB1,NB2).GDIP=zeros(3*ND,NF);
      TRI(1).BOUND(NB1,NB2).GTNS=zeros(3*ND,NF);
%
      fprintf('==================\n Block %2i : Block %2i \n Number of TRI sub-faults : %4i \n',NB1,NB2,NF)
%
      for N=1:NF
        [TRIx,TRIy]=PLTXY(BLK(1).BOUND(NB1,NB2).blat(N,:),BLK(1).BOUND(NB1,NB2).blon(N,:),ALAT,ALON);
        TRIz=-1.*BLK(1).BOUND(NB1,NB2).bdep(N,:);
        F_LOC=[TRIx;TRIy;TRIz];
        [F,DA,NV,ST,DP]=EST_FAULT_TRI(F_LOC);
        NV(3)=-NV(3);DP(3)=-DP(3);
        TRI(1).BOUND(NB1,NB2).clat(N)=mean(BLK(1).BOUND(NB1,NB2).blat(N,:));
        TRI(1).BOUND(NB1,NB2).clon(N)=mean(BLK(1).BOUND(NB1,NB2).blon(N,:));
        TRI(1).BOUND(NB1,NB2).cdep(N)=mean(BLK(1).BOUND(NB1,NB2).bdep(N,:)); % up is plus
        TRI(1).BOUND(NB1,NB2).DA(N)=DA;
        TRI(1).BOUND(NB1,NB2).NV(N,:)=NV;
        TRI(1).BOUND(NB1,NB2).ST(N,:)=ST;
        TRI(1).BOUND(NB1,NB2).DP(N,:)=DP;
        TRI(1).BOUND(NB1,NB2).OXYZ(N,:)=conv2ell_hig(TRI(1).BOUND(NB1,NB2).clat(N),TRI(1).BOUND(NB1,NB2).clon(N),1e3.*TRI(1).BOUND(NB1,NB2).cdep(N));
        U=CalcTriDisps(GRDx,GRDy,GRDz,TRIx,TRIy,TRIz,PR,1,0,0);
        TRI(1).BOUND(NB1,NB2).GSTR(1:3:3*ND,N)=U.x; %E
        TRI(1).BOUND(NB1,NB2).GSTR(2:3:3*ND,N)=U.y; %N
        TRI(1).BOUND(NB1,NB2).GSTR(3:3:3*ND,N)=-U.z; %D
        U=CalcTriDisps(GRDx,GRDy,GRDz,TRIx,TRIy,TRIz,PR,0,1,0);
        TRI(1).BOUND(NB1,NB2).GTNS(1:3:3*ND,N)=U.x; %E
        TRI(1).BOUND(NB1,NB2).GTNS(2:3:3*ND,N)=U.y; %N
        TRI(1).BOUND(NB1,NB2).GTNS(3:3:3*ND,N)=-U.z; %D 
        U=CalcTriDisps(GRDx,GRDy,GRDz,TRIx,TRIy,TRIz,PR,0,0,1);
        TRI(1).BOUND(NB1,NB2).GDIP(1:3:3*ND,N)=U.x; %E
        TRI(1).BOUND(NB1,NB2).GDIP(2:3:3*ND,N)=U.y; %N
        TRI(1).BOUND(NB1,NB2).GDIP(3:3:3*ND,N)=-U.z; %D
        if mod(N,ceil(NF/3)) == 1
          fprintf('MAKE GREEN at TRI sub-faults : %4i / %4i \n',N,NF)
        end
      end
      TRI(1).TNF=TRI(1).TNF+NF;
    end
  end
end
disp('==================')
disp('PASS GREEN_TRI')
disp('==================')
end
%% ESTIMATE FAULT PARAMETERS FOR TRI
function [FLOC,DA,NV,ST,DP]=EST_FAULT_TRI(loc_f)
% Coded by Takeo Ito 2015/11/11 (ver 1.0)
[DA]=AREA_TRI(loc_f(1,:),loc_f(2,:),loc_f(3,:));
FLOC=mean(loc_f,2)';
[NV,ST,DP]=EST_STRDIP_TRI(loc_f(1,:),loc_f(2,:),loc_f(3,:));
end
%% ESTIMATE AREA AT SUB-FAULT FOR TRI
function [DA]=AREA_TRI(X,Y,Z)
% CALC. AREA IN THREE DIMENSION USING HERON'S FOMULA
% Coded by Takeo Ito 2006/03/04 (ver 1.0)
LENG(1)=sqrt((X(1)-X(2)).^2+(Y(1)-Y(2)).^2+(Z(1)-Z(2)).^2);
LENG(2)=sqrt((X(2)-X(3)).^2+(Y(2)-Y(3)).^2+(Z(2)-Z(3)).^2);
LENG(3)=sqrt((X(3)-X(1)).^2+(Y(3)-Y(1)).^2+(Z(3)-Z(1)).^2);
S1=(LENG(1)+LENG(2)+LENG(3))./2;
DA=sqrt(S1*(S1-LENG(1))*(S1-LENG(2))*(S1-LENG(3)));
end
%% ESTIMATE STRKE AND DIP FOR TRI FAULT
function [NV,ST,DP]=EST_STRDIP_TRI(X,Y,Z)
%==========
% CALC. STR AND DIP ON FAULT
% CODE BY T.ITO (2006/03/04)
% Modified by T.ITO (2015/11/13)
% Modified by T.ITO (2016/02/16)
% Modified by Kimura(2017/04/17)
% UP IS MINUS
%==========
NV=cross([X(2);Y(2);Z(2)]-[X(1);Y(1);Z(1)], [X(3);Y(3);Z(3)]-[X(1);Y(1);Z(1)]);
NV=NV./norm(NV);
if (NV(3) < 0); NV = -NV; end; % Enforce clockwise circulation
ST=[-sin(atan2(NV(2),NV(1))) cos(atan2(NV(2),NV(1))) 0];
DP=cross(NV,ST);
end
%% CONVERT TO XYZ FROM ELL
function [x,y,z]=ell2xyz(lat,lon,h)
% ELL2XYZ  Converts ellipsoidal coordinates to cartesian. Vectorized.
% GRS80
% CODE BY T.ITO 2006/12/13     ver0.1
% BUG FIX  BY T.ITO 2015/11/13 ver0.2
% 
a=6378137.0; % m
f=1./298.257222101;
e2=1-(1-f)^2;
%
rad=pi/180;
lat=lat.*rad;
lon=lon.*rad;
clat=cos(lat);
clon=cos(lon);
slat=sin(lat);
slon=sin(lon);
%
v=a./sqrt(1-e2.*slat.*slat);
x=(v+h).*clat.*clon;
y=(v+h).*clat.*slon;
z=(v.*(1-e2)+h).*slat;
end
%% CONVERT TO XYZ FROM ELL AT SURFACE
function [OOxyz]=conv2ell_hig(Olat,Olon,Ohig)
Olat=Olat(:);
Olon=Olon(:);
Ohig=Ohig(:);
deg2rad=pi/180;
[Oxyz(:,1),Oxyz(:,2),Oxyz(:,3)]=ell2xyz(Olat,Olon,Ohig);
Oxyz = Oxyz*1e3;
OOxyz=[Oxyz sin(Olat*deg2rad) sin(Olon*deg2rad) cos(Olat*deg2rad) cos(Olon*deg2rad)];
end
%% MAKE GREEN FUNCTION
%====================================================
function [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds)
% CalcTriDisps.m
%
% Calculates displacements due to slip on a triangular dislocation in an
% elastic half space utilizing the Comninou and Dunders (1975) expressions
% for the displacements due to an angular dislocation in an elastic half
% space.
%
% Arguments
%  sx : x-coordinates of observation points
%  sy : y-coordinates of observation points
%  sz : z-coordinates of observation points
%  x  : x-coordinates of triangle vertices.
%  y  : y-coordinates of triangle vertices.
%  z  : z-coordinates of triangle vertices.
%  pr : Poisson's ratio
%  ss : strike slip displacement
%  ts : tensile slip displacement
%  ds : dip slip displacement
%
% Returns
%  U  : structure containing the displacements (U.x, U.y, U.z)
%
% This paper should and related code should be cited as:
% Brendan J. Meade, Algorithms for the calculation of exact 
% displacements, strains, and stresses for Triangular Dislocation 
% Elements in a uniform elastic half space, Computers & 
% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
%
% Use at your own risk and please let me know of any bugs/errors!
%
% Copyright (c) 2006 Brendan Meade
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
x=x(:);sx=sx(:);
y=y(:);sy=sy(:);
z=z(:);sz=sz(:);

% Calculate the slip vector in XYZ coordinates
normVec                      = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
normVec                      = normVec./norm(normVec);
if (normVec(3) < 0) % Enforce clockwise circulation
   normVec                   = -normVec;
   [x(2),x(3)]               = swap(x(2), x(3));
   [y(2),y(3)]               = swap(y(2), y(3));
   [z(2),z(3)]               = swap(z(2), z(3));
end
strikeVec                    = [-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
dipVec                       = cross(normVec, strikeVec);
slipComp                     = [ss ds ts];
slipVec                      = [strikeVec(:) dipVec(:) normVec(:)] * slipComp(:);
% Solution vectors
U.x                          = zeros(size(sx));
U.y                          = zeros(size(sx));
U.z                          = zeros(size(sx));
% Add a copy of the first vertex to the vertex list for indexing
x(4)                         = x(1);
y(4)                         = y(1);
z(4)                         = z(1);
%
for iTri = 1:3
   % Calculate strike and dip of current leg
   strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
   [rx,~]                   = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
   dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));
%   
   if dip >= 0
      beta                  = pi/180*(90-dip);
      if beta > pi/2
         beta               = pi/2-beta;
      end
   else
      beta                  = -pi/180*(90+dip);
      if beta < -pi/2
         beta = pi/2-abs(beta);
      end
   end
%
   ssVec                    = [ cos(strike/180*pi) sin(strike/180*pi) 0];
   tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
   dsVec                    = cross(ssVec, tsVec);
   lss                      = dot(slipVec, ssVec);
   lts                      = dot(slipVec, tsVec);
   lds                      = dot(slipVec, dsVec);
%
   if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
      % First angular dislocation
      [sx1,sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
      [ux1,uy1,uz1]             = adv(sx1, sy1, sz-z(iTri  ), z(iTri), beta, pr, lss, lts, lds);
                                   
      % Second angular dislocation
      [sx2,sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike); 
      [ux2,uy2,uz2]             = adv(sx2, sy2, sz-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);

      % Rotate vectors to correct for strike
      [uxn,uyn]                 = RotateXyVec(ux1-ux2, uy1-uy2, strike);
      uzn                       = uz1-uz2;
 
      % Add the displacements from current leg
      U.x                       = U.x + uxn;
      U.y                       = U.y + uyn;
      U.z                       = U.z + uzn;
   end
end

% Identify indices for stations under current triangle
inPolyIdx                       = find(inpolygon(sx, sy, x, y) == 1);
underIdx = [];
for iIdx = 1 : numel(inPolyIdx)
   d                            = LinePlaneIntersect(x, y, z, sx(inPolyIdx(iIdx)), sy(inPolyIdx(iIdx)), sz(inPolyIdx(iIdx)));
   if d(3)-sz(inPolyIdx(iIdx)) < 0
      underIdx = [underIdx ; inPolyIdx(iIdx)];
   end
end
% Apply static offset to the points that lie underneath the current triangle
U.x(underIdx)                = U.x(underIdx) - slipVec(1);
U.y(underIdx)                = U.y(underIdx) - slipVec(2);
U.z(underIdx)                = U.z(underIdx) - slipVec(3);
%
end
%====================================================
function d = LinePlaneIntersect(x, y, z, sx, sy, sz)
% Calculate the intersection of a line and a plane using a parametric
% representation of the plane.  This is hardcoded for a vertical line.
numerator                       = [1 1 1 1 ; x(1) x(2) x(3) sx ; y(1) y(2) y(3) sy ; z(1) z(2) z(3) sz];
numerator                       = det(numerator);
denominator                     = [1 1 1 0 ; x(1) x(2) x(3) 0 ; y(1) y(2) y(3) 0 ; z(1) z(2) z(3) -sz];
denominator                     = det(denominator);
if denominator == 0;
   denominator                  = eps;
end
t                               = numerator/denominator; % parametric curve parameter
d                               = [sx sy sz]-([sx sy 0]-[sx sy sz])*t;
end
%====================================================
function [a,b] = swap(a, b)
% Swap two values
temp                            = a;
a                               = b;
b                               = temp;
end
%====================================================
function [xp,yp] = RotateXyVec(x, y, alpha)
% Rotate a vector by an angle alpha
x                             = x(:);
y                             = y(:);
alpha                         = pi/180*alpha;
xp                            = cos(alpha).*x - sin(alpha).*y;
yp                            = sin(alpha).*x + cos(alpha).*y;
end
%====================================================
function [v1,v2,v3] = adv(y1, y2, y3, a, beta, nu, B1, B2, B3)
% These are the displacements in a uniform elastic half space due to slip
% on an angular dislocation (Comninou and Dunders, 1975).  Some of the
% equations for the B2 and B3 cases have been corrected following Thomas
% 1993.  The equations are coded in way such that they roughly correspond
% to each line in original text.  Exceptions have been made where it made 
% more sense because of grouping symbols.

sinbeta           = sin(beta);
cosbeta           = cos(beta);
cotbeta           = cot(beta);
z1                = y1.*cosbeta - y3.*sinbeta;
z3                = y1.*sinbeta + y3.*cosbeta;
R2                = y1.*y1 + y2.*y2 + y3.*y3;
R                 = sqrt(R2);
y3bar             = y3 + 2.*a;
z1bar             = y1.*cosbeta + y3bar.*sinbeta;
z3bar             = -y1.*sinbeta + y3bar.*cosbeta;
R2bar             = y1.*y1 + y2.*y2 + y3bar.*y3bar;
Rbar              = sqrt(R2bar);
F                 = -atan2(y2, y1) + atan2(y2, z1) + atan2(y2.*R.*sinbeta, y1.*z1+(y2.*y2).*cosbeta);
Fbar              = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2.*Rbar.*sinbeta, y1.*z1bar+(y2.*y2).*cosbeta);

% Case I: Burgers vector (B1,0,0)
v1InfB1           = 2.*(1-nu).*(F+Fbar) - y1.*y2.*(1./(R.*(R-y3)) + 1./(Rbar.*(Rbar+y3bar))) - ...
                    y2.*cosbeta.*((R.*sinbeta-y1)./(R.*(R-z3)) + (Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB1           = (1-2.*nu).*(log(R-y3)+log(Rbar+y3bar) - cosbeta.*(log(R-z3)+log(Rbar+z3bar))) - ...
                    y2.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar)) - cosbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar))));
v3InfB1           = y2 .* (1./R - 1./Rbar - cosbeta.*((R.*cosbeta-y3)./(R.*(R-z3)) - (Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1InfB1           = v1InfB1 ./ (8.*pi.*(1-nu));
v2InfB1           = v2InfB1 ./ (8.*pi.*(1-nu));
v3InfB1           = v3InfB1 ./ (8.*pi.*(1-nu));

v1CB1             = -2.*(1-nu).*(1-2.*nu).*Fbar.*(cotbeta.*cotbeta) + (1-2.*nu).*y2./(Rbar+y3bar) .* ((1-2.*nu-a./Rbar).*cotbeta - y1./(Rbar+y3bar).*(nu+a./Rbar)) + ...
                    (1-2.*nu).*y2.*cosbeta.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*(-(1-2.*nu).*cotbeta + y1./(Rbar+y3bar) .* (2.*nu+a./Rbar) + a.*y1./(Rbar.*Rbar)) + ...
                    y2.*(y3bar-a)./(Rbar.*(Rbar+z3bar)).*(cosbeta./(Rbar+z3bar).*((Rbar.*cosbeta+y3bar) .* ((1-2.*nu).*cosbeta-a./Rbar).*cotbeta + 2.*(1-nu).*(Rbar.*sinbeta-y1).*cosbeta) - a.*y3bar.*cosbeta.*cotbeta./(Rbar.*Rbar));
v2CB1             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)-nu).*log(Rbar+y3bar) -(2.*(1-nu).*(cotbeta.*cotbeta)+1-2.*nu).*cosbeta.*log(Rbar+z3bar)) - ...
                    (1-2.*nu)./(Rbar+y3bar).*(y1.*cotbeta.*(1-2.*nu-a./Rbar) + nu.*y3bar - a + (y2.*y2)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
                    (1-2.*nu).*z1bar.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    (y3bar-a)./(Rbar+y3bar).*(-2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta-a) + (y2.*y2)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar)+a.*(y2.*y2)./(Rbar.*Rbar.*Rbar)) + ...
                    (y3bar-a)./(Rbar+z3bar).*((cosbeta.*cosbeta) - 1./Rbar.*((1-2.*nu).*z1bar.*cotbeta+a.*cosbeta) + a.*y3bar.*z1bar.*cotbeta./(Rbar.*Rbar.*Rbar) - 1./(Rbar.*(Rbar+z3bar)) .* ((y2.*y2).*(cosbeta.*cosbeta) - a.*z1bar.*cotbeta./Rbar.*(Rbar.*cosbeta+y3bar)));

v3CB1             = 2.*(1-nu).*(((1-2.*nu).*Fbar.*cotbeta) + (y2./(Rbar+y3bar).*(2.*nu+a./Rbar)) - (y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar))) + ...
                    y2.*(y3bar-a)./Rbar.*(2.*nu./(Rbar+y3bar)+a./(Rbar.*Rbar)) + ...
                    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*(1-2.*nu-(Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta + a./Rbar) - a.*y3bar./(Rbar.*Rbar));

v1CB1             = v1CB1 ./ (4.*pi.*(1-nu));
v2CB1             = v2CB1 ./ (4.*pi.*(1-nu));
v3CB1             = v3CB1 ./ (4.*pi.*(1-nu));

v1B1              = v1InfB1 + v1CB1;
v2B1              = v2InfB1 + v2CB1;
v3B1              = v3InfB1 + v3CB1;


% Case II: Burgers vector (0,B2,0)
v1InfB2           = -(1-2.*nu).*(log(R-y3) + log(Rbar+y3bar)-cosbeta.*(log(R-z3)+log(Rbar+z3bar))) + ...
                    y1.*y1.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) + z1.*(R.*sinbeta-y1)./(R.*(R-z3)) + z1bar.*(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar));
v2InfB2           = 2.*(1-nu).*(F+Fbar) + y1.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) - y2.*(z1./(R.*(R-z3))+z1bar./(Rbar.*(Rbar+z3bar)));
v3InfB2           = -(1-2.*nu).*sinbeta.*(log(R-z3)-log(Rbar+z3bar)) - y1.*(1./R-1./Rbar) + z1.*(R.*cosbeta-y3)./(R.*(R-z3)) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar));
v1InfB2           = v1InfB2 ./ (8.*pi.*(1-nu));
v2InfB2           = v2InfB2 ./ (8.*pi.*(1-nu));
v3InfB2           = v3InfB2 ./ (8.*pi.*(1-nu));

v1CB2             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)+nu).*log(Rbar+y3bar) - (2.*(1-nu).*(cotbeta.*cotbeta)+1).*cosbeta.*log(Rbar+z3bar)) + ...
                    (1-2.*nu)./(Rbar+y3bar).* (-(1-2.*nu).*y1.*cotbeta+nu.*y3bar-a+a.*y1.*cotbeta./Rbar + (y1.*y1)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
                    (1-2.*nu).*cotbeta./(Rbar+z3bar).*(z1bar.*cosbeta - a.*(Rbar.*sinbeta-y1)./(Rbar.*cosbeta)) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    (y3bar-a)./(Rbar+y3bar).*(2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta+a) - (y1.*y1)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar) - a.*(y1.*y1)./(Rbar.*Rbar.*Rbar)) + ...
                    (y3bar-a).*cotbeta./(Rbar+z3bar).*(-cosbeta.*sinbeta+a.*y1.*y3bar./(Rbar.*Rbar.*Rbar.*cosbeta) + (Rbar.*sinbeta-y1)./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta))));
v2CB2             = 2.*(1-nu).*(1-2.*nu).*Fbar.*cotbeta.*cotbeta + (1-2.*nu).*y2./(Rbar+y3bar).*(-(1-2.*nu-a./Rbar).*cotbeta + y1./(Rbar+y3bar).*(nu+a./Rbar)) - ...
                    (1-2.*nu).*y2.*cotbeta./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) - a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*((1-2.*nu).*cotbeta - 2.*nu.*y1./(Rbar+y3bar) - a.*y1./Rbar.*(1./Rbar+1./(Rbar+y3bar))) + ...
                    y2.*(y3bar-a).*cotbeta./(Rbar.*(Rbar+z3bar)).*(-2.*(1-nu).*cosbeta + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) + a.*y3bar./((Rbar.*Rbar).*cosbeta));
v3CB2             = -2.*(1-nu).*(1-2.*nu).*cotbeta .* (log(Rbar+y3bar)-cosbeta.*log(Rbar+z3bar)) - ...
                    2.*(1-nu).*y1./(Rbar+y3bar).*(2.*nu+a./Rbar) + 2.*(1-nu).*z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar) + ...
                   (y3bar-a)./Rbar.*((1-2.*nu).*cotbeta-2.*nu.*y1./(Rbar+y3bar)-a.*y1./(Rbar.*Rbar)) - ...
                   (y3bar-a)./(Rbar+z3bar).*(cosbeta.*sinbeta + (Rbar.*cosbeta+y3bar).*cotbeta./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar)) + a./Rbar.*(sinbeta - y3bar.*z1bar./(Rbar.*Rbar) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1CB2             = v1CB2 ./ (4.*pi.*(1-nu));
v2CB2             = v2CB2 ./ (4.*pi.*(1-nu));
v3CB2             = v3CB2 ./ (4.*pi.*(1-nu));

v1B2              = v1InfB2 + v1CB2;
v2B2              = v2InfB2 + v2CB2;
v3B2              = v3InfB2 + v3CB2;


% Case III: Burgers vector (0,0,B3)
v1InfB3           = y2.*sinbeta.*((R.*sinbeta-y1)./(R.*(R-z3))+(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB3           = (1-2.*nu).*sinbeta.*(log(R-z3)+log(Rbar+z3bar)) - (y2.*y2).*sinbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar)));
v3InfB3           = 2.*(1-nu).*(F-Fbar) + y2.*sinbeta.*((R.*cosbeta-y3)./(R.*(R-z3))-(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar)));
v1InfB3           = v1InfB3 ./ (8.*pi.*(1-nu));
v2InfB3           = v2InfB3 ./ (8.*pi.*(1-nu));
v3InfB3           = v3InfB3 ./ (8.*pi.*(1-nu));

v1CB3             = (1-2.*nu).*(y2./(Rbar+y3bar).*(1+a./Rbar) - y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar)) - ...
                    y2.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) + ...
                    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*((Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v2CB3             = (1-2.*nu).*(-sinbeta.*log(Rbar+z3bar) - y1./(Rbar+y3bar).*(1+a./Rbar) + z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar)) + ...
                    y1.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) - ...
                    (y3bar-a)./(Rbar+z3bar).*(sinbeta.*(cosbeta-a./Rbar) + z1bar./Rbar.*(1+a.*y3bar./(Rbar.*Rbar)) - ...
                    1./(Rbar.*(Rbar+z3bar)).*((y2.*y2).*cosbeta.*sinbeta - a.*z1bar./Rbar.*(Rbar.*cosbeta+y3bar)));
v3CB3             = 2.*(1-nu).*Fbar + 2.*(1-nu).*(y2.*sinbeta./(Rbar+z3bar).*(cosbeta + a./Rbar)) + ...
                    y2.*(y3bar-a).*sinbeta./(Rbar.*(Rbar+z3bar)).*(1 + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v1CB3             = v1CB3 ./ (4.*pi.*(1-nu));
v2CB3             = v2CB3 ./ (4.*pi.*(1-nu));
v3CB3             = v3CB3 ./ (4.*pi.*(1-nu));

v1B3              = v1InfB3 + v1CB3;
v2B3              = v2InfB3 + v2CB3;
v3B3              = v3InfB3 + v3CB3;


% Sum the for each slip component
v1                = B1.*v1B1 + B2.*v1B2 + B3.*v1B3;
v2                = B1.*v2B1 + B2.*v2B2 + B3.*v2B3;
v3                = B1.*v3B1 + B2.*v3B2 + B3.*v3B3;
end
%====================================================
function [V]=volcanic_deformation_mogi(OBS)
[VOL]=READ_VOLCANIC_CONFIGURE;
nu=0.25;           %Poisson's ratio
UX=zeros(1,OBS(1).NOBS);
UY=zeros(1,OBS(1).NOBS);
UZ=zeros(1,OBS(1).NOBS);
for VO=1:VOL.N
  [OBS(1).X,OBS(1).Y]=PLTXY(OBS(1).ALAT,OBS(1).ALON,VOL.lat(VO),VOL.lon(VO));
  OBS(1).X=1e3.*OBS(1).X;OBS(1).Y=1e3.*OBS(1).Y;
  [OBS(1).th,OBS(1).rho]=cart2pol(OBS(1).X,OBS(1).Y);
  [ur,uz]=mogi(OBS(1).rho,VOL.dep(VO),VOL.Vdef(VO),nu);
  [ux,uy]=pol2cart(OBS(1).th,ur);
  UX=UX+ux;UY=UY+uy;UZ=UZ+uz;
end
V=reshape(1e3.*[UX;UY;UZ],OBS(1).NOBS*3,1);

end
%% ====================================================
function [VOL]=READ_VOLCANIC_CONFIGURE
%Parameter estimated by Nishimura et al.[2007]
% Code	Latitude, deg	Longitude, deg	Depth, km	Volume Change Rate, 106 m3/yr	Resolution by Observed Data,%	Volcano
% VA	35.360	138.730	15.0	3.94 +/- 0.26	99	Mount Fuji
% VB	34.969	139.152	10.0	3.93 +/- 0.18	100	Eastern Izu
% VC	34.743	139.405	4.0	1.62 +/- 0.05	100	Oshima
% VD	34.335	139.228	4.9	3.20 +/- 0.15	100	Niijima
% VE	34.238	139.151	1.2	0.42 +/- 0.05	100	Kouzushima
% VF	34.064	139.513	9.5	6.16 +/- 0.25	99	Miyakejima
VOL.code=['VA';'VB';'VC';'VD';'VE';'VF'];
VOL.lat=[35.36;34.969;34.743;34.335;34.238;34.064];
VOL.lon=[138.73;139.152;139.405;139.228;139.151;139.513];
VOL.dep=1e3.*[15;10;4;4.9;1.2;9.5];
VOL.Vdef=1e6.*[3.94;3.93;1.62;3.20;0.42;6.16];
VOL.N=length(VOL.lon);
end
%% ====================================================
function [ur,uz,dt,er,et] = mogi(varargin)
warning('off','all')
%MOGI	Mogi's model (point source in elastic half-space).
%	[Ur,Uz,Dt,Er,Et] = MOGI(R,F,V,nu) or MOGI(R,F,A,P,E,nu) computes radial 
%	and vertical displacements Ur and Uz, ground tilt Dt, radial and 
%	tangential strain Er and Et on surface, at a radial distance R 
%	from the top of the source due to a hydrostatic pressure inside a 
%	sphere of radius A at depth F, in a homogeneous, semi-infinite elastic
%	body and approximation for A << F (center of dilatation).
%
%	MOGI(R,F,V) and MOGI(R,F,A,?�,P) are also allowed for compatibility 
%	(Mogi's original equation considers an isotropic material with Lam?s 
%	constants equal, i.e., lambda = ?�, Poisson's ratio = 0.25).
%
%	Input variables are:
%	   F: depth of the center of the sphere from the surface,
%	   V: volumetric change of the sphere,
%	   A: radius of the sphere,
%	   P: hydrostatic pressure change in the sphere,
%	   E: elasticity (Young's modulus),
%	  nu: Poisson's ratio,
%	   ?�: rigidity (Lam?s constant in case of isotropic material).
%
%	Notes:
%		- Equations are all vectorized, so variables R,F,V,A,?� and P are 
%		  scalar but any of them can be vector or matrix, then outputs 
%		  will be vector or matrix of the same size.
%		- Convention: Uz > 0 = UP, F is depth so in -Z direction.
%		- Units should be constistent, e.g.: R, F, A, Ur and Uz in m imply
%		  V in m3; E, ?� and P in Pa; Dt in rad, Er, Et and nu dimensionless.
%
%	Example for a 3-D plot of exagerated deformed surface due to a 1-bar
%	overpressure in a 10-cm radius sphere at 1-m depth in rock:
%	  [x,y] = meshgrid(-3:.1:3);
%	  [th,rho] = cart2pol(x,y);
%	  [ur,uz] = mogi(rho,1,0.1,1e5,10e9,0.25);
%	  [ux,uy] = pol2cart(th,ur);
%	  ps = 1e8;
%	  surf(x+ux*ps,y+uy*ps,uz*ps), axis equal, light

error(nargchk(3,6,nargin))

for i = 1:nargin
	if ~isnumeric(varargin{i})
		error('All input arguments must be numeric.')
	end
end

% to check if input arguments have compatible sizes, constructs a complex
% vector of sizes, then uses UNIQUE on variables that are not scalar
sz = complex(cellfun('size',varargin,1),cellfun('size',varargin,2));
if length(unique(sz(find(sz~=complex(1,1))))) > 1
	error('All inputs must be scalar or matrix of the same size.')
end

r = varargin{1};
f = varargin{2};

switch nargin
	case 3	% MOGI(R,F,V)
		v = varargin{3};
		nu = 0.25;
	case 4	% MOGI(R,F,V,nu)
		v = varargin{3};
		nu = varargin{4};
	case 5	% MOGI(R,F,A,?�,P)
		a = varargin{3};
		mu = varargin{4};
		p = varargin{5};
		nu = 0.25;
	case 6	% MOGI(R,F,A,P,E,nu)
		a = varargin{3};
		p = varargin{4};		
		nu = varargin{6};
		mu = varargin{5}./(2*(1+nu));
end

if any(nargin==[3,4])
	y = v./pi;
else
	if max(max(a))/min(min(f)) > .1
		warning('Mogi: inaccurate results if F is not much greater than A.')
	end
	y = (a.^3).*p./mu;
end

R = sqrt(f.^2 + r.^2);	% radial distance from source
C = (1-nu).*y;		% intermediate constant

et = C./R.^3;		% tangential horizontal strain
ur = r.*et;		% radial horizontal displacement
uz = f.*et;		% vertical displacement

if nargout > 2
	dt = 3*C.*f.*r./R.^5;	% tilt
end
if nargout > 3
	er = C.*(f.^2 - 2*r.^2)./R.^5;	% radial horizontal strain
end
end
%% ====================================================
function [X,Y]=PLTXY(ALAT,ALON,ALAT0,ALON0)
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN 
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
RLAT = RD.*ALAT;
SLAT = sin(RLAT);
CLAT = cos(RLAT);
V2   = 1.0 + E12.*CLAT.^2;
AL   = ALON-ALON0;
PH1  = ALAT + (V2.*AL.^2.*SLAT.*CLAT)./(2.0*D);
RPH1 = PH1.*RD;
RPH2 = (PH1 + ALAT0).*0.5.*RD;
R    = A.*(1.0-E2)./sqrt((1.0-E2.*sin(RPH2).^2).^3);
AN   = A./sqrt(1.0-E2.*sin(RPH1).^2);
C1   = D./R;
C2   = D./AN;
Y    = (PH1-ALAT0)./C1;
X    = (AL.*CLAT)./C2+(AL.^3.*CLAT.*cos(2.0.*RLAT))./(6.0.*C2.*D.^2);
end
%% CONVERT TO XYZ FROM ELL AT SURFACE
function [OOxyz]=conv2ell(Olat,Olon)
Olat=Olat(:);
Olon=Olon(:);
deg2rad=pi/180;
[Oxyz(:,1),Oxyz(:,2),Oxyz(:,3)]=ell2xyz(Olat,Olon,0);
Oxyz = Oxyz*1e3;
OOxyz=[Oxyz sin(Olat*deg2rad) sin(Olon*deg2rad) cos(Olat*deg2rad) cos(Olon*deg2rad)];
end

%% CALC MOTION BLOCKS
function Est_Motion_BLOCKS(DIR,TCHA,BLK)
% 
exid=exist([DIR,'/rigid']);
if exid~=7; mkdir([DIR,'/rigid']); end
FIDbound   =fopen([DIR,'/rigid/boundary_vector.txt'],'w');
FIDboundval=fopen([DIR,'/rigid/boundary_vector_value.txt'],'w');
FIDlateral   =fopen([DIR,'/rigid/boundary_vector_lateral.txt'],'w');
FIDlateralval=fopen([DIR,'/rigid/boundary_vector_lateral_value.txt'],'w');
FIDdip   =fopen([DIR,'/rigid/boundary_vector_dip.txt'],'w');
FIDdipval=fopen([DIR,'/rigid/boundary_vector_dip_value.txt'],'w');
for NB1=1:BLK(1).NBlock
  BLK(NB1).POL=[TCHA.AVEPOL(3.*NB1-2,:);TCHA.AVEPOL(3.*NB1-1,:);TCHA.AVEPOL(3.*NB1,:)];
  for NB2=NB1+1:BLK(1).NBlock
    BLK(NB2).POL(:)=[TCHA.AVEPOL(3.*NB2-2,:);TCHA.AVEPOL(3.*NB2-1,:);TCHA.AVEPOL(3.*NB2,:)];
    if ~isempty(BLK(1).BOUND(NB1,NB2).LAT) 
      calc_relvelo(BLK,NB1,NB2,FIDbound,FIDboundval,FIDlateral,FIDlateralval,FIDdip,FIDdipval)
    end
  end
end
fclose(FIDbound);
fclose(FIDboundval);
fclose(FIDlateral);
fclose(FIDlateralval);
fclose(FIDdip);
fclose(FIDdipval);
% 
end
%% 
function calc_relvelo(BLK,NB1,NB2,FIDbound,FIDboundval,FIDlateral,FIDlateralval,FIDdip,FIDdipval)
% Calculate relative motion with lateral and normal direction at block boundaries
cBOUNDLON=mean(BLK(1).BOUND(NB1,NB2).LON);
cBOUNDLAT=mean(BLK(1).BOUND(NB1,NB2).LAT);
[BLK1.XY(:,1),BLK1.XY(:,2)]=PLTXY(BLK(NB1).LAT,BLK(NB1).LON,cBOUNDLAT,cBOUNDLON);
[BOUND.XY(:,1),BOUND.XY(:,2)]=PLTXY(BLK(1).BOUND(NB1,NB2).LAT,BLK(1).BOUND(NB1,NB2).LON,cBOUNDLAT,cBOUNDLON);
BOUND.cXY(:,1)=(BOUND.XY(1:end-1,1)+BOUND.XY(2:end,1))./2;
BOUND.cXY(:,2)=(BOUND.XY(1:end-1,2)+BOUND.XY(2:end,2))./2;
% [BOUND.cLAT,BOUND.cLON]=XYTPL(BOUND.cXY(:,1),BOUND.cXY(:,2),cBOUNDLAT,cBOUNDLON);
BOUND.cLAT=(BLK(1).BOUND(NB1,NB2).LAT(1:end-1)+BLK(1).BOUND(NB1,NB2).LAT(2:end))./2;
BOUND.cLON=(BLK(1).BOUND(NB1,NB2).LON(1:end-1)+BLK(1).BOUND(NB1,NB2).LON(2:end))./2;
BOUND.cxyz=conv2ell(BOUND.cLAT,BOUND.cLON);

UV=[0 0 1];
BOUND.vecXY=BOUND.XY(2:end,:)-BOUND.XY(1:end-1,:);
BOUND.norXY=[UV(2).*               0-UV(3).*BOUND.vecXY(:,2)...
             UV(3).*BOUND.vecXY(:,1)-UV(1).*0               ...
             UV(1).*BOUND.vecXY(:,2)-UV(2).*BOUND.vecXY(:,1)];

BOUND.VEL=pole2velo((BLK(NB2).POL(:)-BLK(NB1).POL(:))',BOUND.cxyz);
BOUND.VELx=BOUND.VEL(1:2:end);
BOUND.VELy=BOUND.VEL(2:2:end);
BOUND.VELstr=(BOUND.VELx.*BOUND.vecXY(:,1)+BOUND.VELy.*BOUND.vecXY(:,2)).*BOUND.vecXY...
           ./(BOUND.vecXY(:,1).^2+BOUND.vecXY(:,2).^2);
BOUND.VELdip=(BOUND.VELx.*BOUND.norXY(:,1)+BOUND.VELy.*BOUND.norXY(:,2)).*BOUND.norXY(:,1:2)...
           ./(BOUND.norXY(:,1).^2+BOUND.norXY(:,2).^2);
sclSTR=zeros(size(BOUND.VELx));
sclDIP=zeros(size(BOUND.VELx));
sclVEL=sqrt(BOUND.VELx.^2+BOUND.VELy.^2);
for ii=1:length(BLK(1).BOUND(NB1,NB2).LAT)-1
  ST=[BOUND.VELstr(ii,:) 0]; sclSTR(ii)=sqrt(ST(1).^2+ST(2).^2);
  DP=[BOUND.VELdip(ii,:) 0]; sclDIP(ii)=sqrt(DP(1).^2+DP(2).^2);
  NV=cross(UV,ST)./norm(cross(UV,ST),2);
  DP=DP./norm(DP,2);
  DISTlateral=BOUND.cXY(ii,:)+NV(1:2);
  DISTdipping=BOUND.cXY(ii,:)+DP(1:2);
  inIDst=inpolygon(DISTlateral(1),DISTlateral(2),BLK1.XY(:,1),BLK1.XY(:,2));
  inIDdp=inpolygon(DISTdipping(1),DISTdipping(2),BLK1.XY(:,1),BLK1.XY(:,2));
  if inIDst==1; sclSTR(ii)=-1*sclSTR(ii); end % Left lateral
  if inIDdp~=1; sclDIP(ii)=-1*sclDIP(ii); end % Open
  fprintf(FIDbound,'> -Z %f\n',sclVEL(ii));
  fprintf(FIDbound,'%f %f\n',BLK(1).BOUND(NB1,NB2).LON(ii),  BLK(1).BOUND(NB1,NB2).LAT(ii)  );
  fprintf(FIDbound,'%f %f\n',BLK(1).BOUND(NB1,NB2).LON(ii+1),BLK(1).BOUND(NB1,NB2).LAT(ii+1));
  fprintf(FIDboundval,'%f %f %10.2f\n',BOUND.cLON(ii),BOUND.cLAT(ii),sclVEL(ii));
  fprintf(FIDlateral,'> -Z %f\n',sclSTR(ii));
  fprintf(FIDlateral,'%f %f\n',BLK(1).BOUND(NB1,NB2).LON(ii),  BLK(1).BOUND(NB1,NB2).LAT(ii)  );
  fprintf(FIDlateral,'%f %f\n',BLK(1).BOUND(NB1,NB2).LON(ii+1),BLK(1).BOUND(NB1,NB2).LAT(ii+1));
  fprintf(FIDlateralval,'%f %f %10.2f\n',BOUND.cLON(ii),BOUND.cLAT(ii),sclSTR(ii));
  fprintf(FIDdip,'> -Z %f\n',sclDIP(ii));
  fprintf(FIDdip,'%f %f\n',BLK(1).BOUND(NB1,NB2).LON(ii),  BLK(1).BOUND(NB1,NB2).LAT(ii)  );
  fprintf(FIDdip,'%f %f\n',BLK(1).BOUND(NB1,NB2).LON(ii+1),BLK(1).BOUND(NB1,NB2).LAT(ii+1));
  fprintf(FIDdipval,'%f %f %10.2f\n',BOUND.cLON(ii),BOUND.cLAT(ii),sclDIP(ii));
end
end
%% PLATE MOTION DUE TO EULER POLE (XYZ)
function Vneu=pole2velo(Pxyz,Oxyz)
% pole2velo Convert velocity from Euler pole. Vectorized.
[Nobs,~]=size(Oxyz);
[Npol,~]=size(Pxyz);
Vxyz=zeros(Npol,3,'single'); 
Vneu=zeros(2.*Nobs,Npol,'single');
%
for N=1:Nobs
  Vxyz(:,1) = -Oxyz(N,2).*Pxyz(:,3) + Pxyz(:,2).*Oxyz(N,3);
  Vxyz(:,2) = -Oxyz(N,3).*Pxyz(:,1) + Pxyz(:,3).*Oxyz(N,1);
  Vxyz(:,3) = -Oxyz(N,1).*Pxyz(:,2) + Pxyz(:,1).*Oxyz(N,2);
  Vneu(2.*N-1,:) =            -Oxyz(N,5).*Vxyz(:,1) ...
                              +Oxyz(N,7).*Vxyz(:,2); %E
  Vneu(2.*N,:)   = -Oxyz(N,4).*Oxyz(N,7).*Vxyz(:,1) ...
                   -Oxyz(N,4).*Oxyz(N,5).*Vxyz(:,2) ...
                              +Oxyz(N,6).*Vxyz(:,3); %N
end
end
%% MAKE FIGURES
function [BLK,OBS]=MAKE_FIGS(BLK,OBS)
figure(200); clf(200)
for N=1:BLK(1).NBlock
  if OBS(N).NBLK~=0
    hold on
    quiver(zeros(1,OBS(N).NBLK),zeros(1,OBS(N).NBLK),OBS(N).EVE,OBS(N).NVE,0);
  end
end
%
figure(210); clf(210)
PLON=[];PLAT=[];EVEL=[];NVEL=[];
for N=1:BLK(1).NBlock
  plot(BLK(N).LON,BLK(N).LAT);
  hold on
  PLON=[PLON; OBS(N).LON'];
  PLAT=[PLAT; OBS(N).LAT'];
  EVEL=[EVEL; OBS(N).EEV];
  NVEL=[NVEL; OBS(N).ENV];
end
% text(OBS(1).ALON,OBS(1).ALAT,OBS(1).NAME) 
% hold on
quiver(PLON,PLAT,EVEL,NVEL,'blue');
hold on
quiver(OBS(1).ALON,OBS(1).ALAT,OBS(1).EVEC,OBS(1).NVEC,'green');
%
figure(220); clf(220)
LAT=[];LON=[];VEL=[];
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    if ~isempty(BLK(1).BOUND(NB1,NB2).VEL)
      VEL=[VEL;sqrt(BLK(1).BOUND(NB1,NB2).VEL(1:2:end).^2+...
                    BLK(1).BOUND(NB1,NB2).VEL(2:2:end).^2)];
      LON=[LON; BLK(1).BOUND(NB1,NB2).LON];
      LAT=[LAT; BLK(1).BOUND(NB1,NB2).LAT];
    end
  end
end
for N=1:BLK(1).NBlock
  plot(BLK(N).LON,BLK(N).LAT,'-k');
  hold on
end
scatter(LON,LAT,20,VEL);
hold on
for N=1:BLK(1).NBlock
  text(mean(BLK(N).LON),mean(BLK(N).LAT),num2str(N));
  hold on
end
drawnow
end
%====================================================
function [LAT,LON]=XYTPL(X,Y,ALAT0,ALON0)
%-------------------------------------------------------------------------
%  PLTXY TRANSFORMS (X,Y) TO (ALAT,ALONG)
%  TRANSFORMATION  BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------------------------------------------------------------
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
RLATO = ALAT0.*RD;
SLATO = sin(RLATO);
R     = A.*(1-E2)./sqrt((1-E2.*SLATO.^2).^3);
AN    = A./sqrt(1.0-E2.*SLATO.^2);
V2    = 1 + E12.*cos(RLATO).^2;
C1    = D./R;
C2    = D./AN;
PH1   = ALAT0+C1.*Y;
RPH1  = PH1.*RD;
TPHI1 = tan(RPH1);
CPHI1 = cos(RPH1);
LAT   = PH1-(C2.*X).^2.*V2.*TPHI1./(2.*D);
LON   = ALON0+C2.*X./CPHI1-(C2.*X).^3.*(1.0+2.*TPHI1.^2)./(6.*D.^2.*CPHI1);
end
%% MAKE PARTS OF GREEN FUNCTION AT MESH GRID
function [Dg,Gg,GRD]=ReshapeGreen_PAIR(BLK,GRD,TRIg,PAIR)
% 
NGRD=size(GRD(1).ALON,2);
% 
TMP.P=zeros(3*NGRD,3.*BLK(1).NBlock);
MC=1;
MT=1;
MR=1;
for NB1=1:BLK(1).NBlock
  Dg(1).grdid=zeros(1,NGRD);
  Dg(1).grdid(1,GRD(1).ABLK==NB1)=true;
  Dg(1).GRDID(:,NB1)=reshape(repmat(Dg(1).grdid,3,1),3*NGRD,1);
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(TRIg(1).BOUND(NB1,NB2).clon,2);
    if NF~=0
      FLAG=0;
      for PN=1:size(PAIR,1)
        PAIRID=ismember([NB1 NB2],PAIR(PN,:));
        ISPAIR=sum(PAIRID);
        if ISPAIR==2;FLAG=1;break;end
      end
      if FLAG~=1
        TMP.C(1:3*NGRD,MC     :MC+  NF-1)=zeros(size(TRIg(1).BOUND(NB1,NB2).GSTR));
        TMP.C(1:3*NGRD,MC+  NF:MC+2*NF-1)=zeros(size(TRIg(1).BOUND(NB1,NB2).GDIP));
        TMP.C(1:3*NGRD,MC+2*NF:MC+3*NF-1)=zeros(size(TRIg(1).BOUND(NB1,NB2).GTNS));
      else
        TMP.C(1:3*NGRD,MC     :MC+  NF-1)=TRIg(1).BOUND(NB1,NB2).GSTR;
        TMP.C(1:3*NGRD,MC+  NF:MC+2*NF-1)=TRIg(1).BOUND(NB1,NB2).GDIP;
        TMP.C(1:3*NGRD,MC+2*NF:MC+3*NF-1)=TRIg(1).BOUND(NB1,NB2).GTNS;
      end
      MC=MC+3*NF;
      MT=MT+2*NF;
      MR=MR+  NF;
    end
  end
%   
  IND=GRD(1).ABLK==NB1;  
  NIND=[zeros(size(IND)),IND,zeros(size(IND))]; NIND=logical(reshape(NIND',3*NGRD,1));
  EIND=[IND,zeros(size(IND)),zeros(size(IND))]; EIND=logical(reshape(EIND',3*NGRD,1));
  TMP.P(EIND,3*NB1-2)=-GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,3);
  TMP.P(EIND,3*NB1-1)=-GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,3);
  TMP.P(EIND,3*NB1  )= GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,2)                    +GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,1);
  TMP.P(NIND,3*NB1-2)= GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,3)+GRD(1).AXYZ(IND,6).*GRD(1).AXYZ(IND,2);
  TMP.P(NIND,3*NB1-1)=-GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,3)-GRD(1).AXYZ(IND,6).*GRD(1).AXYZ(IND,1);
  TMP.P(NIND,3*NB1  )= GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,7).*GRD(1).AXYZ(IND,2)-GRD(1).AXYZ(IND,4).*GRD(1).AXYZ(IND,5).*GRD(1).AXYZ(IND,1);
end
% 
Gg(1).C  =TMP.C;
Gg(1).P  =TMP.P;
end
%% MAKE MATRIX
function [GRDvec,GRD]=calc_vector_atmesh_pair(BLK,TCHA,D,G,GRD,TRIg,PAIR)
[~,Gg,GRD]=ReshapeGreen_PAIR(BLK,GRD,TRIg,PAIR);
% 
Mp.SMP=TCHA.AVEPOL;
Mc.SMPMAT=repmat(TCHA.AVEFLT,3,D.CNT);
Mc.SMPMAT=Mc.SMPMAT(D.MID);
GRDvec.RIG=Gg.P*Mp.SMP;
GRDvec.ELA=Gg.C*((G.TB*Mp.SMP).*D(1).CFINV.*Mc.SMPMAT);
GRDvec.SUM=GRDvec.RIG+GRDvec.ELA;
% 
end