function EST_BLOCK_slip
% Estimate BLOCK MOTION
% Code by T.ITO 2016/03/02
%
warning('off','all')
FileOBS='./GNSS_ITRF2008_Colombia_matlab.txt';
%DIRBlock='./colombia_data_set/BLOCK_PB2003/';
%DIRBlock='./colombia_data_set/BLOCK_8/';
%DIRBlock='./colombia_data_set/BLOCK_9/';
%DIRBlock='./colombia_data_set/BLOCK_12/';
%DIRBlock='./colombia_data_set/BLOCK_12_rev2/';
DIRBlock='./BLOCK/';
%DIRBlock='./colombia_data_set/BLOCK_14/';
%DIRBlock='./colombia_data_set/BLOCK_15/';
%DIRBlock='./colombia_data_set/BLOCK_16/';
%DIRBlock='./colombia_data_set/BLOCK_rev2/';
DIRBlock_Interface='./BLOCK_Int/';
%DIRBlock='./colombia_data_set/BLOCK_14/';
%DIRBlock='./colombia_data_set/BLOCK_15/';
%DIRBlock='./colombia_data_set/BLOCK_16/';
%DIRBlock='./colombia_data_set/BLOCK_rev2/';
% READ OBSERVATION FILE
OBS=READ_OBS(FileOBS);
% READ BLOCK BOUNDARY FILE in DIRECTORY
[BLK,OBS]=READ_BLOCK_BOUND(DIRBlock,OBS);
% READ BLOCK INTERFACE BOUNDARY in DIRECTORY 
[BLK]=READ_BLOCK_INTERFACE(BLK,DIRBlock_Interface);
% CALC. ABIC AND BLOCK MOTION
[BLK,OBS]=CALC_AIC(BLK,OBS);
% BLOCK MOTION BETWEEN TWO BLOCKS
[BLK,OBS]=Est_Motion_BLOCKS(BLK,OBS);
% MAKE FIGURES
MAKE_FIGS(BLK,OBS);
%
end
%% READ PLATE INTERFACE
function [BLK]=READ_BLOCK_INTERFACE(BLK,DIRBLK)
%
int_lat=0.25;
int_lon=0.25;
dep_limit=150;
dep_limit_low=10;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    sub_f=fullfile(DIRBLK,['B_',num2str(NB1),'_',num2str(NB2),'.txt']);
    Fid=fopen(sub_f,'r');
    if Fid >= 0
      dep_blk=textscan(Fid,'%f%f%f'); fclose(Fid);
      dep_blk=cell2mat(dep_blk);
      F=scatteredInterpolant(dep_blk(:,1),dep_blk(:,2),dep_blk(:,3),'natural');
      BO_f=fullfile(DIRBLK,['BO_',num2str(NB1),'_',num2str(NB2),'.txt']);
      Fid=fopen(BO_f,'r');
      if Fid >= 0
        bound_blk=textscan(Fid,'%f%f'); fclose(Fid);     
        bound_blk=cell2mat(bound_blk);
        min_lon=min(bound_blk(:,1)); max_lon=max(bound_blk(:,1));
        min_lat=min(bound_blk(:,2)); max_lat=max(bound_blk(:,2));
        [slon,slat]=meshgrid(min_lon:int_lon:max_lon,min_lat:int_lat:max_lat);
        ID=inpolygon(slon,slat,bound_blk(:,1),bound_blk(:,2));
      else
        IDB=boundary(dep_blk(:,1),dep_blk(:,2));
        min_lon=min(dep_blk(:,1)); max_lon=max(dep_blk(:,1));
        min_lat=min(dep_blk(:,2)); max_lat=max(dep_blk(:,2));              
        [slon,slat]=meshgrid(min_lon:int_lon:max_lon,min_lat:int_lat:max_lat);
        ID=inpolygon(slon,slat,dep_blk(IDB,1),dep_blk(IDB,2));
      end
      slon=slon(ID);
      slat=slat(ID);
      sdep=F(slon,slat);
      ID=find(sdep<dep_limit);
      slon=slon(ID);
      slat=slat(ID);
      sdep=F(slon,slat);
      BLK(1).BOUND(NB1,NB2).slon=slon;
      BLK(1).BOUND(NB1,NB2).slat=slat;
      BLK(1).BOUND(NB1,NB2).sdep=sdep;
      BLK(1).BOUND(NB1,NB2).tri =delaunay(slon,slat);
    else
      LENG=length(BLK(1).BOUND(NB1,NB2).LON);
      BLK(1).BOUND(NB1,NB2).slon=[BLK(1).BOUND(NB1,NB2).LON BLK(1).BOUND(NB1,NB2).LON];
      BLK(1).BOUND(NB1,NB2).slat=[BLK(1).BOUND(NB1,NB2).LAT BLK(1).BOUND(NB1,NB2).LAT];
      BLK(1).BOUND(NB1,NB2).sdep=[zeros(size(BLK(1).BOUND(NB1,NB2).LAT)) dep_limit_low.*ones(BLK(1).BOUND(NB1,NB2).LAT)];
      BLK(1).BOUND(NB1,NB2).tri(   1:LENG-1    ,1:3)=[     1:  (LENG-1); 2:LENG; LENG+1:2*LENG-1];
      BLK(1).BOUND(NB1,NB2).tri(LENG:2*(LENG-1),1:3)=[LENG+1:2*(LENG-1); 2:LENG; LENG+2:2*LENG-1];      
    end
  end
end
%====================================================
dep_sub=textscan(Fid,'%f%f%f');
fclose(Fid);
dep_sub=cell2mat(dep_sub);
%====================================================
Fid=fopen(bound_f);
bound=textscan(Fid,'%f%f%f');
fclose(Fid);
bound=cell2mat(bound);
%====================================================
Fid=fopen(inip_f);
ini_point=textscan(Fid,'%f%f%f');
fclose(Fid);
ini_point=cell2mat(ini_point);
[int_size,~]=size(ini_point);
%====================================================
F=scatteredInterpolant(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');
min_lon=min(bound(:,1)); max_lon=max(bound(:,1));
min_lat=min(bound(:,2)); max_lat=max(bound(:,2));
figure(10); clf
plot(bound(:,1),bound(:,2),'r')
hold on
%====================================================
% initial Fixed points -------
s.lon=ini_point(:,1);
s.lat=ini_point(:,2);
s.dep=ini_point(:,3);
%====================================================
n=int_size;
while n<int_mesh
  slat=(max_lat-min_lat).*rand(1)+min_lat;
  slon=(max_lon-min_lon).*rand(1)+min_lon;
  ID=inpolygon(slon,slat,bound(:,1),bound(:,2));
  if ID==1
    n=n+1;
    s.lat(n)=slat;
    s.lon(n)=slon;
    s.dep(n)=F(slon,slat);
    if rem(n,round(int_mesh/10))==1;
      plot3(s.lon,s.lat,s.dep,'.')
      pause(.1)
    end
  end
end
plot3(s.lon,s.lat,s.dep,'.')
%====================================================
tri = delaunay(s.lon,s.lat);
%====================================================
glon=mean(s.lon(tri),2);
glat=mean(s.lat(tri),2);
ID=inpolygon(glon,glat,bound(:,1),bound(:,2));
s.tri=tri(ID==1,:);
%
figure(20); clf
plot(bound(:,1),bound(:,2),'r')
hold on
triplot(s.tri,s.lon,s.lat)
pause(.1)
s.bound=bound;
s.dep_sub=dep_sub;
end
%% MAKE GREEN FUNCTION
function [TRI]=READ_FAULTS_GREEN_TRI(BLK,OBS)
% Coded by Takeo Ito 2015/11/11 (ver 1.0)
PR=0.25;
NF=0;
% TODO: READ MULTIPLE BLOCK VELOCITY
load PRM.TRI_F
%
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    BLK(1).BOUND(NB1,NB2).BLAT
    BLK(1).BOUND(NB1,NB2).BLON
    BLK(1).BOUND(NB1,NB2).BDEP
    BLK(1).BOUND(NB1,NB2).BXYZ=conv2ell(BLK(1).BOUND(NB1,NB2).LAT,BLK(1).BOUND(NB1,NB2).LON);
    
  end
end

while 1
  NF=NF+1;
  loc_f=fscanf(Fid,'%f %f %f \n', [3 3]);
  [~] = fgetl(Fid);
  TRI(1).LON(NF,:)=loc_f(1,:);%Lon
  TRI(1).LAT(NF,:)=loc_f(2,:);%Lat
  TRI(1).HIG(NF,:)=loc_f(3,:);%Hight
  tline = fgetl(Fid); if ~ischar(tline); break; end
end
fclose(Fid);
%
ND=size(OBS(1).ALAT,2);
TRI(1).GSTR=zeros(3*ND,NF);
TRI(1).GDIP=zeros(3*ND,NF);
%
fprintf('==================\nNumber of TRI sub-faults : %i \n',NF)
%
ALAT=mean(OBS(1).ALAT(:));
ALON=mean(OBS(1).ALON(:));
[OBSx,OBSy]=PLTXY(OBS(1).ALAT,OBS(1).ALON,ALAT,ALON);
OBSz=OBS(1).AHIG;
%
for N=1:NF
  [TRIx,TRIy]=PLTXY(TRI(1).LAT(N,:),TRI(1).LON(N,:),ALAT,ALON);
  TRIz=-1.*TRI(1).HIG(N,:);
  F_LOC=[TRI(1).LAT(N,:);TRI(1).LON(N,:);TRI(1).HIG(N,:)]';
  [F,DA,STR,DIP,NV,ST,DP]=EST_FAULT_TRI(F_LOC);
  TRI(1).CLAT(N)=F(1);
  TRI(1).CLON(N)=F(2);
  TRI(1).CHIG(N)=F(3);
  TRI(1).DA(N)=DA;
  TRI(1).STR(N)=STR;
  TRI(1).DIP(N)=DIP;
  TRI(1).NV(N,:)=NV;
  TRI(1).ST(N,:)=ST;
  TRI(1).DP(N,:)=DP;
  U=CalcTriDisps(OBSx,OBSy,OBSz,TRIx,TRIy,TRIz,PR,1,0,0);
  TRI(1).GSTR(1:3:3*ND,N)=U.x; %E
  TRI(1).GSTR(2:3:3*ND,N)=U.y; %N
  TRI(1).GSTR(3:3:3*ND,N)=U.z; %D
  U=CalcTriDisps(OBSx,OBSy,OBSz,TRIx,TRIy,TRIz,PR,0,1,0);
  TRI(1).GTNS(1:3:3*ND,N)=U.x; %E
  TRI(1).GTNS(2:3:3*ND,N)=U.y; %N
  TRI(1).GTNS(3:3:3*ND,N)=U.z; %D 
  U=CalcTriDisps(OBSx,OBSy,OBSz,TRIx,TRIy,TRIz,PR,0,0,1);
  TRI(1).GDIP(1:3:3*ND,N)=U.x; %E
  TRI(1).GDIP(2:3:3*ND,N)=U.y; %N
  TRI(1).GDIP(3:3:3*ND,N)=U.z; %D 
  if mod(N,ceil(NF/5)) == 1
    fprintf('MAKE GREEN at TRI sub-faults : %i / %i \n',N,NF)
  end
end
fprintf('Size of Green Matrix : %i x %i \n',...
             size(TRI(1).GSTR,1),size(TRI(1).GSTR,2))
disp('PASS READ_FAULTS_GREEN_TRI')
end
%% ESTIMATE FAULT PARAMETERS FOR TRI
function [FLOC,DA,STR,DIP,NV,ST,DP]=EST_FAULT_TRI(loc_f)
% Coded by Takeo Ito 2015/11/11 (ver 1.0)
[X,Y]=PLTXY(loc_f(:,2),loc_f(:,1),loc_f(1,2),loc_f(1,1));
[DA]=AREA_TRI(X,Y,loc_f(:,3));
FLOC=mean(loc_f);
[STR,DIP,NV,ST,DP]=EST_STRDIP_TRI(X,Y,loc_f(:,3));
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
function [STR,DIP,NV,ST,DP]=EST_STRDIP_TRI(X,Y,Z)
%==========
% CALC. STR AND DIP ON FAULT
% CODE BY T.ITO (2006/03/04)
% Modified by T.ITO (2015/11/13)
% Modified by T.ITO (2016/02/16)
% DEPTH IS MINUS
%==========
NV=cross([X(2);Y(2);Z(2)]-[X(1);Y(1);Z(1)], [X(3);Y(3);Z(3)]-[X(1);Y(1);Z(1)]);
NV=NV./norm(NV);
if (NV(3) < 0); NV = -NV; end; % Enforce clockwise circulation
ST=[cos(atan2(NV(2),NV(1))) sin(atan2(NV(2),NV(1))) 0];
DP=cross(NV,ST);
STR=(180./pi).*atan2(NV(2),NV(1));
DIP=(180./pi).*acos(NV(3));
end
%% MAKE FIGURES
function [BLK,OBS]=Est_Motion_BLOCKS(BLK,OBS)
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    if ~isempty(BLK(1).BOUND(NB1,NB2).LAT) 
      BLK(1).BOUND(NB1,NB2).VEL=pole2velo((BLK(NB2).POL(:)-BLK(NB1).POL(:))',BLK(1).BOUND(NB1,NB2).BXYZ);
    end
  end
end
end
%% MAKE FIGURES
function [BLK,OBS]=MAKE_FIGS(BLK,OBS)
figure(100); clf(100)
for N=1:BLK(1).NBlock
  if OBS(N).NBLK~=0
    hold on
    quiver(zeros(1,OBS(N).NBLK),zeros(1,OBS(N).NBLK),OBS(N).EVE,OBS(N).NVE,0);
  end
end
%
figure(200); clf(200)
PLON=[];PLAT=[];EVEL=[];NVEL=[];
for N=1:BLK(1).NBlock
  plot(BLK(N).LON,BLK(N).LAT);
  hold on
  PLON=[PLON; OBS(N).LON'];
  PLAT=[PLAT; OBS(N).LAT'];
  EVEL=[EVEL; OBS(N).EEV];
  NVEL=[NVEL; OBS(N).ENV];
end
text(OBS(1).ALON,OBS(1).ALAT,OBS(1).NAME) 
hold on
quiver(PLON,PLAT,EVEL,NVEL);
hold on
quiver(OBS(1).ALON,OBS(1).ALAT,OBS(1).EVEC,OBS(1).NVEC);
%
figure(300); clf(300)
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
end
%% CALCLATION AIC AND BLOCK MOTION
function [BLK,OBS]=CALC_AIC(BLK,OBS)
TSig=0; NumB=0;
for N=1:BLK(1).NBlock
  Sig=0;EVne=[];POLE=[0; 0; 0];
  OBS(N).EEV=OBS(N).Vne(1:2:end);
  OBS(N).ENV=OBS(N).Vne(2:2:end);
  if OBS(N).NBLK~=0
    Sig=0;
    EVne=[0 0];
    if OBS(N).NBLK>=1
      NumB=NumB+1;
      [POLE,EVne,Sig]=est_pole_w(OBS(N).OXYZ,OBS(N).Vne,OBS(N).Vww);
      TSig=TSig+Sig.*2.*OBS(N).NBLK;
    end
  end
  BLK(N).SIG=Sig;
  BLK(N).POL=POLE;
  OBS(N).EEV=EVne(1:2:end);
  OBS(N).ENV=EVne(2:2:end);
  fprintf('BLOCK=%2d NUM_OBS=%2d Sigma^2=%5.2f \n',N,OBS(N).NBLK,Sig)
  if OBS(N).NBLK>=2 
    fprintf('OBS(E,N) ')
    fprintf('%5.2f ',OBS(N).Vne);fprintf('\n')
    fprintf('EST(E,N) ')
    fprintf('%5.2f ',EVne)      ;fprintf('\n')
  fprintf('\n')
  end
end
AIC=(OBS(1).NOBS.*2).*log(TSig./(OBS(1).NOBS.*2))+2.*NumB.*3;
cAIC=AIC+2.*NumB.*3.*(NumB.*3+1)./(OBS(1).NOBS.*2-NumB.*3-1);
fprintf('Sigma^2=%8.3f AIC=%7.3f cAIC=%7.3f K=%2d\n',TSig./(OBS(1).NOBS.*2),AIC,cAIC,NumB.*3)
%
end
%% READ BLOCK BOUNDARY DATA
function [BLK,OBS]=READ_BLOCK_BOUND(DIR,OBS)
EXT='*.txt';
file=dir([DIR,EXT]);
[NBlock,~]=size(file);
BLK(1).NBlock=NBlock;
for NB=1:BLK(1).NBlock
  tmp=load(fullfile(DIR,file(NB).name));
  BLK(NB).name=file(NB).name;
  BLK(NB).LON=tmp(:,1);
  BLK(NB).LAT=tmp(:,2);
end
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    [~,ialon,~]=intersect(BLK(NB1).LON,BLK(NB2).LON);
    [~,ialat,~]=intersect(BLK(NB1).LAT,BLK(NB2).LAT);
    [Ca,~,~]=intersect(ialat,ialon);
    BLK(1).BOUND(NB1,NB2).LAT=BLK(NB1).LAT(Ca);
    BLK(1).BOUND(NB1,NB2).LON=BLK(NB1).LON(Ca);
    BLK(1).BOUND(NB1,NB2).BXYZ=conv2ell(BLK(1).BOUND(NB1,NB2).LAT,BLK(1).BOUND(NB1,NB2).LON);
  end
end
for N=1:BLK(1).NBlock
  IND=inpolygon(OBS(1).ALON,OBS(1).ALAT,BLK(N).LON,BLK(N).LAT);
  OBS(1).ABLK(IND)=N;
  OBS(N).NBLK=sum(IND);
  OBS(N).LAT=OBS(1).ALAT(IND);
  OBS(N).LON=OBS(1).ALON(IND);
  OBS(N).HIG=OBS(1).AHIG(IND);
  OBS(N).EVE=OBS(1).EVEC(IND);
  OBS(N).NVE=OBS(1).NVEC(IND);
  OBS(N).HVE=OBS(1).HVEC(IND);
  OBS(N).EER=OBS(1).EERR(IND);
  OBS(N).NER=OBS(1).NERR(IND);
  OBS(N).HER=OBS(1).HERR(IND);
  OBS(N).OXYZ=conv2ell(OBS(N).LAT,OBS(N).LON);
  OBS(N).Vne=reshape([OBS(1).EVEC(IND); OBS(1).NVEC(IND)],OBS(N).NBLK.*2,1);
  OBS(N).Vww=reshape([OBS(1).EERR(IND); OBS(1).NERR(IND)],OBS(N).NBLK.*2,1);
end
end
%% READ OBSERVATION DATA
function OBS=READ_OBS(FileOBS)
%-------------------
% INPUT format Observations:
% site_name lon lat EW_comp. NS_comp. UD_comp. ERR_EW ERR_NS ERR_UD 
%-------------------
Fid_OBS=fopen(FileOBS,'r');
N=0;
while 1
  tline=fgetl(Fid_OBS);
  if ~ischar(tline); break; end
  str=strsplit(tline);
  N=N+1;
  OBS(1).NAME(N)=cellstr(str(1));
  OBS(1).ALON(N) =str2double(cellstr(str(2))); %LON
  OBS(1).ALAT(N) =str2double(cellstr(str(3))); %LAT
  OBS(1).AHIG(N) =str2double(cellstr(str(4))); %HIG
  OBS(1).EVEC(N) =str2double(cellstr(str(5))); %E-W
  OBS(1).NVEC(N) =str2double(cellstr(str(6))); %N-S
  OBS(1).HVEC(N) =str2double(cellstr(str(7))); %U-D
  OBS(1).EERR(N) =str2double(cellstr(str(8))); %E-W
  OBS(1).NERR(N) =str2double(cellstr(str(9))); %N-S
  OBS(1).HERR(N) =str2double(cellstr(str(10))); %U-D
end
OBS(1).NOBS=N;
OBS(1).ABLK=zeros(OBS(1).NOBS,1);
end
%% Estimate BLOCK Motion
function [PL,EVne,Sigma]=est_pole_w(Oxyz,Vne,w)
[Nobs,~]=size(Oxyz);
R=zeros(Nobs.*2,3);
%R(:,1) = -Oxyz(:,2).*pvec(3) + pvec(2).*Oxyz(:,3);
%R(:,2) = -Oxyz(:,3).*pvec(1) + pvec(3).*Oxyz(:,1);
%R(:,3) = -Oxyz(:,1).*pvec(2) + pvec(1).*Oxyz(:,2);
for N=1:Nobs
  R(2.*N-1,1)=-Oxyz(N,7).*Oxyz(N,3);
  R(2.*N-1,2)=-Oxyz(N,5).*Oxyz(N,3);
  R(2.*N-1,3)= Oxyz(N,5).*Oxyz(N,2)+Oxyz(N,7).*Oxyz(N,1);
  R(2.*N,1)  = Oxyz(N,4).*Oxyz(N,5).*Oxyz(N,3)+Oxyz(N,6).*Oxyz(N,2);
  R(2.*N,2)  =-Oxyz(N,4).*Oxyz(N,7).*Oxyz(N,3)-Oxyz(N,6).*Oxyz(N,1);
  R(2.*N,3)  = Oxyz(N,4).*Oxyz(N,7).*Oxyz(N,2)-Oxyz(N,4).*Oxyz(N,5).*Oxyz(N,1);
end
[PL,~,Sigma]=lscov(R,Vne,w);
EVne=R*PL;
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
function [OOxyz]=conv2ell(Olat,Olon)
Olat=Olat(:);
Olon=Olon(:);
deg2rad=pi/180;
[Oxyz(:,1),Oxyz(:,2),Oxyz(:,3)]=ell2xyz(Olat,Olon,0);
Oxyz = Oxyz*1e3;
OOxyz=[Oxyz sin(Olat*deg2rad) sin(Olon*deg2rad) cos(Olat*deg2rad) cos(Olon*deg2rad)];
end
%% Make GREEN FUNCTION
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
   segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
   [rx,ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
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
   ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
   tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
   dsVec                    = cross(ssVec, tsVec);
   lss                      = dot(slipVec, ssVec);
   lts                      = dot(slipVec, tsVec);
   lds                      = dot(slipVec, dsVec);
%
   if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
      % First angular dislocation
      [sx1,sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
      [ux1,uy1,uz1]             = adv(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
                                   
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
%====================================================