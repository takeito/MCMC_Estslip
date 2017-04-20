function EST_BLOCK
% Estimate BLOCK MOTION
% Code by T.ITO 2016/03/02
%
warning('off','all')
%FileOBS='./colombia_data_set/GNSS_ITRF2008_Colombia_matlab.txt'; %DATOS 201603
%FileOBS='./colombia_data_set/GNSS_segal.ABSO.txt'; %DATOS 201608 63 estaciones
FileOBS='./colombia_data_set/GNSS_segal.ABSO53.txt'; %DATOS 201608 53 estaciones
%Archivos de bloques
DIRBlock='./colombia_data_set/BLOCK_PB2003/';
%DIRBlock='./colombia_data_set/BLOCK_8/';
%DIRBlock='./colombia_data_set/BLOCK_9/';
%DIRBlock='./colombia_data_set/BLOCK_12/';
%DIRBlock='./colombia_data_set/BLOCK_12_rev2/';
%DIRBlock='./colombia_data_set/BLOCK_12_rev3/';
%DIRBlock='./colombia_data_set/BLOCK_14/';
%DIRBlock='./colombia_data_set/BLOCK_15/';
%DIRBlock='./colombia_data_set/BLOCK_16/';
%DIRBlock='./colombia_data_set/BLOCK_rev2/';
%DIRBlock='./colombia_data_set/BLOCK/';
% READ OBSERVATION FILE
OBS=READ_OBS(FileOBS);
% READ BLOCK BOUNDARY FILE in DIRECTORY 
[BLK,OBS]=READ_BLOCK_BOUND(DIRBlock,OBS);
% CALC. ABIC AND BLOCK MOTION
[BLK,OBS]=CALC_AIC(BLK,OBS);
% BLOCK MOTION BETWEEN TWO BLOCKS
[BLK,OBS]=Est_Motion_BLOCKS(BLK,OBS);
% MAKE FIGURES
MAKE_FIGS(BLK,OBS);
%
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
    title(['\bf\fontsize{16}Feather Plot'])
   %Peter Bird Model Legend
   %labels={'Caribbean Plate','North Andean Block', 'Nazca Plate', 'Panamá Block', 'South America Plate'};
   %legend(labels,'Location', 'southeastoutside')
   %(8) 9 Block Model Legend
   %labels={'Caribbean Plate','Chocó Block', 'Gran Maracaibo Block','Guajira Block', 'North Andean Block','Nazca Plate', 'Panamá Block', 'South America Plate'};
   %legend(labels,'Location', 'southeastoutside')
   %9 Block Model Legend
   %labels={'Caribbean Plate','Chocó Block', 'Gran Maracaibo Block','Guajira Block', 'North Andean Block','Nazca Plate', 'Panamá Block', 'South America Plate'};
   %legend(labels,'Location', 'southeastoutside')
   %12 Block Model Legend
   %labels={'Baudó Block','Caribbean Plate','Chocó Block', 'Gran Maracaibo Block', 'Guajira Block', 'North Andean Block','Nazca Plate', 'Panamá Block', 'Romeral Block', 'Sinú-San Jacinto Block', 'South America Plate'};
   %legend(labels,'Location', 'southeastoutside')
   %14 Block Model Legend
   %labels={'Guajira-Gran Maracaibo Block', 'Sinú-San Jacinto Block', 'East Flank Eastern Cordillera Block', 'Chocó-Cauca-Nariño Block','Septentrional Central Cordillera Block', 'Serranía San Lucas-Valle Medio del Magdalena Block', 'West Flank Eastern Cordillera Block','Romeral-Cauca-Patia Block','Valle Superior del Magdalena Block','Caribbean Plate', 'Nazca Plate', 'Panamá Block', 'South America Plate'};
   %legend(labels,'Location', 'southeastoutside')
   %16 Block Model Legend
   %labels={'Guajira-Gran Maracaibo Block', 'Sinú-San Jacinto Block','Valle Inferior del Magdalena Block', 'East Flank Eastern Cordillera Block', 'Chocó-Cauca-Nariño Block', 'Septentrional Central Cordillera Block', 'Serranía San Lucas-Valle Medio del Magdalena Block', 'West Flank Eastern Cordillera Block','Romeral Block', 'Cauca-Patia Block','Valle Superior del Magdalena Block','Caribbean Plate', 'Nazca Plate', 'Panamá Block', 'South America Plate'};
   %legend(labels,'Location', 'southeastoutside')
   %25 Block Model Legend
   %labels={'Alta Guajira Block', 'Baja Guajira Block', 'Sierra Nevada de Santa Marta-Cesar-Rancheria Block', 'Serranía Perija Block', 'Maracaibo Block', 'Sinú-San Jacinto Block', 'Valle Inferior del Magdalena Block', 'Macizo de Santander Block', 'Pamplona Block', 'Chocó Block', 'Septentrional Central Cordillera Block', 'Serranía de San Lucas Block', 'Valle Medio del Magdalena Block', 'West Flank Eastern Cordillera Block', 'East Flank Eastern Cordillera Block', 'Cauca-Nariño Block','Romeral Block', 'Cauca-Patia Block', 'Valle Superior del Magdalena Block', 'Macizo de Garzón Block','Caribbean Plate', 'Nazca Plate', 'South America Plate','Panamá Block'};
   %legend(labels,'Location', 'southeastoutside')
   %12 Block #1 Model Legend
   %labels={'Guajira-Gran Maracaibo Block', 'Sinú-San Jacinto Block','Chocó-Cauca-Nariño Block', 'Septentrional Central Cordillera Block','Romeral-Cauca-Patia Block', 'Serranía San Lucas-Valle Medio del Magdalena Block','Valle Superior del Magdalena-East Flank Eastern Cordillera Block','Caribbean Plate', 'Nazca Plate', 'Panamá Block',  'South America Plate'};
   %legend(labels,'Location', 'southeastoutside')
   %12 Block #2 Model Legend
   labels={'Guajira-Gran Maracaibo Block', 'Sinú-San Jacinto Block','Chocó-Cauca-Nariño Block', 'Septentrional Central Cordillera Block','Romeral-Cauca-Patia Block', 'Serranía San Lucas-Valle Medio del Magdalena Block','Valle Superior del Magdalena-East Flank Eastern Cordillera Block','Caribbean Plate', 'Nazca Plate', 'Panamá Block',  'South America Plate'};
   legend(labels,'Location', 'southeastoutside')
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
a=quiver(PLON,PLAT,EVEL,NVEL);
hold on
b=quiver(OBS(1).ALON,OBS(1).ALAT,OBS(1).EVEC,OBS(1).NVEC);
title(['\bf\fontsize{16}Estimated velocity vectors VS Observed velocity vectors'])
labels1={'Estimated Velocities', 'Observed Velocities'};
legend([a,b], labels1,'Location','southeastoutside')
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
colorbar
hold on
title(['\bf\fontsize{16}Slip rate in block boundaries'])
for N=1:BLK(1).NBlock
  text(mean(BLK(N).LON),mean(BLK(N).LAT),num2str(N));
  hold on
end
end
% CALCLATION AIC AND BLOCK MOTION
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
  %tmp
  BLK(NB).name=file(NB).name;
  BLK(NB).LON=tmp(:,1);
  BLK(NB).LAT=tmp(:,2);
  %BLK(NB)
end
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    [~,ialon,~]=intersect(BLK(NB1).LON,BLK(NB2).LON);
    [~,ialat,~]=intersect(BLK(NB1).LAT,BLK(NB2).LAT);
    [Ca,~,~]=intersect(ialat,ialon);
    BLK(1).BOUND(NB1,NB2).LAT=BLK(NB1).LAT(Ca);
    %BLK(1).BOUND(NB1,NB2).LAT
    BLK(1).BOUND(NB1,NB2).LON=BLK(NB1).LON(Ca);
    %BLK(1).BOUND(NB1,NB2).LON
    BLK(1).BOUND(NB1,NB2).BXYZ=conv2ell(BLK(1).BOUND(NB1,NB2).LAT,BLK(1).BOUND(NB1,NB2).LON);
    %BLK(1).BOUND(NB1,NB2).BXYZ
  end
end
for N=1:BLK(1).NBlock
  IND=inpolygon(OBS(1).ALON,OBS(1).ALAT,BLK(N).LON,BLK(N).LAT);
  %IND
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
  %OBS(N).OXYZ
  OBS(N).Vne=reshape([OBS(1).EVEC(IND); OBS(1).NVEC(IND)],OBS(N).NBLK.*2,1);
  %OBS(N).Vne
  OBS(N).Vww=reshape([OBS(1).EERR(IND); OBS(1).NERR(IND)],OBS(N).NBLK.*2,1);
end
end
%% READ OBSERVATION DATA
function OBS=READ_OBS(FileOBS)
%-------------------
% INPUT format Observations:
% site_name lon lat h EW_comp. NS_comp. UD_comp. ERR_EW ERR_NS ERR_UD 
%-------------------
Fid_OBS=fopen(FileOBS,'r');
N=0;
%FileOBS
while 1
  tline=fgetl(Fid_OBS);
  %tline
  if ~ischar(tline); break; end
  str=strsplit(tline);
  %str
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
%R
%R(:,1) = -Oxyz(:,2).*pvec(3) + pvec(2).*Oxyz(:,3);
%R(:,2) = -Oxyz(:,3).*pvec(1) + pvec(3).*Oxyz(:,1);
%R(:,3) = -Oxyz(:,1).*pvec(2) + pvec(1).*Oxyz(:,2);
for N=1:Nobs
  R(2.*N-1,1)=-Oxyz(N,7).*Oxyz(N,3);
  %Oxyz(N,7)
  %R(2.*N-1,1)
  R(2.*N-1,2)=-Oxyz(N,5).*Oxyz(N,3);
  R(2.*N-1,3)= Oxyz(N,5).*Oxyz(N,2)+Oxyz(N,7).*Oxyz(N,1);
  R(2.*N,1)  = Oxyz(N,4).*Oxyz(N,5).*Oxyz(N,3)+Oxyz(N,6).*Oxyz(N,2);
  R(2.*N,2)  =-Oxyz(N,4).*Oxyz(N,7).*Oxyz(N,3)-Oxyz(N,6).*Oxyz(N,1);
  R(2.*N,3)  = Oxyz(N,4).*Oxyz(N,7).*Oxyz(N,2)-Oxyz(N,4).*Oxyz(N,5).*Oxyz(N,1);
  %R
end
[PL,~,Sigma]=lscov(R,Vne,w);
%Donde PL corresponde con la tasa de rotación del polo en rad/yr en
%coordenadas cartesianas geocéntricas
%PL
%R
%Vne
%w
%PL
EVne=R*PL;
%EVne
end
%% PLATE MOTION DUE TO EULER POLE (XYZ)
function Vneu=pole2velo(Pxyz,Oxyz)
% pole2velo Convert velocity from Euler pole. Vectorized.
[Nobs,~]=size(Oxyz);
[Npol,~]=size(Pxyz);
%size(Pxyz)
%size(Oxyz)
%Oxyz
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
%Vxyz
%Vneu
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