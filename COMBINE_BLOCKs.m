%% MAIN PROGRAM
function COMBINE_BLOCKs(BDIR,BNO1,BNO2)
% ----------------------
% THIS PROGRAM COMBINE 2 BLOCKS AND SAVE THEM TO NEW FOLDER
% USAGE)
% COMBINE_BLOCKs(BDIR,BNO1,BNO2)
% BDIR:BLOCK FOLDER
% BNO1:BLOCK NUMBER THAT YOU WANT TO COMBINE WITH ANOTER ONE
% BNO2:ANOTHER ONE
% EXAMPLE)
% COMBINE_BLOCKs('BLOCK_sw_japan',10,11)
% COMBINED BLOCK FILE IS SAVED IN ./BLOCK_MODEL/MODEL_??/
% '??' IS CONFIGURED AUTOMATICALLY.
% ----------------------
% INPUT_SET='./PARAMETER/parameter.txt';
% READ PARAMETER FOR MCMC Inversion 
% [PRM]=READ_PARAMETERS(INPUT_SET);
% READ BLOCK BOUNDARY FILE in DIRECTORY
PRM.DIRBlock=BDIR;
[BLK]=READ_BLOCK_BOUND(PRM.DIRBlock);
% COMBINE BLOCKS
[BLK]=COMBINE_BOUND(BLK,BNO1,BNO2);
% SAVE COMBINED BLOCK FILE TO 'BLOCK_MODEL' FOLDER
[B,DIR]=WRITE_FILE(BLK,BNO1,BNO2);
% MAKE_FIGS(BLK);
MAKE_FIGS(BLK,B,BNO1,BNO2,DIR);
end
%% Write blocks in txt format
function [B,DIR]=WRITE_FILE(BLK,BNO1,BNO2)
DIR1='BLOCK_MODEL';
for DN=1:Inf
  DIR=['MODEL_',num2str(DN,'%02i')];
  FDIR=fullfile(DIR1,DIR);
  EXID=exist(FDIR);
  if EXID~=7; mkdir(FDIR); fprintf('%s\n',['Export folder: ',FDIR]); break; end
end
NB=1;
for FN=1:BLK(1).NBlock
  namesplit=strsplit(BLK(FN).name,{'_','.'});
  B(FN).NO=str2num(char(namesplit(1)));
  B(FN).NAME=cell2mat(namesplit(2:end-1));
  if B(FN).NO==BNO1; BNO1=FN; continue; end
  if B(FN).NO==BNO2; BNO2=FN; continue; end
  rename=[num2str(NB,'%02i'),'_',B(FN).NAME,'.txt'];
  copyfile(BLK(FN).fullname,fullfile(FDIR,rename));
  NB=NB+1;
end
B(BLK(1).NBlock+1).NO=B(BLK(1).NBlock).NO+1;
B(BLK(1).NBlock+1).NAME=[B(BNO1).NAME,'-',B(BNO2).NAME];
newfile=[num2str(NB,'%02i'),'_',B(BLK(1).NBlock+1).NAME,'.txt'];
NEWFILE=fullfile(FDIR,newfile);
LOGFILE=fullfile(FDIR,'combine.log');
FID=fopen(NEWFILE,'w');
fprintf(FID,'%15.9f %15.9f \n',[BLK(BLK(1).NBlock+1).LON BLK(BLK(1).NBlock+1).LAT]');
fclose(FID);
FIDl=fopen(LOGFILE,'w');
fprintf(FID,'%s\n',['BLOCK ',num2str(B(BNO1).NO,'%02i'),' and BLOCK ',num2str(B(BNO2).NO,'%02i'),' were combined!!']);
fprintf(FID,'%s\n',['And then BLOCK ',num2str(NB,'%02i'),' ( ',newfile, ' ) ' 'was created!!']);
fclose(FIDl);
end
%% MAKE FIGURES
function MAKE_FIGS(BLK,B,BNO1,BNO2,DIR)
% figure('Name','OBS_vector'); clf
% for N=1:BLK(1).NBlock
%   if OBS(N).NBLK~=0
%     hold on
%     quiver(zeros(1,OBS(N).NBLK),zeros(1,OBS(N).NBLK),OBS(N).EVE,OBS(N).NVE,0);
%   end
% end
% %
% figure('Name','OBS_sites, OBS_vector, Computed_vector'); clf
% PLON=[];PLAT=[];EVEL=[];NVEL=[];
% for N=1:BLK(1).NBlock
%   plot(BLK(N).LON,BLK(N).LAT);
%   hold on
%   PLON=[PLON; OBS(N).LON'];
%   PLAT=[PLAT; OBS(N).LAT'];
%   EVEL=[EVEL; OBS(N).EEV];
%   NVEL=[NVEL; OBS(N).ENV];
% end
% text(OBS(1).ALON,OBS(1).ALAT,OBS(1).NAME) 
% hold on
% quiver(PLON,PLAT,EVEL,NVEL);
% hold on
% quiver(OBS(1).ALON,OBS(1).ALAT,OBS(1).EVEC,OBS(1).NVEC);
%
figure(90); clf(90)
% LAT=[];LON=[];VEL=[];
% for NB1=1:BLK(1).NBlock
%   if NB1==BNO1||NB1==BNO2; continue; end
%   for NB2=NB1+1:BLK(1).NBlock
%     if NB1==BNO1||NB1==BNO2; continue; end
%     if ~isempty(BLK(1).BOUND(NB1,NB2).VEL)
%       VEL=[VEL;sqrt(BLK(1).BOUND(NB1,NB2).VEL(1:2:end).^2+...
%                     BLK(1).BOUND(NB1,NB2).VEL(2:2:end).^2)];
%       LON=[LON; BLK(1).BOUND(NB1,NB2).LON];
%       LAT=[LAT; BLK(1).BOUND(NB1,NB2).LAT];
%     end
%   end
% end
for N=1:BLK(1).NBlock+1
  if B(N).NO==BNO1 || B(N).NO==BNO2; continue; end
  plot(BLK(N).LON,BLK(N).LAT,'-k');
  hold on
end
% scatter(LON,LAT,20,VEL);
hold on
NB=1;
for N=1:BLK(1).NBlock+1
  if B(N).NO==BNO1 || B(N).NO==BNO2; continue; end
  text(mean(BLK(N).LON),mean(BLK(N).LAT),num2str(NB));
  NB=NB+1;
  hold on
end
title(DIR);
end
%% READ PARAMETER FILE 
function [PRM]=READ_PARAMETERS(INPUT_SET)
%
Fid=fopen(INPUT_SET,'r');
PRM.HOME_D=pwd;
FileOBS=fscanf(Fid,'%s \n',[1,1]);
PRM.FileOBS=fullfile(PRM.HOME_D,FileOBS);
[~]=fgetl(Fid);
DIRBlock=fscanf(Fid,'%s \n',[1,1]);
PRM.DIRBlock=fullfile(PRM.HOME_D,DIRBlock);
fclose(Fid);
%====================================================
fprintf('==================\nINPUT PARAMETERS\n==================\n') 
fprintf('HOME_D             : %s \n',PRM.HOME_D) 
fprintf('FileOBS            : %s \n',PRM.FileOBS) 
fprintf('DIRBlock           : %s \n',PRM.DIRBlock)
%====================================================
disp('PASS READ_PARAMETERS')
end
%% READ BLOCK BOUNDARY DATA
function [BLK]=READ_BLOCK_BOUND(DIR)
EXT='*.txt';
file=dir([DIR,'/',EXT]);
[NBlock,~]=size(file);
BLK(1).NBlock=NBlock;
for NB=1:BLK(1).NBlock
  BLK(NB).fullname=fullfile(DIR,file(NB).name);
  tmp=load(BLK(NB).fullname);
  BLK(NB).name=file(NB).name;
  BLK(NB).LON=tmp(:,1);
  BLK(NB).LAT=tmp(:,2);
end
fprintf('READ BLOCK FILES : %4i \n',BLK(1).NBlock)
% figure('Name','BLOCK_BOUNDARY_LINE')
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    BLK(1).BOUND(NB1,NB2).LAT=[];
    BLK(1).BOUND(NB1,NB2).LON=[];
    LCa=inpolygon(BLK(NB1).LON,BLK(NB1).LAT,BLK(NB2).LON,BLK(NB2).LAT);
    Ca=find(LCa);
    if ~isempty(Ca)
      if and(LCa(1),LCa(end))
        Ca0=find(LCa~=true,1,'last')+1:length(LCa)-1;
        Ca1=1:find(LCa~=true,1,'first')-1;
        Ca=[Ca0 Ca1];
      end
      BLK(1).BOUND(NB1,NB2).LAT=BLK(NB1).LAT(Ca);
      BLK(1).BOUND(NB1,NB2).LON=BLK(NB1).LON(Ca);
      BLK(1).BOUND(NB1,NB2).BXYZ=conv2ell(BLK(1).BOUND(NB1,NB2).LAT,BLK(1).BOUND(NB1,NB2).LON);
      fprintf('BLOCK BOUNDARY : %2i %2i \n',NB1,NB2)
%       plot(BLK(1).BOUND(NB1,NB2).LON,BLK(1).BOUND(NB1,NB2).LAT)
%       hold on
    end
  end
end
end
%% Find shared point between blocks
function [Ca,LCa]=mach_bo(BLK,NB1,NB2)
LCa=inpolygon(BLK(NB1).LON,BLK(NB1).LAT,BLK(NB2).LON,BLK(NB2).LAT);
Ca=find(LCa);
if ~isempty(Ca)
  if and(LCa(1),LCa(end))
    Ca0=find(LCa~=true,1,'last')+1:length(LCa)-1;
    Ca1=1:find(LCa~=true,1,'first')-1;
    Ca=[Ca0 Ca1];
  end
end
end
%% Combine two blocks
function BLK=COMBINE_BOUND(BLK,NB1,NB2)
%% TODO: Need to combine BLOCK Interface file corresponding to combined BLOCK
for FN=1:BLK(1).NBlock
  namesplit=strsplit(BLK(FN).name,{'_','.'});
  Btmp(FN).NO=str2num(char(namesplit(1)));
  Btmp(FN).NAME=cell2mat(namesplit(2:end-1));
  if Btmp(FN).NO==NB1; NB1=FN; continue; end
  if Btmp(FN).NO==NB2; NB2=FN; end
end

[B.I(1).IND,B.I(1).AIND]=mach_bo(BLK,NB1,NB2);
[B.I(2).IND,B.I(2).AIND]=mach_bo(BLK,NB2,NB1);
B.I(1).NAIND=~B.I(1).AIND;
B.I(2).NAIND=~B.I(2).AIND;
if and(B.I(1).AIND(1),B.I(1).AIND(end))
  B.I(1).NAIND(find(B.I(1).NAIND==true,1,'first')-1)=true;
  B.I(1).NAIND(find(B.I(1).NAIND==true,1,'last')+1)=true;
  BLK(NB1).NLON=BLK(NB1).LON(B.I(1).NAIND);
  BLK(NB1).NLAT=BLK(NB1).LAT(B.I(1).NAIND);
else
  NCa0=find(B.I(1).AIND==true,1,'last'):length(B.I(1).AIND)-1;
  NCa1=1:find(B.I(1).AIND==true,1,'first');
  NCa=[NCa0 NCa1];
  BLK(NB1).NLON=BLK(NB1).LON(NCa);
  BLK(NB1).NLAT=BLK(NB1).LAT(NCa);
end
if and(B.I(2).AIND(1),B.I(2).AIND(end))
  BLK(NB2).NLON=BLK(NB2).LON(B.I(2).NAIND);
  BLK(NB2).NLAT=BLK(NB2).LAT(B.I(2).NAIND);
else
  NCa0=find(B.I(2).AIND==true,1,'last')+1:length(B.I(2).AIND)-1;
  NCa1=1:find(B.I(2).AIND==true,1,'first')-1;
  NCa=[NCa0 NCa1];
  BLK(NB2).NLON=BLK(NB2).LON(NCa);
  BLK(NB2).NLAT=BLK(NB2).LAT(NCa);
end
if isequal(BLK(NB1).LON(B.I(1).IND),BLK(NB2).LON(B.I(2).IND)) && isequal(BLK(NB1).LAT(B.I(1).IND),BLK(NB2).LAT(B.I(2).IND))  % Normal direction
  BLK(BLK(1).NBlock+1).LON=[BLK(NB1).NLON; flip(BLK(NB2).NLON); BLK(NB1).NLON(1)];
  BLK(BLK(1).NBlock+1).LAT=[BLK(NB1).NLAT; flip(BLK(NB2).NLAT); BLK(NB1).NLAT(1)];
else  % Reverse direction
  BLK(BLK(1).NBlock+1).LON=[BLK(NB1).NLON; BLK(NB2).NLON; BLK(NB1).NLON(1)];
  BLK(BLK(1).NBlock+1).LAT=[BLK(NB1).NLAT; BLK(NB2).NLAT; BLK(NB1).NLAT(1)];
end
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
