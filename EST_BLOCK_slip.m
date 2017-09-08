function EST_BLOCK_slip
% Estimate BLOCK MOTION
% Code by T.ITO 2016/03/02
%
warning('off','all')
INPUT.Parfile='./PARAMETER/parameter.txt';
INPUT.Optfile='./PARAMETER/opt_bound_par.txt';
% READ PARAMETER FOR MCMC Inversion 
[PRM]=READ_PARAMETERS(INPUT);
% READ OBSERVATION FILE
[OBS]=READ_OBS(PRM.FileOBS);
% READ BLOCK BOUNDARY FILE in DIRECTORY
[BLK,OBS]=READ_BLOCK_BOUND(PRM.DIRBlock,OBS);
% READ BLOCK INTERFACE BOUNDARY in DIRECTORY 
[BLK]=READ_BLOCK_INTERFACE(BLK,PRM);
% SHOW BLOCK BOUNDARY MAP
SHOW_BLOCK_BOUND(BLK)
% READ FIX EULER POLES
[POL,PRM]=READ_EULER_POLES(BLK,PRM);
% READ RIGID BLOCK BOUNDARY
[BLK,PRM]=READ_RIGID_BOUND(BLK,PRM);
% CALC. GREEN FUNCTION
[TRI,OBS]=GREEN_TRI(BLK,OBS);
% Combain to Green function
[D,G]=COMB_GREEN(BLK,OBS,TRI);
% CALC. ABIC AND BLOCK MOTION
[BLK,OBS]=CALC_AIC(BLK,OBS);
% BLOCK MOTION BETWEEN TWO BLOCKS
[BLK,OBS]=Est_Motion_BLOCKS(BLK,OBS);
% MAKE FIGURES
MAKE_FIGS(BLK,OBS);
% CAL Markov chain Monte Calro
[CHA]=MH_MCMC(D,G,BLK,PRM,OBS,POL);
% MAKE FIGURES
%MAKE_FIG(CHA,BLK,OBS,PRM);
OUTPUT.DIR='./Result/';
WRITE_CHA(CHA,BLK,TRI,PRM,OUTPUT,D,G)
%
end
%% WIRTE OUTPUT FILE
function WRITE_CHA(CHA,BLK,TRI,PRM,OUTPUT,D,G)
%
for DN=1:Inf
  DDIR=['Test_',num2str(DN,'%02i')];
  ADIR=fullfile(OUTPUT.DIR,DDIR);
  EXID=exist(ADIR);
  if EXID~=7; mkdir(ADIR); break; end
end
% 
fprintf('Write OUTPUT FILE: %s \n',OUTPUT.DIR)
dlmwrite(fullfile(ADIR,'Mp.txt'),single(CHA.Mp));
dlmwrite(fullfile(ADIR,'Mc.txt'),single(CHA.Mc));
dlmwrite(fullfile(ADIR,'La.txt'),single(CHA.La));
%
save(fullfile(ADIR,'CHA.mat'),'CHA','-v7.3')
save(fullfile(ADIR,'BLK.mat'),'BLK','-v7.3')
save(fullfile(ADIR,'TRI.mat'),'TRI','-v7.3')
save(fullfile(ADIR,'PRM.mat'),'PRM','-v7.3')
save(fullfile(ADIR,'GRN.mat'),'D','G','-v7.3')
% 
movefile('./Result/CHA_test*.mat',ADIR)
% 
end
%% UNIFORM MESH GENERATION
function [p,t]=mesh2D_uni(bou,int_bo,p_fix)
% bo(lon,lat)   : boundary 
% int_bo        : spaceing distance (km)
% pfix(lon,lat) : fixed points
dptol=0.01; ttol=0.01; delt=0.01; deps=0.01*int_bo; 
%
fprintf('Optimal interval %5.1f \n',int_bo)
ALON=mean(bou(:,1));ALAT=mean(bou(:,2));
[bo(:,2),bo(:,1)]=PLTXY(bou(:,2),bou(:,1),ALAT,ALON);
[pfix(:,2),pfix(:,1)]=PLTXY(p_fix(:,2),p_fix(:,1),ALAT,ALON);
[x,y]=meshgrid(min(bo(:,1)):0.7*int_bo:max(bo(:,1)),min(bo(:,2)):0.7*int_bo:max(bo(:,2)));
x(2:2:end,:)=x(2:2:end,:)+int_bo/2;
p=[x(:),y(:)];
p=p(dist_bo(p,bo)<deps,:);
if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end 
pfix=unique(pfix,'rows'); nfix=size(pfix,1);
p=[pfix; p];
np=size(p,1);
count=0;
pold=inf;
while 1
  count=count+1;
  if max(sqrt(sum((p-pold).^2,2))/int_bo)>ttol
    pold=p;
    t=delaunayn(p);
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
    t=t(dist_bo(pmid,bo)<-deps,:);
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
    bars=unique(sort(bars,2),'rows');
  end
  barvec=p(bars(:,1),:)-p(bars(:,2),:);
  L=sqrt(sum(barvec.^2,2));
  F=max(int_bo-L,0);
  Fvec=F./L*[1,1].*barvec;
  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],np,2));
  Ftot(1:size(pfix,1),:)=0;
  p=p+delt*Ftot;
  d=dist_bo(p,bo); ix=d>0; 
  if sum(ix)>1
    dgx=(dist_bo([p(ix,1)+deps,p(ix,2)],bo)-d(ix))/deps;
    dgy=(dist_bo([p(ix,1),p(ix,2)+deps],bo)-d(ix))/deps;
    dXY=d(ix)./(dgx.^2+dgy.^2);
    p(ix,:)=p(ix,:)-[dXY dXY].*[dgx dgy];
  end
  max_change=max(sqrt(sum(delt*Ftot(d<-deps,:).^2,2))/int_bo);
  if mod(count,100)==0
    fprintf('Count %4i Maxchange %4.1f Out side %3i Ave L %5.1f\n',count,100*max_change,sum(d>0),mean(L))
  end
  if abs(mean(L)-int_bo)/int_bo > 0.05 && mod(count,5)==0
    [~,ind_min]=min(L);
    p(setdiff(reshape(bars(ind_min,:),[],1),1:nfix),:)=[];
    pold=inf;  np=size(p,1);
    continue;
  end
  if count > 100*np | (abs(mean(L)-int_bo)/int_bo < 0.05 & max_change < dptol)
    fprintf('Count %4i Maxchange %4.1f Out side %3i Ave L %5.1f\n',count,100*max_change,sum(d>0),mean(L))
    break;
  end
end
[lat,lon]=XYTPL(p(:,2),p(:,1),ALAT,ALON);
p=[lon lat];
end
%====================================================
function d=dist_bo(po,bo)
% Distance of bo (line) and po (point)  
npo=length(po);
nbo=length(bo);
d=inf(npo,1);
for n=1:npo
  [~,ind]=sort(sqrt((bo(:,1)-po(n,1)).^2+(bo(:,2)-po(n,2)).^2));
  in1=ind(1);
  if in1 == 1; in2=2; in3=nbo;
  elseif in1 == nbo; in2=1; in3=nbo-1;
  else
    in2=in1+1; in3=in2-1;
  end
  if vec2ang(po(n,:),bo(in1,:),bo(in2,:)) < vec2ang(po(n,:),bo(in1,:),bo(in3,:))
    in2=in3;
  end
  a=bo(in2,2)-bo(in1,2);b=bo(in2,1)-bo(in1,1);
  d(n)=abs(a.*po(n,1)-b.*po(n,2)-a.*bo(in1,1)+b.*bo(in1,2))./sqrt(a.^2+b.^2);
end
d=(-1).^(inpolygon(po(:,1),po(:,2),bo(:,1),bo(:,2))).*d;
end
%====================================================
function ang=vec2ang(a,b1,b2)
b1b2=sqrt(sum((b2-b1).^2));
Lb1a=sqrt(sum((b1-a).^2));
Lb2a=sqrt(sum((b2-a).^2));
ang=min([dot(a-b1,b2-b1)/(b1b2.*Lb1a);dot(a-b2,b1-b2)/(b1b2.*Lb2a)]);
end
%====================================================
%% SHOW BLOCK BOUNDARY MAP
function SHOW_BLOCK_BOUND(BLK)
%
minlat=[];minlon=[];maxlat=[];maxlon=[];
minlatc=[];minlonc=[];maxlatc=[];maxlonc=[];
for NB=1:BLK(1).NBlock
  minlatc=min([(min(BLK(NB).LAT)+max(BLK(NB).LAT))/2;minlatc]); 
  minlonc=min([(min(BLK(NB).LON)+max(BLK(NB).LON))/2;minlonc]); 
  maxlatc=max([(min(BLK(NB).LAT)+max(BLK(NB).LAT))/2;maxlatc]); 
  maxlonc=max([(min(BLK(NB).LON)+max(BLK(NB).LON))/2;maxlonc]); 
  minlat=min([BLK(NB).LAT;minlat]); 
  minlon=min([BLK(NB).LON;minlon]); 
  maxlat=max([BLK(NB).LAT;maxlat]); 
  maxlon=max([BLK(NB).LON;maxlon]); 
end
%
figure(310); clf(310)
%figure(310,'Name','BLOCK_AND_BOUNDARY_MAP WIDE')
h = worldmap([minlat,maxlat],[minlon,maxlon]);
getm(h, 'MapProjection');
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
for NB=1:BLK(1).NBlock
  hold on
  plotm(BLK(NB).LAT,BLK(NB).LON,'red')
  hold on
  textm(mean(BLK(NB).LAT),mean(BLK(NB).LON),int2str(NB))
end
drawnow
%
figure(320); clf(320)
%figure('Name','BLOCK_AND_BOUNDARY_MAP REGIONAL')
h = worldmap([minlatc,maxlatc],[minlonc,maxlonc]);
getm(h, 'MapProjection');
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
for NB1=1:BLK(1).NBlock
  hold on
  plotm(BLK(NB1).LAT,BLK(NB1).LON,'red')
  hold on
  textm(mean(BLK(NB).LAT),mean(BLK(NB).LON),int2str(NB))
  for NB2=NB1+1:BLK(1).NBlock
    if ~isempty(BLK(1).BOUND(NB1,NB2).LAT)
      hold on
      plotm(BLK(1).BOUND(NB1,NB2).LAT,BLK(1).BOUND(NB1,NB2).LON,'o')
    end
  end
end
drawnow
%
figure(330); clf(330)
%figure('Name','BLOCK_AND_BOUNDARY_GEOMETRY')
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
    if NF~=0
      patch(BLK(1).BOUND(NB1,NB2).blon',BLK(1).BOUND(NB1,NB2).blat',BLK(1).BOUND(NB1,NB2).bdep',BLK(1).BOUND(NB1,NB2).bdep');
      hold on
    end
  end
end
drawnow 
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
FilePole=fscanf(Fid,'%s \n',[1,1]);
PRM.FilePole=fullfile(PRM.HOME_D,FilePole);
[~]=fgetl(Fid);
FileRigb=fscanf(Fid,'%s \n',[1,1]);
PRM.FileRigb=fullfile(PRM.HOME_D,FileRigb);
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
fprintf('File fixed epole   : %s \n',PRM.FilePole) 
fprintf('File Rigid boundary: %s \n',PRM.FileRigb) 
fprintf('GPUdev (CPU:99)    : %i \n',PRM.GPU) 
fprintf('ITR(Max_Nitr)      : %i \n',PRM.ITR) 
fprintf('CHA(Chain)         : %i \n',PRM.CHA) 
fprintf('KEP(KEEP)          : %i \n',PRM.KEP) 
fprintf('RWD(Walk_dis)      : %4.2f \n',PRM.RWD) 
fprintf('==================\n') 
%====================================================
disp('PASS READ_PARAMETERS')
end
%% MAKE MATRIX
function [D,G]=COMB_GREEN(BLK,OBS,TRI)
% Coded by Takeo Ito 2017/01/02 (ver 1.1)
% pole unit is mm
NOBS=length(OBS(1).EVEC);
TMP.OBS(1:3:3*NOBS)=OBS.EVEC;
TMP.OBS(2:3:3*NOBS)=OBS.NVEC;
TMP.OBS(3:3:3*NOBS)=OBS.HVEC;
TMP.ERR(1:3:3*NOBS)=OBS.EERR;
TMP.ERR(2:3:3*NOBS)=OBS.NERR;
TMP.ERR(3:3:3*NOBS)=OBS.HERR;
%
D(1).IND=find(TMP.ERR~=0)';
D(1).OBS=TMP.OBS(D(1).IND)';
D(1).ERR=TMP.ERR(D(1).IND)';
D(1).MID=[];
D(1).OBSID=zeros(3*NOBS,BLK(1).NBlock);
D(1).TRA=zeros(TRI.TNF,BLK(1).NBlock);
D(1).CFID=false(3*TRI(1).TNF,1);
D(1).CFDIPID=true(3*TRI(1).TNF,1);
D(1).CNT=0;
%
% (G(1).C * (( G(1).T * ( G(1).B1 - G(1).B2 ) * Mp)*Mc ) + G(1).P * Mp
%
G(1).T =zeros(3*TRI(1).TNF,2.*TRI(1).TNF);
G(1).Tt=zeros(3*TRI(1).TNF,2.*TRI(1).TNF);
G(1).B =zeros(2*TRI(1).TNF,3.*BLK(1).NBlock);
G(1).B1=zeros(3*TRI(1).TNF,3.*BLK(1).NBlock);
G(1).B2=zeros(3*TRI(1).TNF,3.*BLK(1).NBlock);
TMP.P=zeros(3*NOBS,3.*BLK(1).NBlock);
%
MC=1;
MT=1;
MR=1;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(TRI(1).BOUND(NB1,NB2).clon,2);
    if NF~=0
      D(1).CNT=D(1).CNT+1;
      D(1).mID=zeros(BLK(1).NB,1);
      D(1).mID(MR:MR+NF-1)=1;
      D(1).MID=[D(1).MID D(1).mID];
      if BLK(1).BOUND(NB1,NB2).type==5;D(1).CFID(MC:MC+3*NF-1)=true;end
      D(1).CFDIPID(MC+NF:MC+2*NF-1,1)=false;
      D(1).CFDIPID(MC+NF:MC+2*NF-1,1)=false;
% Need Project to direction of relative plate motion estimated from Pole
%       TMP.LD=sqrt(TRI(1).BOUND(NB1,NB2).DP(:,1).^2+TRI(1).BOUND(NB1,NB2).DP(:,2).^2);
%       TMP.LD(TMP.LD==0)=1;
      TMP.C(1:3*NOBS,MC     :MC+  NF-1)=TRI(1).BOUND(NB1,NB2).GSTR;
      TMP.C(1:3*NOBS,MC+  NF:MC+2*NF-1)=TRI(1).BOUND(NB1,NB2).GDIP;
      TMP.C(1:3*NOBS,MC+2*NF:MC+3*NF-1)=TRI(1).BOUND(NB1,NB2).GTNS;
      G(1).T(MC   :MC+  NF-1,MT   :MT+  NF-1)=diag(TRI(1).BOUND(NB1,NB2).ST(:,1));
      G(1).T(MC+NF:MC+2*NF-1,MT   :MT+  NF-1)=diag(TRI(1).BOUND(NB1,NB2).DP(:,1));
      G(1).T(MC   :MC+  NF-1,MT+NF:MT+2*NF-1)=diag(TRI(1).BOUND(NB1,NB2).ST(:,2));
      G(1).T(MC+NF:MC+2*NF-1,MT+NF:MT+2*NF-1)=diag(TRI(1).BOUND(NB1,NB2).DP(:,2));
      G(1).B(MT   :MT+  NF-1,3*NB1-2)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3));
      G(1).B(MT   :MT+  NF-1,3*NB1-1)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3));
      G(1).B(MT   :MT+  NF-1,3*NB1  )=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2)...
                                      +TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B(MT+NF:MT+2*NF-1,3*NB1-2)=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3)...
                                      +TRI(1).BOUND(NB1,NB2).OXYZ(:,6).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2));
      G(1).B(MT+NF:MT+2*NF-1,3*NB1-1)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3)...
                                      -TRI(1).BOUND(NB1,NB2).OXYZ(:,6).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B(MT+NF:MT+2*NF-1,3*NB1  )=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2)...
                                      -TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));

%       G(1).B(MT   :MT+  NF-1,3*NB1-2)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B(MT   :MT+  NF-1,3*NB1-2);
%       G(1).B(MT   :MT+  NF-1,3*NB1-1)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B(MT   :MT+  NF-1,3*NB1-1);
%       G(1).B(MT   :MT+  NF-1,3*NB1  )=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B(MT   :MT+  NF-1,3*NB1  );
%       G(1).B(MT+NF:MT+2*NF-1,3*NB1-2)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B(MT+NF:MT+2*NF-1,3*NB1-2);
%       G(1).B(MT+NF:MT+2*NF-1,3*NB1-1)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B(MT+NF:MT+2*NF-1,3*NB1-1);
%       G(1).B(MT+NF:MT+2*NF-1,3*NB1  )=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B(MT+NF:MT+2*NF-1,3*NB1  );
      
      G(1).B(MT   :MT+  NF-1,3*NB2-2)=-G(1).B(MT   :MT+  NF-1,3*NB1-2);
      G(1).B(MT   :MT+  NF-1,3*NB2-1)=-G(1).B(MT   :MT+  NF-1,3*NB1-1);
      G(1).B(MT   :MT+  NF-1,3*NB2  )=-G(1).B(MT   :MT+  NF-1,3*NB1  );
      G(1).B(MT+NF:MT+2*NF-1,3*NB2-2)=-G(1).B(MT+NF:MT+2*NF-1,3*NB1-2);
      G(1).B(MT+NF:MT+2*NF-1,3*NB2-1)=-G(1).B(MT+NF:MT+2*NF-1,3*NB1-1);
      G(1).B(MT+NF:MT+2*NF-1,3*NB2  )=-G(1).B(MT+NF:MT+2*NF-1,3*NB1  );           
% Transrate Matrix for dE and dN
      G(1).Tt(MC+NF:MC+2*NF-1,MT   :MT+  NF-1)=diag(TRI(1).BOUND(NB1,NB2).ST(:,1));
      G(1).Tt(MC   :MC+  NF-1,MT   :MT+  NF-1)=diag(TRI(1).BOUND(NB1,NB2).DP(:,1));
      G(1).Tt(MC+NF:MC+2*NF-1,MT+NF:MT+2*NF-1)=diag(TRI(1).BOUND(NB1,NB2).ST(:,2));
      G(1).Tt(MC   :MC+  NF-1,MT+NF:MT+2*NF-1)=diag(TRI(1).BOUND(NB1,NB2).DP(:,2));
%       
      G(1).B1(MC     :MC+  NF-1,3*NB1-2)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3));
      G(1).B1(MC     :MC+  NF-1,3*NB1-1)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3));
      G(1).B1(MC     :MC+  NF-1,3*NB1  )=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2)...
                                      +TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-2)=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3)...
                                      +TRI(1).BOUND(NB1,NB2).OXYZ(:,6).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2));
      G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-1)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3)...
                                      -TRI(1).BOUND(NB1,NB2).OXYZ(:,6).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B1(MC+  NF:MC+2*NF-1,3*NB1  )=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2)...
                                      -TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B1(MC+2*NF:MC+3*NF-1,3*NB1-2)=0;
      G(1).B1(MC+2*NF:MC+3*NF-1,3*NB1-1)=0;
      G(1).B1(MC+2*NF:MC+3*NF-1,3*NB1  )=0;
%       G(1).B1(MC     :MC+  NF-1,3*NB1-2)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B1(MC     :MC+  NF-1,3*NB1-2);
%       G(1).B1(MC     :MC+  NF-1,3*NB1-1)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B1(MC     :MC+  NF-1,3*NB1-1);
%       G(1).B1(MC     :MC+  NF-1,3*NB1  )=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B1(MC     :MC+  NF-1,3*NB1  );
%       G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-2)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-2);
%       G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-1)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-1);
%       G(1).B1(MC+  NF:MC+2*NF-1,3*NB1  )=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B1(MC+  NF:MC+2*NF-1,3*NB1  );
%       G(1).B1(MC+2*NF:MC+3*NF-1,3*NB1-2)=0;
%       G(1).B1(MC+2*NF:MC+3*NF-1,3*NB1-1)=0;
%       G(1).B1(MC+2*NF:MC+3*NF-1,3*NB1  )=0;
      G(1).B1(MC     :MC+  NF-1,3*NB2-2)=-G(1).B1(MC     :MC+  NF-1,3*NB1-2);
      G(1).B1(MC     :MC+  NF-1,3*NB2-1)=-G(1).B1(MC     :MC+  NF-1,3*NB1-1);
      G(1).B1(MC     :MC+  NF-1,3*NB2  )=-G(1).B1(MC     :MC+  NF-1,3*NB1  );
      G(1).B1(MC+  NF:MC+2*NF-1,3*NB2-2)=-G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-2);
      G(1).B1(MC+  NF:MC+2*NF-1,3*NB2-1)=-G(1).B1(MC+  NF:MC+2*NF-1,3*NB1-1);
      G(1).B1(MC+  NF:MC+2*NF-1,3*NB2  )=-G(1).B1(MC+  NF:MC+2*NF-1,3*NB1  );
      G(1).B1(MC+2*NF:MC+3*NF-1,3*NB2-2)=0;
      G(1).B1(MC+2*NF:MC+3*NF-1,3*NB2-1)=0;
      G(1).B1(MC+2*NF:MC+3*NF-1,3*NB2  )=0;
% 
      G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-2)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3));
      G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-1)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3));
      G(1).B2(MC+  NF:MC+2*NF-1,3*NB1  )=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2)...
                                      +TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B2(MC     :MC+  NF-1,3*NB1-2)=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3)...
                                      +TRI(1).BOUND(NB1,NB2).OXYZ(:,6).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2));
      G(1).B2(MC     :MC+  NF-1,3*NB1-1)=-1.*(-TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,3)...
                                      -TRI(1).BOUND(NB1,NB2).OXYZ(:,6).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B2(MC     :MC+  NF-1,3*NB1  )=-1.*( TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,7).*TRI(1).BOUND(NB1,NB2).OXYZ(:,2)...
                                      -TRI(1).BOUND(NB1,NB2).OXYZ(:,4).*TRI(1).BOUND(NB1,NB2).OXYZ(:,5).*TRI(1).BOUND(NB1,NB2).OXYZ(:,1));
      G(1).B2(MC+2*NF:MC+3*NF-1,3*NB1-2)=0;
      G(1).B2(MC+2*NF:MC+3*NF-1,3*NB1-1)=0;
      G(1).B2(MC+2*NF:MC+3*NF-1,3*NB1  )=0;
%       G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-2)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-2);
%       G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-1)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-1);
%       G(1).B2(MC+  NF:MC+2*NF-1,3*NB1  )=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B2(MC+  NF:MC+2*NF-1,3*NB1  );
%       G(1).B2(MC     :MC+  NF-1,3*NB1-2)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B2(MC     :MC+  NF-1,3*NB1-2);
%       G(1).B2(MC     :MC+  NF-1,3*NB1-1)=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B2(MC     :MC+  NF-1,3*NB1-1);
%       G(1).B2(MC     :MC+  NF-1,3*NB1  )=TRI(1).BOUND(NB1,NB2).SDTINV.*G(1).B2(MC     :MC+  NF-1,3*NB1  );
%       G(1).B2(MC+2*NF:MC+3*NF-1,3*NB1-2)=0;
%       G(1).B2(MC+2*NF:MC+3*NF-1,3*NB1-1)=0;
%       G(1).B2(MC+2*NF:MC+3*NF-1,3*NB1  )=0;
      G(1).B2(MC+  NF:MC+2*NF-1,3*NB2-2)=-G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-2);
      G(1).B2(MC+  NF:MC+2*NF-1,3*NB2-1)=-G(1).B2(MC+  NF:MC+2*NF-1,3*NB1-1);
      G(1).B2(MC+  NF:MC+2*NF-1,3*NB2  )=-G(1).B2(MC+  NF:MC+2*NF-1,3*NB1  );
      G(1).B2(MC     :MC+  NF-1,3*NB2-2)=-G(1).B2(MC     :MC+  NF-1,3*NB1-2);
      G(1).B2(MC     :MC+  NF-1,3*NB2-1)=-G(1).B2(MC     :MC+  NF-1,3*NB1-1);
      G(1).B2(MC     :MC+  NF-1,3*NB2  )=-G(1).B2(MC     :MC+  NF-1,3*NB1  );
      G(1).B2(MC+2*NF:MC+3*NF-1,3*NB2-2)=0;
      G(1).B2(MC+2*NF:MC+3*NF-1,3*NB2-1)=0;
      G(1).B2(MC+2*NF:MC+3*NF-1,3*NB2  )=0;
%       
      MC=MC+3*NF;
      MT=MT+2*NF;
      MR=MR+  NF;
    end
  end
%   
  IND=OBS(1).ABLK==NB1;  
  NIND=[zeros(size(IND)),IND,zeros(size(IND))]; NIND=logical(reshape(NIND',3*NOBS,1));
  EIND=[IND,zeros(size(IND)),zeros(size(IND))]; EIND=logical(reshape(EIND',3*NOBS,1));
  TMP.P(EIND,3*NB1-2)=-OBS(1).AXYZ(IND,7).*OBS(1).AXYZ(IND,3);
  TMP.P(EIND,3*NB1-1)=-OBS(1).AXYZ(IND,5).*OBS(1).AXYZ(IND,3);
  TMP.P(EIND,3*NB1  )= OBS(1).AXYZ(IND,5).*OBS(1).AXYZ(IND,2)                    +OBS(1).AXYZ(IND,7).*OBS(1).AXYZ(IND,1);
  TMP.P(NIND,3*NB1-2)= OBS(1).AXYZ(IND,4).*OBS(1).AXYZ(IND,5).*OBS(1).AXYZ(IND,3)+OBS(1).AXYZ(IND,6).*OBS(1).AXYZ(IND,2);
  TMP.P(NIND,3*NB1-1)=-OBS(1).AXYZ(IND,4).*OBS(1).AXYZ(IND,7).*OBS(1).AXYZ(IND,3)-OBS(1).AXYZ(IND,6).*OBS(1).AXYZ(IND,1);
  TMP.P(NIND,3*NB1  )= OBS(1).AXYZ(IND,4).*OBS(1).AXYZ(IND,7).*OBS(1).AXYZ(IND,2)-OBS(1).AXYZ(IND,4).*OBS(1).AXYZ(IND,5).*OBS(1).AXYZ(IND,1);
end
% 
for BL=1:BLK(1).NBlock
  D(1).obsid=zeros(1,NOBS);
  D(1).obsid(1,OBS(1).ABLK==BL)=true;
  D(1).OBSID(:,BL)=reshape(repmat(D(1).obsid,3,1),3*NOBS,1);
  MC=1;
  MT=1;
  MR=1;
  rigblkID=BLK(1).RGPAIR(:,1)==BL;
  FLAG1=sum(rigblkID);
  for NB1=1:BLK(1).NBlock
    for NB2=NB1+1:BLK(1).NBlock
      NF=size(TRI(1).BOUND(NB1,NB2).clon,2);
      if NF~=0
        if FLAG1==0
          if NB2==BL && NB1<NB2
            D(1).TRA(MC:MC+3*NF-1,BL)=-1;
          else
            D(1).TRA(MC:MC+3*NF-1,BL)= 1;
          end
        else
          FLAG2=0;
          RGPAIRtmp=BLK(1).RGPAIR(rigblkID,:);
          for RB=1:size(RGPAIRtmp,1)
            RGPAIRID=ismember([NB1 NB2],RGPAIRtmp(RB,2:3));
            ISPAIR=sum(RGPAIRID);
            if ISPAIR==2;FLAG2=1;break;end
          end
          if FLAG2==1
            D(1).TRA(MC:MC+3*NF-1,BL)=0;
          elseif NB2==BL && NB1<NB2
            D(1).TRA(MC:MC+3*NF-1,BL)=-1;
          else
            D(1).TRA(MC:MC+3*NF-1,BL)= 1;
          end
        end
        MC=MC+3*NF;
        MT=MT+2*NF;
        MR=MR+  NF;
      end
    end
  end
end
% 
G(1).C  =TMP.C(D(1).IND,:);
G(1).P  =TMP.P(D(1).IND,:);
G(1).T  =   sparse(G(1).T);
G(1).TB =    G(1).T*G(1).B;
G(1).Tt =  sparse(G(1).Tt);
G(1).TtB=   G(1).Tt*G(1).B;
D(1).MID=logical(repmat(D(1).MID,3,1));
end
%% Markov chain Monte Calro
function [CHA]=MH_MCMC(D,G,BLK,PRM,OBS,POL)
% Markov chain Monte Calro
RR=(D(1).OBS./D(1).ERR)'*(D(1).OBS./D(1).ERR);
fprintf('Residual=%9.3f \n',RR);
%
% if PRM.GPU==99
%   precision='double';
% else
%   precision='single';
% end
precision='double';
RWD=PRM.RWD;
Mc.INT=1e-2;
Mp.INT=1e-10;
La.INT=1e+1;
Mc.N=BLK(1).NB;
Mp.N=3.*BLK(1).NBlock;
La.N=1;
Mc.STD=Mc.INT.*ones(Mc.N,1,precision);
Mp.STD=Mp.INT.*ones(Mp.N,1,precision);
La.STD=La.INT.*ones(La.N,1,precision);
Mc.OLD=   -0.5+rand(Mc.N,1,precision);
Mp.OLD= double(BLK(1).POLE);
La.OLD= zeros(La.N,1,precision);
CHA.Mc= zeros(Mc.N,PRM.KEP,precision);
CHA.Mp= zeros(Mp.N,PRM.KEP,precision);
CHA.La= zeros(La.N,PRM.KEP,precision);
% Set FIX POLES if POL.FIXflag=1
% MpScale=Mp.INT.*ones(Mp.N,1,precision);
if POL.FIXflag==1
  Mp.OLD(POL.ID)=0; Mp.OLD=Mp.OLD+POL.FIXw;
  Mp.STD(POL.ID)=0;
%   MpScale(POL.ID)=0;
end
%
RES.OLD=inf(1,1,precision);
PRI.OLD=inf(1,1,precision);
McScale=0.13;
MpScale=(1.3E-9).*ones(Mp.N,1,precision).*~POL.ID;
% McScale=0.05;
% MpScale=3E-10.*ones(Mp.N,1,precision).*~POL.ID;
LO_Mc=-1;
UP_Mc=1;
% GPU Initialize 
if PRM.GPU~=99
  g=gpuDevice(PRM.GPU);
  reset(g);
  g_men=g.TotalMemory;
  r_men=(Mc.N+Mp.N+La.N).*(PRM.KEP+2).*4;
  res_m=g_men-r_men;
  fprintf('USE GPU Max Chain=%4d Nitr=%2d Mc=%4d Mp=%3d res_Memory(GB)=%6.3f\n',...
           PRM.CHA,PRM.ITR,Mc.N,Mp.N,res_m./(1024.*1024.*1024));      
%   CHA.Mc=gpuArray(CHA.Mc);
%   CHA.Mp=gpuArray(CHA.Mp);
%   CHA.La=gpuArray(CHA.La);
  Mc.STD=gpuArray(Mc.STD);
  Mp.STD=gpuArray(Mp.STD);
  La.STD=gpuArray(La.STD);
  Mc.OLD=gpuArray(Mc.OLD);
  Mp.OLD=gpuArray(Mp.OLD);
  La.OLD=gpuArray(La.OLD);
  D(1).OBS=gpuArray(D(1).OBS);
  D(1).ERR=gpuArray(D(1).ERR);
  G.TB=gpuArray(G.TB);
  G.C=gpuArray(G.C);
  G.P=gpuArray(G.P);
  McScale=gpuArray(McScale);
  MpScale=gpuArray(MpScale);
  LO_Mc=gpuArray(LO_Mc);
  UP_Mc=gpuArray(UP_Mc);
  RWD=gpuArray(RWD);
  G.B1=gpuArray(G.B1);
  G.B2=gpuArray(G.B2);
  G.TtB=gpuArray(G.TtB);
  D.TRA=gpuArray(D.TRA);
  D.OBSID=gpuArray(D.OBSID);
else
  fprintf('USE CPU Max Chain=%4d Nitr=%2d Mc=%4d Mp=%3d \n',...
            PRM.CHA,PRM.ITR,Mc.N,Mp.N);
end
%
RT=0;
COUNT=0;
%
while not(COUNT==20)
  RT  =RT+1;
  NACC=0;tic
  if PRM.GPU~=99
    logU=log(rand(PRM.CHA,1,precision,'gpuArray'));
    rMc = rand(Mc.N,PRM.CHA,precision,'gpuArray')-0.5;
    rMp = rand(Mp.N,PRM.CHA,precision,'gpuArray')-0.5;
    rLa = rand(La.N,PRM.CHA,precision,'gpuArray')-0.5;
  else
    logU=log(rand(PRM.CHA,1,precision));
    rMc =rand(Mc.N,PRM.CHA,precision)-0.5;
    rMp =rand(Mp.N,PRM.CHA,precision)-0.5;
    rLa =rand(La.N,PRM.CHA,precision)-0.5;
  end
  for iT=1:PRM.CHA
% SAMPLE SECTION
%     McUp=min(UP_Mc,Mc.OLD+0.5.*RWD.*Mc.STD);
%     McLo=max(LO_Mc,Mc.OLD-0.5.*RWD.*Mc.STD);
%     Mc.SMP=McLo+(McUp-McLo).*rMc(:,iT);
    McTMP=Mc.OLD+0.5.*RWD.*McScale.*rMc(:,iT);
    McREJID=McTMP>UP_Mc | McTMP<LO_Mc;
    McTMP(McREJID)=Mc.OLD(McREJID);
    Mc.SMP=McTMP;
%     Mc.SMP=max(min(McTMP,UP_Mc),LO_Mc);
%     Mp.SMP=Mp.OLD+RWD.*Mp.STD.*rMp(:,iT);
    Mp.SMP=Mp.OLD+RWD.*MpScale.*rMp(:,iT);
    La.SMP=La.OLD+RWD.*La.STD.*rLa(:,iT);
% MAKE Mc.SMPMAT
    Mc.SMPMAT=repmat(Mc.SMP,3,D.CNT);
    Mc.SMPMAT=Mc.SMPMAT(D.MID);
% Calc GPU memory free capacity
    if PRM.GPU~=99
      Byte1=whos('G');
      Byte2=whos('Mp');
      b=waitGPU(Byte1.bytes+Byte2.bytes);
    end
% Calc Correction factor of subducting rate for DIP direction.
% VE^2+VN^2 = Vst^2+(CF*Vdp)^2 <=> (G.B1*Mp.SMP).^2+(G.B2*Mp.SMP).^2 = (G.TtB*Mp.SMP).^2+(CF*G.TB*Mp.SMP).^2
    CFsq=((G.B1*Mp.SMP).^2+(G.B2*Mp.SMP).^2-(G.TtB*Mp.SMP).^2)./((G.TB*Mp.SMP).^2);
    CFsq(CFsq<0)=0;
%     if sum(((G.B1*Mp.SMP).^2+(G.B2*Mp.SMP).^2-(G.TtB*Mp.SMP).^2)./((G.TB*Mp.SMP).^2)<0)~=0
%         keyboard
%     end
    CF=sqrt(CFsq);
    CF(or(D.CFDIPID,or(isnan(CF),D.CFID)))=1;
% CALC APRIORI AND RESIDUAL COUPLING RATE SECTION
    CAL.RIG=G.P*Mp.SMP;
    CAL.ela=G.C*(repmat((G.TB*Mp.SMP),1,BLK(1).NBlock).*D.TRA.*repmat(Mc.SMPMAT,1,BLK(1).NBlock).*repmat(CF,1,BLK(1).NBlock));
    CAL.ELA=sum(CAL.ela.*D.OBSID,2);
    CAL.SMP=CAL.RIG+CAL.ELA;
    if PRM.GPU~=99
      clear('CAL.RIG','CAL.ela','CAL,ELA','CF','CFsq');
    end
%   CAL.SMP=G.C*((G.TB*Mp.SMP).*Mc.SMPMAT)+G.P*Mp.SMP;       
%   CAL.SMP=G.P*Mp.SMP;
% CALC RESIDUAL SECTION
    RES.SMP=sum(((D(1).OBS-CAL.SMP)./D(1).ERR).^2,1);
% Mc is better Zero 
    PRI.SMP=sum(abs(Mc.SMP),1);   
%% MAKE Probably Density Function
% $$ PDF_{post}=\frac{\frac{1}{\sqrt{2\pi\exp(L)}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}}{2}}\exp{\frac{-M^{2}}{2\times\exp{L}}}{\frac{1}{\sqrt{2\pi\exp(L_{old})}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}_{old}}{2}}\exp{\frac{-M^{2}_{old}}{2\times\exp{L_{old}}}} $$%%
%  log(x(x>0));
%   q1 = logproppdf(x0,y);
%   q2 = logproppdf(y,x0);
% this is a generic formula.
%   rho = (q1+logpdf(y))-(q2+logpdf(x0));  
    Pdf = -0.5.*...
         ((RES.SMP+La.SMP+exp(-La.SMP).*PRI.SMP)...
         -(RES.OLD+La.OLD+exp(-La.OLD).*PRI.OLD));
%   Pdf = -0.5.*(RES.SMP-RES.OLD);
    ACC=Pdf > logU(iT);
    if ACC
      Mc.OLD  = Mc.SMP;
      Mp.OLD  = Mp.SMP;
      La.OLD  = La.SMP;
      RES.OLD = RES.SMP;
      PRI.OLD = PRI.SMP;
    end
% KEEP SECTION
    if iT >= PRM.CHA-PRM.KEP
      if PRM.GPU~=99
        CHA.Mc(:,iT-(PRM.CHA-PRM.KEP)+1)=gather(Mc.SMP);
        CHA.Mp(:,iT-(PRM.CHA-PRM.KEP)+1)=gather(Mp.SMP);
        CHA.La(:,iT-(PRM.CHA-PRM.KEP)+1)=gather(La.SMP);
      else
        CHA.Mc(:,iT-(PRM.CHA-PRM.KEP)+1)=Mc.SMP;
        CHA.Mp(:,iT-(PRM.CHA-PRM.KEP)+1)=Mp.SMP;
        CHA.La(:,iT-(PRM.CHA-PRM.KEP)+1)=La.SMP;
      end
      if ACC; NACC=NACC+1; end;
    end
  end
  COMPRESS_DATA(CHA,PRM,RT,NACC);
%
  CHA.AJR=NACC./PRM.CHA;
%
  Mc.STD=std(CHA.Mc,1,2);
  Mp.STD=std(CHA.Mp,1,2); 
  La.STD=std(CHA.La,1,2);
%
  fprintf('T=%3d Res=%6.3f Accept=%5.1f RWD=%5.2f Time=%5.1fsec\n',...
           RT,1-RES.OLD./RR,100*CHA.AJR,RWD,toc)
%
  for BK=1:BLK(1).NBlock
    [latp,lonp,ang]=xyzp2lla(CHA.Mp(3.*BK-2,:),CHA.Mp(3.*BK-1,:),CHA.Mp(3.*BK,:));
    fprintf('POLE OF BLOCK %2i = lat:%7.2f deg. lon:%8.2f deg. ang:%9.2e deg./m.y. \n',...
      BK,mean(latp),mean(lonp),mean(ang));
  end
  fprintf('Lamda = %7.2f \n',mean(CHA.La));
%
  if CHA.AJR > 0.24
    RWD=RWD*1.1;
    COUNT=0;
  elseif CHA.AJR < 0.22
    RWD=RWD*0.9;
    COUNT=0;
  else
    COUNT=COUNT+1;
  end
  CHA.SMP=CAL.SMP;
  % debug-----------
  Mpmean=mean(CHA.Mp,2);
  Mcmean=mean(CHA.Mc,2);
  Mcmeanrep=repmat(Mcmean,3,D.CNT);Mcmeanrep=Mcmeanrep(D.MID);
  CFsq=((G.B1*Mpmean).^2+(G.B2*Mpmean).^2-(G.TtB*Mpmean).^2)./((G.TB*Mpmean).^2);
  CFsq(CFsq<0)=0;
  CF=sqrt(CFsq);
  CF(or(D.CFDIPID,or(isnan(CF),D.CFID)))=1;
  VEC.RIG=G.P*Mpmean;
  VEC.ela=G.C*(repmat((G.TB*Mpmean),1,BLK(1).NBlock).*D.TRA.*repmat(Mcmeanrep,1,BLK(1).NBlock).*repmat(CF,1,BLK(1).NBlock));
  VEC.ELA=sum(VEC.ela.*D.OBSID,2);
%   vec.rel=G.C*((G.TB*poltmp).*CF);
  % debug-----------
  if PRM.GPU~=99
    cCHA.Mc=gather(CHA.Mc);
    cCHA.Mp=gather(CHA.Mp);
    cCHA.La=gather(CHA.La);
    cCHA.SMP=gather(CHA.SMP);
    MAKE_FIG(cCHA,BLK,OBS,RT,gather(LO_Mc),gather(UP_Mc),VEC)
  else
    MAKE_FIG(CHA,BLK,OBS,RT,LO_Mc,UP_Mc,VEC)
  end
  if RT > PRM.ITR; break; end;
end
if PRM.GPU~=99
  CHA.Mc=gather(CHA.Mc);
  CHA.Mp=gather(CHA.Mp);
  CHA.La=gather(CHA.La);
  CHA.SMP=gather(CHA.SMP);
end
CHA.Res=RES.SMP;
fprintf('RMS=: %8.3f\n',CHA.Res)
fprintf('=== FINISHED MH_MCMC ===\n')
end
%% Compress CHA sampling
function COMPRESS_DATA(CHA,PRM,ITR,NACC)
% 
% load('./Result_red/Test_07/CHA.mat'); % test
% 
CHA.Mc=single(CHA.Mc);
CHA.Mp=single(CHA.Mp);
% if PRM.GPU==99&&gpuDeviceCount==0
if PRM.GPU==99
  MEANMc=mean(CHA.Mc,2);
  MEANMp=mean(CHA.Mp,2);
  COVMc=cov(CHA.Mc');
  COVMp=cov(CHA.Mp');
else
  gCHA.Mc=gpuArray(CHA.Mc);
  gCHA.Mp=gpuArray(CHA.Mp);
  MEANMc=mean(gCHA.Mc,2);
  MEANMp=mean(gCHA.Mp,2);
  COVMc=cov(gCHA.Mc');
  COVMp=cov(gCHA.Mp');
  MEANMc=gather(MEANMc);
  MEANMp=gather(MEANMp);
  COVMc=gather(COVMc);
  COVMp=gather(COVMp);
end
% 
McMAX=max(CHA.Mc,[],2);
McMIN=min(CHA.Mc,[],2);
MpMAX=max(CHA.Mp,[],2);
MpMIN=min(CHA.Mp,[],2);
% 
Mcscale=100./(McMAX-McMIN);
McBASE=(CHA.Mc-McMIN).*Mcscale*2.55-128;
Mcint8=int8(McBASE);
Mpscale=100./(MpMAX-MpMIN);
MpBASE=(CHA.Mp-MpMIN).*Mpscale*2.55-128;
Mpint8=int8(MpBASE);
% 
binedge=int8(-128:127);
% 
for ii=1:size(Mcint8,1)
  cha.McCOMPRESS.NFLT(ii).Mcscale=Mcscale(ii);
  cha.McCOMPRESS.NFLT(ii).McMAX=McMAX(ii);
  cha.McCOMPRESS.NFLT(ii).McMIN=McMIN(ii);
  cha.McCOMPRESS.NFLT(ii).McHIST=histcounts(Mcint8(ii,:),binedge);
end
cha.McCOMPRESS.COVMc=COVMc;
cha.McCOMPRESS.MEANMc=MEANMc;
cha.McCOMPRESS.SMPMc=int8(McBASE);
% 
for ii=1:size(Mpint8,1)
  cha.MpCOMPRESS.NPOL(ii).Mpscale=Mpscale(ii);
  cha.MpCOMPRESS.NPOL(ii).MpMAX=MpMAX(ii);
  cha.MpCOMPRESS.NPOL(ii).MpMIN=MpMIN(ii);
  cha.MpCOMPRESS.NPOL(ii).MpHIST=histcounts(Mpint8(ii,:),binedge);
end
cha.MpCOMPRESS.COVMp=COVMp;
cha.MpCOMPRESS.MEANMp=MEANMp;
cha.MpCOMPRESS.SMPMp=int8(MpBASE);
% 
cha.AJR=NACC./PRM.CHA;
save(['./Result/CHA_test',num2str(ITR,'%03i'),'.mat'],'cha','-v7.3'); % test
% 
end
%% Show results for makeing FIGURES
function MAKE_FIG(CHA,BLK,OBS,RT,LO_Mc,UP_Mc,VEC)
% Color palette(POLAR)
red=[0:1/32:1 ones(1,32)]';
green=[0:1/32:1 1-1/32:-1/32:0]';
blue=[ones(1,32) 1:-1/32:0]';
rwb=[red green blue];
rw =rwb(33:end,:);
if LO_Mc==-1; cmap=rwb;
else; cmap=rw; end
% 
figure(100);clf(100)
% BUG to wait zero
NN=1;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
    if NF~=0
      patch(BLK(1).BOUND(NB1,NB2).blon',BLK(1).BOUND(NB1,NB2).blat',BLK(1).BOUND(NB1,NB2).bdep',mean(CHA.Mc(NN:NN+NF-1,:),2));
      NN=NN+NF;
      hold on
    end
  end
end
ax=gca;
ax.CLim=[LO_Mc UP_Mc];
colormap(cmap)
colorbar
%
figure(110);clf(110)
% BUG to wait zero
NN=1;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
    if NF~=0
      patch(BLK(1).BOUND(NB1,NB2).blon',BLK(1).BOUND(NB1,NB2).blat',BLK(1).BOUND(NB1,NB2).bdep',std(CHA.Mc(NN:NN+NF-1,:),0,2));
      NN=NN+NF;
      hold on
    end
  end
end
colormap(parula)
colorbar
%
figure(120);clf(120)
for NB=1:BLK(1).NBlock
  plot(BLK(NB).LON,BLK(NB).LAT,'red')
  hold on
  text(mean(BLK(NB).LON),mean(BLK(NB).LAT),int2str(NB))
  hold on
  [latp,lonp,~]=xyzp2lla(CHA.Mp(3.*NB-2,:),CHA.Mp(3.*NB-1,:),CHA.Mp(3.*NB,:));
  minlon=min(lonp); maxlon=max(lonp); 
  minlat=min(latp); maxlat=max(latp); 
  if maxlon-minlon < 0.5; binlon=[minlon maxlon]; else; binlon=minlon:0.5:maxlon; end  
  if maxlat-minlat < 0.5; binlat=[minlat maxlat]; else; binlat=minlat:0.5:maxlat; end  
  histogram2(lonp,latp,binlon,binlat,'Normalization','probability','FaceColor','flat')
  hold on
  text(double(mean(lonp)),double(mean(latp)),int2str(NB))
  hold on
end
quiver(OBS(1).ALON,OBS(1).ALAT,OBS(1).EVEC,OBS(1).NVEC,'green')
quiver(OBS(1).ALON,OBS(1).ALAT,CHA.SMP(1:3:end)',CHA.SMP(2:3:end)','blue')
colorbar
hold on
%
figure(130);clf(130)
quiver(OBS(1).ALON,OBS(1).ALAT,OBS(1).EVEC,OBS(1).NVEC,'green')
hold on
quiver(OBS(1).ALON,OBS(1).ALAT,CHA.SMP(1:3:end)',CHA.SMP(2:3:end)','blue')
hold on
axis([OBS(1).LONMIN-1,OBS(1).LONMAX+1,OBS(1).LATMIN-1,OBS(1).LATMAX+1]);
title(['Iteration Number: ',num2str(RT)]);
% 
% debug----------
figure(140);clf(140)
quiver(OBS(1).ALON,OBS(1).ALAT,VEC.RIG(1:3:end)',VEC.RIG(2:3:end)','k')
hold on
quiver(OBS(1).ALON,OBS(1).ALAT,VEC.ELA(1:3:end)',VEC.ELA(2:3:end)','r')
% hold on
% quiver(OBS(1).ALON,OBS(1).ALAT,vec.rel(1:3:end)',vec.rel(2:3:end)','m')
axis([OBS(1).LONMIN-1,OBS(1).LONMAX+1,OBS(1).LATMIN-1,OBS(1).LATMAX+1]);
title(['Iteration Number: ',num2str(RT)]);
% debug----------
drawnow
end
%% READ PLATE INTERFACE
function [BLK]=READ_BLOCK_INTERFACE(BLK,PRM)
% Coded by Takeo Ito 2016/12/21 (ver 1.0)
%
int_tri=50;
dep_limit=-100;
dep_limit_low=-20;
DIRBLK=PRM.DIRBlock_Interface;
BLK(1).NB=0;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    BLK(1).BOUND(NB1,NB2).type=1;
    pre_tri_f=fullfile(DIRBLK,['triB_',num2str(NB1),'_',num2str(NB2),'.txt']); 
    Fid=fopen(pre_tri_f,'r');
    if Fid >= 0
      fprintf('BLOCK INTERFACE: %2i  %2i \n',NB1,NB2)
      fprintf('READ INTERFACE TRI BOUDARY FILE : %s \n',pre_tri_f)
      NF=0;
      blon=zeros(1,3);blat=zeros(1,3);bdep=zeros(1,3);
      while 1
        NF=NF+1;
        loc_f=fscanf(Fid,'%f %f %f \n', [3 3]);
        [~] = fgetl(Fid);
        blon(NF,:)=loc_f(1,:);%Lon
        blat(NF,:)=loc_f(2,:);%Lat
        bdep(NF,:)=loc_f(3,:);%Hight
        tline = fgetl(Fid); if ~ischar(tline); break; end
      end
      fclose(Fid);
      BO_tri_f=fullfile(DIRBLK,['triBO_',num2str(NB1),'_',num2str(NB2),'.txt']); 
      Fid=fopen(BO_tri_f,'r');
      if Fid >= 0
        bound_blk=textscan(Fid,'%f%f'); fclose(Fid);     
        bound_blk=cell2mat(bound_blk);
        Bslon=mean(blon,2);
        Bslat=mean(blat,2);
        ID=inpolygon(Bslon,Bslat,bound_blk(:,1),bound_blk(:,2));
        blon=blon(ID,:);
        blat=blat(ID,:);
        bdep=bdep(ID,:);
      end
      BLK(1).BOUND(NB1,NB2).blon=blon;%Lon
      BLK(1).BOUND(NB1,NB2).blat=blat;%Lat
      BLK(1).BOUND(NB1,NB2).bdep=bdep;%Hight      
    else
      sub_f=fullfile(DIRBLK,['B_',num2str(NB1),'_',num2str(NB2),'.txt']);
      Fid=fopen(sub_f,'r');
      if Fid >= 0
        fprintf('BLOCK INTERFACE: %2i  %2i \n',NB1,NB2)
        fprintf('READ INTERFACE BOUDARY SHAPE FILE : %s \n',sub_f)
        dep_blk=textscan(Fid,'%f%f%f'); fclose(Fid);
        dep_blk=cell2mat(dep_blk);
        F=scatteredInterpolant(dep_blk(:,1),dep_blk(:,2),dep_blk(:,3));
        BO_f=fullfile(DIRBLK,['BO_',num2str(NB1),'_',num2str(NB2),'.txt']);
        Fid=fopen(BO_f,'r');
        if Fid >= 0
          fprintf('READ INTERFACE BOUDARY DEFINITION FILE : %s \n',BO_f)
          bound_blk=textscan(Fid,'%f%f'); fclose(Fid);     
          bound_blk=cell2mat(bound_blk);
        else
          IDB=boundary(dep_blk(:,1),dep_blk(:,2));
          bound_blk=dep_blk(IDB,:);
        end
        iNB=intersect(find(PRM.OptB1==NB1),find(PRM.OptB2==NB2));
        if isempty(iNB)
          int_bo=int_tri;
        else
          int_bo=PRM.OptINT(iNB);
        end
        [p,Bstri]=mesh2D_uni(bound_blk,int_bo,bound_blk);
        Bslon=p(:,1);
        Bslat=p(:,2);
        Bsdep=F(Bslon,Bslat);
       else
        Bstri=[];
        LENG=length(BLK(1).BOUND(NB1,NB2).LON);
        if LENG~=0
          Bslon=[BLK(1).BOUND(NB1,NB2).LON;BLK(1).BOUND(NB1,NB2).LON(1);(BLK(1).BOUND(NB1,NB2).LON(1:LENG-1)+BLK(1).BOUND(NB1,NB2).LON(2:LENG))./2;BLK(1).BOUND(NB1,NB2).LON(LENG)];
          Bslat=[BLK(1).BOUND(NB1,NB2).LAT;BLK(1).BOUND(NB1,NB2).LAT(1);(BLK(1).BOUND(NB1,NB2).LAT(1:LENG-1)+BLK(1).BOUND(NB1,NB2).LAT(2:LENG))./2;BLK(1).BOUND(NB1,NB2).LAT(LENG)];
          Bsdep=[zeros(LENG,1)            ; dep_limit_low.*ones(LENG+1,1)];
          Bstri(1:LENG-1     ,1:3)=[1     :LENG-1;     2:LENG    ;LENG+2:2*LENG]';
          Bstri(LENG:2*LENG-1,1:3)=[LENG+1:2*LENG;LENG+2:2*LENG+1;     1:  LENG]';
          fprintf('BLOCK INTERFACE: %2i  %2i AUTO SET %4i \n',NB1,NB2,(LENG-1)*2+1)
          BLK(1).BOUND(NB1,NB2).type=5;
        else
          BLK(1).BOUND(NB1,NB2).blon=[];
          BLK(1).BOUND(NB1,NB2).blat=[];
          BLK(1).BOUND(NB1,NB2).bdep=[];          
        end
      end
      if ~isempty(Bstri)
        BLK(1).BOUND(NB1,NB2).blon=[Bslon(Bstri(:,1)),Bslon(Bstri(:,2)),Bslon(Bstri(:,3))];
        BLK(1).BOUND(NB1,NB2).blat=[Bslat(Bstri(:,1)),Bslat(Bstri(:,2)),Bslat(Bstri(:,3))];
        BLK(1).BOUND(NB1,NB2).bdep=[Bsdep(Bstri(:,1)),Bsdep(Bstri(:,2)),Bsdep(Bstri(:,3))];
      end
    end
    BLK(1).NB=BLK(1).NB+size(BLK(1).BOUND(NB1,NB2).blon,1);
  end
end
end
%% READ FIX EULER POLES
function [POL,PRM]=READ_EULER_POLES(BLK,PRM)
% Fix euler poles at the block which has no observation site.
% BLID  : Block ID that includes fix POLE
% OMEGA : unit is deg/Myr
POL.ID=false;
if exist(PRM.FilePole,'file')~=2; POL.FIXflag=0; return; end
% 
POL.FIXflag=1;
FID=fopen(PRM.FilePole,'r');
TMP=fscanf(FID,'%d %d %f %f %f\n',[5,Inf]);
POL.ID=zeros(1,BLK(1).NBlock);
FIXw=zeros(BLK(1).NBlock,3);
POL.FLAG =TMP(1,:);
POL.BLID =TMP(2,:);
POL.LAT  =TMP(3,:);
POL.LON  =TMP(4,:);
POL.OMEGA=TMP(5,:);
POL.FLAG =logical(POL.FLAG);   % use or not
POL.BLID =POL.BLID(POL.FLAG) ;
POL.LAT  =POL.LAT(POL.FLAG)  ;
POL.LON  =POL.LON(POL.FLAG)  ;
POL.OMEGA=POL.OMEGA(POL.FLAG);
POL.LAT  =deg2rad(POL.LAT)        ;
POL.LON  =deg2rad(POL.LON)        ;
POL.OMEGA=deg2rad(POL.OMEGA.*1e-6);
POL.wx=POL.OMEGA.*cos(POL.LAT).*cos(POL.LON);
POL.wy=POL.OMEGA.*cos(POL.LAT).*sin(POL.LON);
POL.wz=POL.OMEGA.*sin(POL.LAT)              ;
POL.ID(POL.BLID)=true;
FIXw(POL.BLID,1)=POL.wx;
FIXw(POL.BLID,2)=POL.wy;
FIXw(POL.BLID,3)=POL.wz;
POL.ID=logical(reshape(repmat(POL.ID,3,1),3*BLK(1).NBlock,1));
POL.FIXw=reshape(FIXw',3*BLK(1).NBlock,1);
% 
PRM.APRIORIPOLE=TMP';
% 
end
%% READ RIGID BLOCK BOUNDARY PAIR
function [BLK,PRM]=READ_RIGID_BOUND(BLK,PRM)
if exist(PRM.FileRigb,'file')~=2; return; end
FID=fopen(PRM.FileRigb,'r');
TMP=fscanf(FID,'%d %d %d\n',[3 Inf]);
BLK(1).RGPAIR=TMP';
end
%% MAKE GREEN FUNCTION
function [TRI,OBS]=GREEN_TRI(BLK,OBS)
% Coded by Takeo Ito 2017/01/02 (ver 1.1)
PR=0.25;
ND=size(OBS(1).ALAT,2);
%
ALAT=mean(OBS(1).ALAT(:));
ALON=mean(OBS(1).ALON(:));
OBSMTR=[OBS(1).AXYZ(:,1:3) ones(OBS(1).NOBS,1)];
[OBSx,OBSy]=PLTXY(OBS(1).ALAT,OBS(1).ALON,ALAT,ALON);
OBSz=1e-3.*OBS(1).AHIG;
%
TRI(1).OBSDIS=[];
% TRI(1).AXYZ=[];
% TRI(1).NORMXYZ=[];
% TRI(1).PLANED=[];
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
        TRIEDGE=zeros(3,7);
        for Nn=1:3
          TRIEDGE(Nn,:)=conv2ell(BLK(1).BOUND(NB1,NB2).blat(N,Nn),BLK(1).BOUND(NB1,NB2).blon(N,Nn),BLK(1).BOUND(NB1,NB2).bdep(N,Nn));
        end
        TRI(1).BOUND(NB1,NB2).NORMXYZ(N,:)=cross(TRIEDGE(1,1:3),TRIEDGE(2,1:3));
        TRI(1).BOUND(NB1,NB2).PLANED(N,1)=-TRI(1).BOUND(NB1,NB2).NORMXYZ(N,1)*TRIEDGE(1,1)-TRI(1).BOUND(NB1,NB2).NORMXYZ(N,2)*TRIEDGE(1,2)-TRI(1).BOUND(NB1,NB2).NORMXYZ(N,3)*TRIEDGE(1,3);
        TRI(1).BOUND(NB1,NB2).OBSDIS(:,N)=abs(OBSMTR*[TRI(1).BOUND(NB1,NB2).NORMXYZ(N,:) TRI(1).BOUND(NB1,NB2).PLANED(N,1)]')./sqrt(TRI(1).BOUND(NB1,NB2).NORMXYZ(N,1)^2+TRI(1).BOUND(NB1,NB2).NORMXYZ(N,2)^2+TRI(1).BOUND(NB1,NB2).NORMXYZ(N,3)^2);
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
        TRI(1).BOUND(NB1,NB2).OXYZ(N,:)=conv2ell(TRI(1).BOUND(NB1,NB2).clat(N),TRI(1).BOUND(NB1,NB2).clon(N),1e3.*TRI(1).BOUND(NB1,NB2).cdep(N));
        U=CalcTriDisps(OBSx,OBSy,OBSz,TRIx,TRIy,TRIz,PR,1,0,0);
        TRI(1).BOUND(NB1,NB2).GSTR(1:3:3*ND,N)=U.x; %E
        TRI(1).BOUND(NB1,NB2).GSTR(2:3:3*ND,N)=U.y; %N
        TRI(1).BOUND(NB1,NB2).GSTR(3:3:3*ND,N)=-U.z; %D
        U=CalcTriDisps(OBSx,OBSy,OBSz,TRIx,TRIy,TRIz,PR,0,1,0);
        TRI(1).BOUND(NB1,NB2).GTNS(1:3:3*ND,N)=U.x; %E
        TRI(1).BOUND(NB1,NB2).GTNS(2:3:3*ND,N)=U.y; %N
        TRI(1).BOUND(NB1,NB2).GTNS(3:3:3*ND,N)=-U.z; %D 
        U=CalcTriDisps(OBSx,OBSy,OBSz,TRIx,TRIy,TRIz,PR,0,0,1);
        TRI(1).BOUND(NB1,NB2).GDIP(1:3:3*ND,N)=U.x; %E
        TRI(1).BOUND(NB1,NB2).GDIP(2:3:3*ND,N)=U.y; %N
        TRI(1).BOUND(NB1,NB2).GDIP(3:3:3*ND,N)=-U.z; %D
        if mod(N,ceil(NF/3)) == 1
          fprintf('MAKE GREEN at TRI sub-faults : %4i / %4i \n',N,NF)
        end
      end
%       [BLK,TRI]=DISCRIMINATE_DIRECTION(BLK,TRI,NB1,NB2);
      TRI(1).TNF=TRI(1).TNF+NF;
      TRI(1).OBSDIS=[TRI(1).OBSDIS TRI(1).BOUND(NB1,NB2).OBSDIS(:,N)];
      OBS(1).Gw=(1./min(TRI(1).OBSDIS,[],2))./max(1./min(TRI(1).OBSDIS,[],2));
    end
  end
end
disp('==================')
disp('PASS GREEN_TRI')
disp('==================')
end
%% TODO: DISCRIMINATE BOUNDARY TYPE AND SUBFAULT SURFACE DIRECTION
function [BLK,TRI]=DISCRIMINATE_DIRECTION(BLK,TRI,NB1,NB2)
% Coded by H.Kimura 2017/4/28 (test ver.)
% BLK(1).BOUND(NB1,NB2).type=5; %flag
switch BLK(1).BOUND(NB1,NB2).type
  case 1
    SFID1=inpolygon(TRI(1).BOUND(NB1,NB2).clon,TRI(1).BOUND(NB1,NB2).clat,BLK(NB1).LON,BLK(NB1).LAT);
    SFID2=inpolygon(TRI(1).BOUND(NB1,NB2).clon,TRI(1).BOUND(NB1,NB2).clat,BLK(NB2).LON,BLK(NB2).LAT);
    if sum(SFID1)>sum(SFID2)
      TRI(1).BOUND(NB1,NB2).SDTINV=    ones(size(TRI(1).BOUND(NB1,NB2).clat,2),1);
    elseif sum(SFID1)<sum(SFID2)
      TRI(1).BOUND(NB1,NB2).SDTINV=-1.*ones(size(TRI(1).BOUND(NB1,NB2).clat,2),1);
    end
  case 5
%     ORTHO=TRI(1).BOUND(NB1,NB2).NV(:,3)==0;
%     CTRI =[TRI(1).BOUND(NB1,NB2).clon' TRI(1).BOUND(NB1,NB2).clat' zeros(size(TRI(1).BOUND(NB1,NB2).clat,2),1)];
%     COUT =inpolygon(CTRI(:,1),CTRI(:,2),BLK(NB2).LON,BLK(NB2).LAT)&~inpolygon(CTRI(:,1),CTRI(:,2),BLK(NB1).LON,BLK(NB1).LAT);
%     DPEND=CTRI+1e-3.*TRI(1).BOUND(NB1,NB2).DP;
%     IDOUT=inpolygon(DPEND(:,1),DPEND(:,2),BLK(NB2).LON,BLK(NB2).LAT);
%     TMPID=or(COUT,and(ORTHO,IDOUT));
%     TRI(1).BOUND(NB1,NB2).SDTINV=-1.*TMPID+~TMPID;
    TRI(1).BOUND(NB1,NB2).SDTINV=ones(size(TRI(1).BOUND(NB1,NB2).clat,2),1);
  otherwise
    fprintf('%s\n','No fault.')
end
% 
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
%% CALC MOTION BLOCKS
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
%% CALCLATION AIC AND BLOCK MOTION
function [BLK,OBS]=CALC_AIC(BLK,OBS)
TSig=0; NumB=0;
BLK(1).POLE=[];
for N=1:BLK(1).NBlock
  Sig=0;EVne=[];POLE=[0; 0; 0];
  OBS(N).EEV=zeros(OBS(N).NBLK,1);
  OBS(N).ENV=zeros(OBS(N).NBLK,1);
  if OBS(N).NBLK~=0
    Sig=0;
    EVne=[0 0];
    if OBS(N).NBLK>=1
      NumB=NumB+1;
      OBS(N).GRweight=OBS(1).Gw(OBS(1).ABLK==N);
      OBS(N).GRweight=reshape(repmat(OBS(1).Gw(OBS(1).ABLK==N),1,2)',2*size(OBS(N).GRweight,1),1);
      [POLE,EVne,Sig]=est_pole_w(OBS(N).OXYZ,OBS(N).Vne,OBS(N).Vww.*OBS(N).GRweight);
      TSig=TSig+Sig.*2.*OBS(N).NBLK;
    end
  end
  BLK(N).SIG=Sig;
  BLK(N).POL=POLE;
  BLK(1).POLE=[BLK(1).POLE;BLK(N).POL];
  OBS(N).EEV=EVne(1:2:end);
  OBS(N).ENV=EVne(2:2:end);
  fprintf('BLOCK=%2d NUM_OBS=%2d Sigma^2=%5.2f ',N,OBS(N).NBLK,Sig)
  [latp,lonp,ang]=xyzp2lla(POLE(1),POLE(2),POLE(3));
  fprintf('Lat:%7.2f deg. Lon:%8.2f deg. Ang:%9.2e deg./m.y. \n',latp,lonp,ang);    
%   if OBS(N).NBLK>=2 
%     fprintf('OBS(E,N) ')
%     fprintf('%5.2f ',OBS(N).Vne);fprintf('\n')
%     fprintf('EST(E,N) ')
%     fprintf('%5.2f ',EVne)      ;fprintf('\n')
%   fprintf('\n')
%   end
end
AIC=(OBS(1).NOBS.*2).*log(TSig./(OBS(1).NOBS.*2))+2.*NumB.*3;
cAIC=AIC+2.*NumB.*3.*(NumB.*3+1)./(OBS(1).NOBS.*2-NumB.*3-1);
fprintf('Sigma^2=%8.3f AIC=%7.3f cAIC=%7.3f K=%2d\n',TSig./(OBS(1).NOBS.*2),AIC,cAIC,NumB.*3)
%
end
%% READ BLOCK BOUNDARY DATA
function [BLK,OBS]=READ_BLOCK_BOUND(DIR,OBS)
EXT='*.txt';
file=dir([DIR,'/',EXT]);
[NBlock,~]=size(file);
BLK(1).NBlock=NBlock;
for NB=1:BLK(1).NBlock
  tmp=load(fullfile(DIR,file(NB).name));
  BLK(NB).name=file(NB).name;
  BLK(NB).LON=tmp(:,1);
  BLK(NB).LAT=tmp(:,2);
end
fprintf('READ BLOCK FILES : %4i \n',BLK(1).NBlock)
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    BLK(1).BOUND(NB1,NB2).LAT=[];
    BLK(1).BOUND(NB1,NB2).LON=[];
    LCa=inpolygon(BLK(NB1).LON,BLK(NB1).LAT,BLK(NB2).LON,BLK(NB2).LAT);
    Ca=find(LCa);
    if ~isempty(Ca) && sum(LCa)~=1
      if and(LCa(1),LCa(end))
        Ca0=find(LCa~=true,1,'last')+1:length(LCa)-1;
        Ca1=1:find(LCa~=true,1,'first')-1;
        Ca=[Ca0 Ca1];
      end
      BLK(1).BOUND(NB1,NB2).LAT=BLK(NB1).LAT(Ca);
      BLK(1).BOUND(NB1,NB2).LON=BLK(NB1).LON(Ca);
      BLK(1).BOUND(NB1,NB2).BXYZ=conv2ell(BLK(1).BOUND(NB1,NB2).LAT,BLK(1).BOUND(NB1,NB2).LON,zeros(size(BLK(1).BOUND(NB1,NB2).LON)));
    end
  end
end
%
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
  OBS(N).OXYZ=conv2ell(OBS(N).LAT,OBS(N).LON,OBS(N).HIG);
  OBS(N).Vne=reshape([OBS(1).EVEC(IND); OBS(1).NVEC(IND)],OBS(N).NBLK.*2,1);
  OBS(N).Vww=reshape([OBS(1).EERR(IND); OBS(1).NERR(IND)],OBS(N).NBLK.*2,1);
end
end
%% READ OBSERVATION DATA
function OBS=READ_OBS(FileOBS)
%-------------------
% INPUT format Observations:
% unit is mm/yr
% site_name lon lat EW_comp. NS_comp. UD_comp. ERR_EW ERR_NS ERR_UD 
%-------------------
Fid_OBS=fopen(FileOBS,'r');
N=0;
while 1
  tline=fgetl(Fid_OBS);
  if ~ischar(tline); break; end
  str=strsplit(tline);
  N=N+1;
  OBS(1).NAME(N) =cellstr(str(1));
  OBS(1).ALON(N) =str2double(cellstr(str(2))); %LON
  OBS(1).ALAT(N) =str2double(cellstr(str(3))); %LAT
  OBS(1).AHIG(N) =str2double(cellstr(str(4))); %HIG
  OBS(1).EVEC(N) =str2double(cellstr(str(5))); %E-W
  OBS(1).NVEC(N) =str2double(cellstr(str(6))); %N-S
  OBS(1).HVEC(N) =str2double(cellstr(str(7))); %U-D
  OBS(1).EERR(N) =str2double(cellstr(str(8))); %E-W
  OBS(1).NERR(N) =str2double(cellstr(str(9))); %N-S
  OBS(1).HERR(N) =str2double(cellstr(str(10))); %U-D
  OBS(1).AXYZ(N,:)=conv2ell(OBS(1).ALAT(N),OBS(1).ALON(N),OBS(1).AHIG(N));
end
OBS(1).NOBS=N;
OBS(1).ABLK=zeros(OBS(1).NOBS,1);
OBS(1).LATMAX=max(OBS(1).ALAT);
OBS(1).LATMIN=min(OBS(1).ALAT);
OBS(1).LONMAX=max(OBS(1).ALON);
OBS(1).LONMIN=min(OBS(1).ALON);
fprintf('==================\n') 
fprintf('Number of GNSS site %4d \n',N)
fprintf('==================\n') 
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
[PL,~,Sigma]=lscov(R,Vne,1./w);
EVne=R*PL;
end
%% PLATE MOTION DUE TO EULER POLE (XYZ)
function Vneu=pole2velo(Pxyz,Oxyz)
% pole2velo Convert velocity from Euler pole. Vectorized.
[Nobs,~]=size(Oxyz);
[Npol,~]=size(Pxyz);
Vxyz=zeros(Npol,3,'double'); 
Vneu=zeros(2.*Nobs,Npol,'double');
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
%% CONVERT TO XYZ FROM ELL AT SURFACE
function [OOxyz]=conv2ell(Olat,Olon,Ohig)
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
%====================================================
function free=waitGPU(varargin)
a=true;
d=gpuDevice;
if isempty(varargin)
    limit=30;
else
    limit=varargin{1};
end
% tic
while a
    if limit<=100
        free=d.FreeMemory/d.TotalMemory*100;
        if free>limit
            break
        end
        pause(0.5);
    elseif limit>100
        free=d.FreeMemory;
        if free>limit
            break
        end
        pause(0.5)
    end
end
% waittime=toc;
% if waittime>0.5
%     disp(['waiting time' num2str(waittime)])
% end
end
%====================================================