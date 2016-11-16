%% MCMC Inversion for Geodetic 
function MCMC_SlipBlock(INPUT_FILE)
% Coded    by Takeo Ito 2011/11/08 (ver 1.0)
% Modified by Takeo Ito 2012/10/26 (ver 1.1)
% Modified by Takeo Ito 2015/11/11 (ver 1.2)
%% READ PARAMETERS
PRM=READ_PARAMETERS(INPUT_FILE);
%% READ OBSERVATIONS
OBS=READ_OBS(PRM);
%% READ FAULTS AND MAKE GREEN FUNCTION FOR TRI
TRI=READ_FAULTS_GREEN_TRI(PRM,OBS);
%% READ FAULTS AND MAKE GREEN FUNCTION FOR RECT
RECT=READ_FAULTS_GREEN_RECT(PRM,OBS);
%save GG.mat
%load GG.mat
%% MCMC COMPUTATION 
[MCRE]=GPU_CALC_MCMC_XYr(PRM,OBS,TRI,RECT,1);
%% MCMC OUT_PUT
for N=1:size(MCRE.SLIP,1);
  VAL(N)=mean(MCRE.SLIP(N,:));
end
OUTPUT_GMT(TRI(1) ,VAL(1:length(TRI(1).LAT))    ,'test_out_GMT_tri.txt',3);
OUTPUT_GMT(RECT(1),VAL(length(TRI(1).LAT)+1:end),'test_out_GMT_rect.txt',4);
end
%% MCMC CALC. FOR GPU
function [MCRE]=GPU_CALC_MCMC_XYr(PRM,OBS,TRI,RECT,devGPU)
RR=(OBS(1).YYY./OBS(1).EEE)*(OBS(1).YYY./OBS(1).EEE)';
fprintf('Residual=%9.3f \n',RR);
% TODO: CHECK GPU etc.
% GPU Initialize 
%
%g=gpuDevice(devGPU);
% TODO: CHECK GPU MEMORY
%g_men=g.TotalMemory; %byte
%reset(g);
%
RWD=PRM.RWD;
%
% TODO: NB FROM READ OBS SECTION
% BLOCK1:GPS DATA REF AM
MC(1).NOBS=0;
NB=2;
for B=1:NB
  MC(B).OLAT=(single(OBS(B).LAT)); %gpuArray
  MC(B).OLON=(single(OBS(B).LON)); %gpuArray
  MC(B).OXYZ=conv2ell(MC(B).OLAT,MC(B).OLON); %
  MC(B).NO=length(OBS(B).LON);
  MC(1).NOBS=MC(1).NOBS+3*MC(B).NO;
end
% 
GEE=(single(repmat(1./OBS(1).EEE',1,PRM.NPL)));%gpuArray
GYY=(single(repmat(   OBS(1).YYY',1,PRM.NPL)));%gpuArray
%
% TODO: NFT FROM READ FUALTS SECTION OF TRI
MC(1).NFLT=0;
F=0;
NFT=1;
MC(1).SOURCE=1;
for FT=1:NFT
  F=F+1;
  MC(F).FLAT=(single(TRI(FT).CLAT));%gpuArray
  MC(F).FLON=(single(TRI(FT).CLON));%gpuArray
  MC(F).NP=length(MC(F).FLAT);
  MC(F).FXYZ=conv2ell(MC(F).FLAT,MC(F).FLON);
  MC(F).STCO=(single(repmat(cosd(TRI(FT).STR)',1,PRM.NPL)));%gpuArray E %CHECK
  MC(F).STSI=(single(repmat(sind(TRI(FT).STR)',1,PRM.NPL)));%gpuArray N %CHECK
  MC(F).GSTR=(single(TRI(FT).GSTR./repmat(OBS(1).EEE',1,MC(F).NP)));%gpuArray
  MC(F).GDIP=(single(TRI(FT).GDIP./repmat(OBS(1).EEE',1,MC(F).NP)));%gpuArray
  if MC(F).SOURCE==1; MC(F).GTNS=(single(TRI(FT).GTNS./repmat(OBS(1).EEE',1,MC(F).NP)));end;
  MC(1).NFLT=MC(1).NFLT+length(MC(F).FLAT);
end
% TODO: NFR FROM READ FAULTS SECTION OF RECT
NFR=1;
for FR=1:NFR
  F=F+1;
  MC(F).FLAT=(single(RECT(FR).CLAT));%gpuArray
  MC(F).FLON=(single(RECT(FR).CLON));%gpuArray
  MC(F).NP=length(MC(F).FLAT);
  MC(F).FXYZ=conv2ell(MC(F).FLAT,MC(F).FLON);
  MC(F).STCO=(single(repmat(cosd(RECT(FR).STR)',1,PRM.NPL)));%gpuArray
  MC(F).STSI=(single(repmat(sind(RECT(FR).STR)',1,PRM.NPL)));%gpuArray
  MC(F).GSTR=(single(RECT(FR).GSTR./repmat(OBS(1).EEE',1,MC(F).NP)));%gpuArray
  MC(F).GDIP=(single(RECT(FR).GDIP./repmat(OBS(1).EEE',1,MC(F).NP)));%gpuArray
  MC(1).NFLT=MC(1).NFLT+length(MC(F).FLAT);
end
% TODO: POLE
NP=2;
%MC(1).POLE=[ 63.17 -122.82 0.297]; % AM  63.17 -122.82 0.297 NNR-MORVEL56
%MC(2).POLE=[-46.02  -31.36 0.910]; % PS -46.02  -31.36 0.910 NNR-MORVEL56
%for P=1:NP
%  [Px,Py,Pz]=ell2xyz(MC(P).POLE(1),MC(P).POLE(2),0);
%  MC(P).PXYZ=(1e-6.*pi./180).*MC(P).POLE(3).*[Px,Py,Pz]./norm([Px,Py,Pz]);
%  MC(P).PXYZ=single(repmat(MC(P).PXYZ',1,PRM.NPL));
%end
%MC(1).BVEL=pole2velo(MC(P).PXYZ',MC(1).OXYZ);
%MC(1).BVEL=zeros(size(MC(1).BVEL));
%
LDIM=PRM.NPL.*PRM.KEP;
% TODO: DETERMINED NP (NUMBER OF ESTIMATED POLE)
NEP=3;
SLIP.INT=1e-2;
POLE.INT=1e-10;
LAMD.INT=1e+1;
SLIP.STD=SLIP.INT.*ones(MC(1).NFLT,PRM.NPL,'single');%gpuArray
POLE.STD=POLE.INT.*ones(     3*NEP,PRM.NPL,'single');%gpuArray
LAMD.STD=LAMD.INT.*ones(         1,PRM.NPL,'single');%gpuArray
SLIP.OLD=0.5.*ones(MC(1).NFLT,PRM.NPL,'single');%gpuArray
POLE.OLD=zeros(     3*NEP,PRM.NPL,'single');%gpuArray
LAMD.OLD=zeros(         1,PRM.NPL,'single');%gpuArray
SLIP.CHA=zeros(MC(1).NFLT,LDIM,'single');%gpuArray
POLE.CHA=zeros(     3*NEP,LDIM,'single');%gpuArray
LAMD.CHA=zeros(         1,LDIM,'single');%gpuArray
%
if MC(1).SOURCE==1; 
  TENS.INT=1e-2;
  TENS.STD=TENS.INT.*ones(MC(1).NFLT,PRM.NPL,'single');%gpuArray
  TENS.OLD=zeros(MC(1).NFLT,PRM.NPL,'single');%gpuArray
  TENS.CHA=zeros(MC(1).NFLT,LDIM,'single');%gpuArray
end
%
RES.OLD=inf(1,PRM.NPL,'single');%gpuArray
PRI.OLD=inf(1,PRM.NPL,'single');%gpuArray
%
RT=0;
COUNT=0;
fprintf('USE GPU Max Chain=%4d PP=%5d Nitr=%2d M=%3d \n',...
           PRM.CHA,PRM.NPL,PRM.ITR,MC(1).NFLT);
%
LO_LIMIT=0;
UP_LIMIT=1;
while not(COUNT==2)
  RT=RT+1;
  NACC=0;tic
  for iT=1:PRM.CHA
% SAMPLE SECTION
    SLIP.SMP=SLIP.OLD+RWD.*SLIP.STD.*(rand(MC(1).NFLT,PRM.NPL,'single')-0.5);
    POLE.SMP=POLE.OLD+RWD.*POLE.STD.*(rand(     3*NEP,PRM.NPL,'single')-0.5);
    LAMD.SMP=LAMD.OLD+RWD.*LAMD.STD.*(rand(         1,PRM.NPL,'single')-0.5);
    if MC(1).SOURCE==1;
      TENS.SMP=TENS.OLD+RWD.*TENS.STD.*(rand(MC(1).NFLT,PRM.NPL,'single')-0.5);
    end
% RESAMPLE SECTION
    IND_S=find(SLIP.SMP<LO_LIMIT | SLIP.SMP>UP_LIMIT);
    while isempty(IND_S)==0
      SLIP.SMP(IND_S)=SLIP.OLD(IND_S)+...
               RWD.*SLIP.STD(IND_S).*(rand(length(IND_S),1,'single')-0.5);
      IND_S=find(SLIP.SMP<LO_LIMIT | SLIP.SMP>UP_LIMIT);
      if isempty(IND_S)==1; break; end
    end
% CORRECTION FOR PDF DUE TO RESAMPLE EFFECT
%   WD=RWD.*SLIP.STD;
%   Q_CORR=(min(SLIP.OLD-LO_LIMIT,WD)+min(UP_LIMIT-SLIP.OLD,WD))./...
%          (min(SLIP.SMP-LO_LIMIT,WD)+min(UP_LIMIT-SLIP.SMP,WD));
% TODO:MORE SYSTEMATIC
% BLOCK MOTION SECTION
    MC(1).BVEL=pole2velo((POLE.SMP(1:3,:))',MC(1).OXYZ);
    MC(2).BVEL=pole2velo((POLE.SMP(4:6,:))',MC(2).OXYZ);
% FAULT MOTION SECTION
    MC(1).PVEL=pole2velo((POLE.SMP(7:9,:)-POLE.SMP(4:6,:))',MC(1).FXYZ); %TRENCH
    MC(2).PVEL=pole2velo((POLE.SMP(4:6,:)-POLE.SMP(1:3,:))',MC(2).FXYZ); %MTL
% TODO:NEED MORE SYSTEMATIC AND EFFICIENT
% CALC APRIORI AND RESIDUAL COUPLING RATE SECTION 
    EN=0;
    CAL=zeros(MC(1).NOBS,PRM.NPL,'single');%gpuArray
    for F=1:NP
      SN=EN+1;
      EN=SN+MC(F).NP-1;
      MC(F).ASTR=SLIP.SMP(SN:EN,:).*( MC(F).STCO.*MC(F).PVEL(1:2:end,:)+MC(F).STSI.*MC(F).PVEL(2:2:end,:));
      MC(F).ADIP=SLIP.SMP(SN:EN,:).*(-MC(F).STSI.*MC(F).PVEL(1:2:end,:)+MC(F).STCO.*MC(F).PVEL(2:2:end,:));
      CAL=CAL-(MC(F).GSTR*MC(F).ASTR+MC(F).GDIP*MC(F).ADIP);
      if MC(F).SOURCE==1;
        CAL=CAL-MC(F).GTNS*TENS.SMP(SN:EN,:);
      end;
    end
    PRI.SMP=sum(MC(2).ADIP.^2);
% CALC RESIDUAL SECTION BLOCK MOTION
% TODO: NEED TEST FOR ESTIMATING CAL
    SN=1;
    for B=1:NB
      EN=SN+3*MC(B).NO-1;
      CAL(SN  :3:EN,:)=CAL(SN  :3:EN,:)+MC(B).BVEL(1:2:end,:); %EAST
      CAL(SN+1:3:EN,:)=CAL(SN+1:3:EN,:)+MC(B).BVEL(2:2:end,:); %NORTH
      SN=EN+1;
    end
    RES.SMP=sum((GEE.*(GYY-CAL)).^2,1);
%% MAKE Probably Density Function
% $$ PDF_{post}=\frac{\frac{1}{\sqrt{2\pi\exp(L)}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}}{2}}\exp{\frac{-M^{2}}{2\times\exp{L}}}{\frac{1}{\sqrt{2\pi\exp(L_{old})}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}_{old}}{2}}\exp{\frac{-M^{2}_{old}}{2\times\exp{L_{old}}}} $$
%%
%    Pdf=expm1(0.5.*(-RES.SMP+RES.OLD))+1;
    Pdf=expm1(-0.5.*...
            ((RES.SMP+LAMD.SMP+exp(-LAMD.SMP).*PRI.SMP)...
            -(RES.OLD+LAMD.OLD+exp(-LAMD.OLD).*PRI.OLD)))+1;
% TODO:‚¤[‚ñ‚â‚Á‚Ï‚èƒ_ƒB
%    IND_M=(Pdf.*Q_CORR)>rand(1,PRM.NPL,'single');
    IND_M=Pdf>rand(1,PRM.NPL,'single');
% REVISE SECTION
    SLIP.OLD(:,IND_M) = SLIP.SMP(:,IND_M);
    POLE.OLD(:,IND_M) = POLE.SMP(:,IND_M);
    LAMD.OLD(:,IND_M) = LAMD.SMP(:,IND_M);
    if MC(1).SOURCE==1; TENS.OLD(:,IND_M) = TENS.SMP(:,IND_M); end;
    RES.OLD(IND_M)    = RES.SMP(IND_M);
    PRI.OLD(IND_M)    = PRI.SMP(IND_M);
% KEEP SECTION
    if iT > PRM.CHA-PRM.KEP
      SN=(iT-(PRM.CHA-PRM.KEP)-1)*PRM.NPL+1;
      EN=(iT-(PRM.CHA-PRM.KEP))  *PRM.NPL;
      SLIP.CHA(:,SN:EN)=SLIP.SMP;
      POLE.CHA(:,SN:EN)=POLE.SMP;
      LAMD.CHA(:,SN:EN)=LAMD.SMP;
      if MC(1).SOURCE==1; TENS.CHA(:,SN:EN) = TENS.SMP; end;
      NACC=NACC+sum(IND_M);
    end
  end
%
  AJR=NACC./LDIM;
%
  SLIP.STD=repmat(std(SLIP.CHA,1,2),1,PRM.NPL);
  POLE.STD=repmat(std(POLE.CHA,1,2),1,PRM.NPL);
  LAMD.STD=repmat(std(LAMD.CHA,1,2),1,PRM.NPL);
  if MC(1).SOURCE==1; TENS.STD=repmat(std(TENS.CHA,1,2),1,PRM.NPL); end;
%
%  SLIP.OLD=repmat(mean(SLIP.CHA,2),1,PRM.NPL);
%  POLE.OLD=repmat(mean(POLE.CHA,2),1,PRM.NPL);
%  LAMD.OLD=repmat(mean(LAMD.CHA,2),1,PRM.NPL);
%
  fprintf('T=%3d MaxRes=%6.3f MinRes=%6.3f Accept=%5.1f RWD=%5.2f Time=%5.1fsec\n',...
           RT,1-max(RES.OLD)./RR,1-min(RES.OLD)./RR,100*AJR,RWD,toc)
%TODO: NUMBER OF ESTIMATED POLE
  fprintf('POLE=%9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n',mean(POLE.OLD,2));
  fprintf('VEL_T(E,N)=%9.2e %9.2e VEL_M(E,N)=%9.2e %9.2e \n',...
  mean(mean(MC(1).PVEL(1:2:end,:))),mean(mean(MC(1).PVEL(2:2:end,:))),...
  mean(mean(MC(2).PVEL(1:2:end,:))),mean(mean(MC(2).PVEL(2:2:end,:))));
%  for N=1:5
%    fprintf('SLIP=%4.2f ',mean(SLIP.OLD(N,:)))
%  end
  if AJR > 0.24
    RWD=RWD*1.1;
    COUNT=0;
  elseif AJR < 0.22
    RWD=RWD*0.618;
    COUNT=0;
  else
    COUNT=COUNT+1;
  end
  if RT > PRM.ITR; break; end;
end
%% FIGURES
figure(100);clf(100)
for NPL=1:PRM.NPL
  quiver(OBS(1).ALON,OBS(1).ALAT,CAL(1:3:end,NPL)',CAL(2:3:end,NPL)','blue')
  hold on
end
quiver(OBS(1).ALON,OBS(1).ALAT,GYY(1:3:end,1)',GYY(2:3:end,1)','red')
hold on
patch( TRI(1).LON', TRI(1).LAT', TRI(1).HIG',repmat(mean(SLIP.CHA(1:MC(1).NP,:),2)'    ,3,1))
hold on
patch(RECT(1).LON',RECT(1).LAT',RECT(1).HIG',repmat(mean(SLIP.CHA(MC(1).NP+1:end,:),2)',4,1))
%
figure(105);clf(105)
for NPL=1:PRM.NPL
  quiver(OBS(1).ALON,OBS(1).ALAT,CAL(1:3:end,NPL)',CAL(2:3:end,NPL)','blue')
  hold on
end
quiver(OBS(1).ALON,OBS(1).ALAT,GYY(1:3:end,1)',GYY(2:3:end,1)','red')
hold on
patch( TRI(1).LON', TRI(1).LAT', TRI(1).HIG',repmat(mean(TENS.CHA(1:MC(1).NP,:),2)'    ,3,1))
hold on
patch(RECT(1).LON',RECT(1).LAT',RECT(1).HIG',repmat(mean(TENS.CHA(MC(1).NP+1:end,:),2)',4,1))
%
figure(110);clf(110)
patch( TRI(1).LON', TRI(1).LAT', TRI(1).HIG',repmat(std(SLIP.CHA(1:MC(1).NP,:),0,2)'    ,3,1))
hold on
patch(RECT(1).LON',RECT(1).LAT',RECT(1).HIG',repmat(std(SLIP.CHA(MC(1).NP+1:end,:),0,2)',4,1))
%
figure(120);clf(120)
PVEL_T=pole2velo((POLE.CHA(7:9,:)-POLE.CHA(4:6,:))',MC(1).FXYZ); %TRENCH
PVEL_M=pole2velo((POLE.CHA(4:6,:)-POLE.CHA(1:3,:))',MC(2).FXYZ); %MTL
patch( TRI(1).LON', TRI(1).LAT', TRI(1).HIG',repmat(mean(sqrt(PVEL_T(1:2:end,:).^2+PVEL_T(2:2:end,:).^2),2)',3,1))
hold on
patch(RECT(1).LON',RECT(1).LAT',RECT(1).HIG',repmat(mean(sqrt(PVEL_M(1:2:end,:).^2+PVEL_M(2:2:end,:).^2),2)',4,1))
hold on
quiver3(double( TRI(1).CLON),double( TRI(1).CLAT),double( TRI(1).CHIG),double(mean(PVEL_T(1:2:end,:),2)'),double(mean(PVEL_T(2:2:end,:),2)'),zeros(size( TRI(1).CHIG)))
hold on
quiver3(double(RECT(1).CLON),double(RECT(1).CLAT),double(RECT(1).CHIG),double(mean(PVEL_M(1:2:end,:),2)'),double(mean(PVEL_M(2:2:end,:),2)'),zeros(size(RECT(1).CHIG)))
%
figure(130);clf(130)
BVEL_1=pole2velo((POLE.CHA(1:3,:))',MC(1).OXYZ);
BVEL_2=pole2velo((POLE.CHA(4:6,:))',MC(2).OXYZ);
for NPL=1:PRM.NPL
  quiver(MC(1).OLON,MC(1).OLAT,BVEL_1(1:2:end,NPL)',BVEL_1(2:2:end,NPL)','blue')
  hold on
  quiver(MC(2).OLON,MC(2).OLAT,BVEL_2(1:2:end,NPL)',BVEL_2(2:2:end,NPL)','red')
  hold on
end
%
figure(140);clf(140)
hist(LAMD.CHA(:),100)
%
%MCRE.SLIP=gather(SLIP.CHA);
%MCRE.POLE=gather(POLE.CHA);
%MCRE.LAMD=gather(LAMD.CHA);
%
MCRE.SLIP=SLIP.CHA;
MCRE.POLE=POLE.CHA;
MCRE.LAMD=LAMD.CHA;
%
%KEEP.SLIPMAX=max(MCRE.SLIP);
%KEEP.SLIPMIN=min(MCRE.SLIP);
%KEEP.SLIPINT=(KEEP.SLIPMAX-KEEP.SLIPMIN)./254;
%KEEP.POLEMAX=max(MCRE.POLE);
%KEEP.POLEMIN=min(MCRE.POLE);
%KEEP.POLEINT=(KEEP.POLEMAX-KEEP.POLEMIN)./254;
%KEEP.LAMDMAX=max(MCRE.LAMD);
%KEEP.LAMDMIN=min(MCRE.LAMD);
%KEEP.LAMDINT=(KEEP.LAMDMAX-KEEP.LAMDMIN)./254;
%TODO: NOT CLEAN
%KEEP.SLIP=uint8(round((MCRE.SLIP-repmat(LDIM,KEEP.SLIPMIN,1))./repmat(KEEP.SLIPINT,LDIM,1))); %CHECK
%KEEP.POLE=uint8(round((MCRE.POLE-repmat(LDIM,KEEP.POLEMIN,1))./repmat(KEEP.POLEINT,LDIM,1))); %CHECK
%KEEP.LAMD=uint8(round((MCRE.LAMD-repmat(LDIM,KEEP.LAMDMIN,1))./repmat(KEEP.LAMDINT,LDIM,1))); %CHECK
fprintf('FINISH GPU_CALC_MCMC_XYr \n==================\n')
end
%% READ TRI FAULT AND MAKE GREEN FUNCTION
function [TRI]=READ_FAULTS_GREEN_TRI(PRM,OBS)
% Coded by Takeo Ito 2015/11/11 (ver 1.0)
PR=0.25;
NF=0;
% TODO: READ MULTIPLE BLOCK VELOCITY
load PRM.TRI_F
%
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
%% READ RECT FAULT AND MAKE GREEN FUNCTION
function [RECT]=READ_FAULTS_GREEN_RECT(PRM,OBS)
% Coded    by Takeo Ito 2015/11/12 (ver 1.0)
PR=0.25;
Fid=fopen(PRM.RECT_F,'r');
NF=0;
while 1
  NF=NF+1;
  loc_f=fscanf(Fid,'%f %f %f \n', [3 4]);
  [~] = fgetl(Fid);
  RECT(1).LON(NF,:)=loc_f(1,:);%LON
  RECT(1).LAT(NF,:)=loc_f(2,:);%LAT
  RECT(1).HIG(NF,:)=loc_f(3,:);%HIG
  tline = fgetl(Fid); if ~ischar(tline); break; end
end
fclose(Fid);
%
ND=size(OBS(1).ALAT,2);
RECT(1).GSTR=zeros(3*ND,NF);
RECT(1).GDIP=zeros(3*ND,NF);
%
fprintf('==================\nNumber of RECT sub-faults : %i \n',NF)
%
for N=1:NF
  F_LOC=[RECT(1).LAT(N,:);RECT(1).LON(N,:);RECT(1).HIG(N,:)]';
  [F,DA,STR,DIP,NV,ST,DP,AL,AW]=EST_FAULT_RECT(F_LOC);
  RECT(1).CLAT(N)=F(1);
  RECT(1).CLON(N)=F(2);
  RECT(1).CHIG(N)=F(3);
  RECT(1).DA(N)=DA;
  RECT(1).STR(N)=STR;
  RECT(1).DIP(N)=DIP;
  RECT(1).NV(N,:)=NV;
  RECT(1).ST(N,:)=ST;
  RECT(1).DP(N,:)=DP;
  U=DISLOC_RECT(OBS(1).ALAT,OBS(1).ALON,F(1),F(2),-F(3),STR+90,DIP,0,AL,AW,PR);
  RECT(1).GSTR(1:3:3*ND,N)=U.E; %E
  RECT(1).GSTR(2:3:3*ND,N)=U.N; %N
  RECT(1).GSTR(3:3:3*ND,N)=U.U; %D
  U=DISLOC_RECT(OBS(1).ALAT,OBS(1).ALON,F(1),F(2),-F(3),STR+90,DIP,90,AL,AW,PR);
  RECT(1).GDIP(1:3:3*ND,N)=U.E; %E
  RECT(1).GDIP(2:3:3*ND,N)=U.N; %N
  RECT(1).GDIP(3:3:3*ND,N)=U.U; %D
end
fprintf('Size of Green Matrix : %i x %i \n',size(RECT(1).GSTR,1),size(RECT(1).GSTR,2))
fprintf('PASS READ_FAULTS_GREEN_RECT \n==================\n')
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
%% ESTIMATE FAULT PARAMETERS FOR RECT
function [FLOC,DA,STR,DIP,NV,ST,DP,AL,AW]=EST_FAULT_RECT(loc_f)
% Coded by Takeo Ito 2015/11/11 (ver 1.0)
[X,Y]=PLTXY(loc_f(:,2),loc_f(:,1),loc_f(1,2),loc_f(1,1));
[DA]=AREA_RECT(X,Y,loc_f(:,3));
[AL,AW]=F_LENG_RECT(X,Y,loc_f(:,3),DA);
FLOC=mean(loc_f);
[STR,DIP,NV,ST,DP]=EST_STRDIP_RECT(X,Y,loc_f(:,3));
end
%% ESTIMATE AREA AT SUB-FAULT FOR RECT
function [DA]=AREA_RECT(X,Y,Z)
% CALC. AREA IN THREE DIMENSION USING HERON'S FOMULA
% Coded by Takeo Ito 2006/03/04 (ver 1.0)
LENG(1)=sqrt((X(1)-X(2)).^2+(Y(1)-Y(2)).^2+(Z(1)-Z(2)).^2);
LENG(2)=sqrt((X(3)-X(2)).^2+(Y(3)-Y(2)).^2+(Z(3)-Z(2)).^2);
LENG(3)=sqrt((X(3)-X(4)).^2+(Y(3)-Y(4)).^2+(Z(3)-Z(4)).^2);
LENG(4)=sqrt((X(1)-X(4)).^2+(Y(1)-Y(4)).^2+(Z(1)-Z(4)).^2);
LENG(5)=sqrt((X(2)-X(4)).^2+(Y(2)-Y(4)).^2+(Z(2)-Z(4)).^2);
S1=(LENG(1)+LENG(2)+LENG(5))./2;
S2=(LENG(3)+LENG(4)+LENG(5))./2;
DA1=sqrt(S1*(S1-LENG(1))*(S1-LENG(2))*(S1-LENG(5)));
DA2=sqrt(S2*(S2-LENG(3))*(S2-LENG(4))*(S2-LENG(5)));
DA=real(DA1+DA2);
end
%% ESTIMATE FAULT LENGTH AND WIDTH FOR RECT
function [AL,AW]=F_LENG_RECT(X,Y,Z,DA)
% Coded by Takeo Ito 2006/03/04 (ver 1.0)
ALt(1)=sqrt((X(1)-X(2)).^2+(Y(1)-Y(2)).^2+(Z(1)-Z(2)).^2);
ALt(2)=sqrt((X(3)-X(4)).^2+(Y(3)-Y(4)).^2+(Z(3)-Z(4)).^2);
AL=(ALt(1)+ALt(2))./2;
AW=DA./AL;
end
%% ESTIMATE STRKE AND DIP FOR RECT FAULT
function [STR,DIP,NV,ST,DP]=EST_STRDIP_RECT(X,Y,Z)
%==========
% CALC. STR AND DIP ON FAULT
% CODE BY T.ITO (2006/03/04)
% Modified by T.ITO (2016/02/16)
% DEPTH IS MINUS
%==========
[STR1,DIP1,NV1,ST1,DP1]=EST_STRDIP_TRI(X(1:3),Y(1:3),Z(1:3));
[STR2,DIP2,NV2,ST2,DP2]=EST_STRDIP_TRI(X(2:4),Y(2:4),Z(2:4));
STR=(STR1+STR2)./2;
DIP=(DIP1+DIP2)./2;
NV=(NV1+NV2)./2;
ST=(ST1+ST2)./2;
DP=(DP1+DP2)./2;
end
%% READ PARAMETER FOR MCMC Inversion 
function [PRM]=READ_PARAMETERS(INPUT_SET)
% MCMC Inversion for Geodetic 
% Coded    by Takeo Ito 2011/11/08 (ver 1.0)
% Modified by Takeo Ito 2012/10/26 (ver 1.1)
% Modified by Takeo Ito 2015/11/11 (ver 1.2)
% Modified by Takeo Ito 2016/07/06 (ver 1.3)
%
Fid=fopen(INPUT_SET,'r');
PRM.HOME_D=pwd;
TRI_F=fscanf(Fid,'%s \n',[1,1]);
PRM.TRI_F=fullfile(PRM.HOME_D,TRI_F);
[~]=fgetl(Fid);
RECT_F=fscanf(Fid,'%s \n',[1,1]);
PRM.RECT_F=fullfile(PRM.HOME_D,RECT_F);
[~]=fgetl(Fid);
OBS_F=fscanf(Fid,'%s \n',[1,1]);
PRM.OBS_F=fullfile(PRM.HOME_D,OBS_F);
[~]=fgetl(Fid);
%
PRM.NPL=fscanf(Fid,'%d \n',[1,1]);
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
fprintf('==================\nINPUT PARAMETERS\n==================\n') 
fprintf('HOME_D       : %s \n',PRM.HOME_D) 
fprintf('TRI_F        : %s \n',PRM.TRI_F) 
fprintf('RECT_F       : %s \n',PRM.RECT_F)
fprintf('OBS_F        : %s \n',PRM.OBS_F) 
fprintf('PP           : %i \n',PRM.NPL) 
fprintf('Max_Nitr     : %i \n',PRM.ITR) 
fprintf('Chain        : %i \n',PRM.CHA) 
fprintf('KEEP         : %i \n'   ,PRM.KEP) 
fprintf('Walk_dis     : %4.2f \n',PRM.RWD) 
fprintf('==================\n') 
%====================================================
disp('PASS READ_PARAMETERS')
end
%% OUT PUT RESULT FOR GMT FORMAT
function OUTPUT_GMT(FLOC,VAL,OUT_F,NUM)
% Coded by Takeo Ito 2015/11/11 (ver 1.0)
% 
Fid=fopen(OUT_F,'w');
NF=length(VAL); 
for N=1:NF
  fprintf(Fid,'> -Z%7.5f \n',VAL(N));
  for ND=1:NUM
    fprintf(Fid,'%7.3f %7.4f \n',FLOC.LON(N,ND),FLOC.LAT(N,ND));
  end
    fprintf(Fid,'%7.3f %7.4f \n',FLOC.LON(N,1),FLOC.LAT(N,1));
end
fclose(Fid);
fprintf('OUT PUT GMT FORMAT :%s \n',OUT_F)

end
%% READ OBSERVATION DATA
function [OBS]=READ_OBS(PRM)
% moditied 2015/11/11 by T.ITO
% moditied 2016/07/07 by T.ITO
%-------------------
% format:
% site_name lon lat EW_comp. NS_comp. UD_comp. ERR_EW ERR_NS ERR_UD N_BK
%-------------------
Fid_OBS=fopen(PRM.OBS_F,'r');
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
  OBS(1).YYY(3*N-2) =str2double(cellstr(str(5))); %E-W
  OBS(1).YYY(3*N-1) =str2double(cellstr(str(6))); %N-S
  OBS(1).YYY(3*N)   =str2double(cellstr(str(7))); %U-D
  OBS(1).EEE(3*N-2) =str2double(cellstr(str(8))); %E-W
  OBS(1).EEE(3*N-1) =str2double(cellstr(str(9))); %N-S
  OBS(1).EEE(3*N)   =str2double(cellstr(str(10))); %U-D
end
plot(OBS(1).ALON,OBS(1).ALAT,'.')
% TODO: FOR TEST READ TWO BLOCK REGION
MB=zeros(max(OBS(1).BLK),1);
for NOBS=1:N
  MB(OBS(1).BLK(NOBS))=MB(OBS(1).BLK(NOBS))+1;
  OBS(OBS(1).BLK(N)).LAT(MB(OBS(1).BLK(NOBS)))=OBS(1).ALAT(N);
  OBS(OBS(1).BLK(N)).LON(MB(OBS(1).BLK(NOBS)))=OBS(1).ALON(N);
  OBS(OBS(1).BLK(N)).HIG(MB(OBS(1).BLK(NOBS)))=OBS(1).ALAT(N);
end
fprintf('==================\nNumber of observation site : %i \n',N)
end
%% PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
function [X,Y]=PLTXY(ALAT,ALON,ALAT0,ALON0)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN 
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
A=6.378160e3;   % km e3 --> m e6 --> mm e9
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
%% PLATE MOTION DUE TO EULER POLE (ELL)
function [Vneu]=EULER2VEL(olat,olon,plat,plon,pomega)
% This function computes the linear velocity Vneu(vx,vy,vz) 
% due to angular rotation (pomega[DEG/MYR]) 
% at a set of input points P(plat,plon).
% See Cox and Hart (1986), p. 155.
%
% INPUT:
%    plat, plon, pomega = euler pole[deg], omega
%    olat, olon         = calc points
% OUTPUT:
%    Vneu = local velocities (mm/yr)
%
% euler pole : convert (lat lon omega) => vector
%               deg/Myr --> rad/yr
%
% ph-sunda -58.03 -11.47 1.076
%
% NNR-MORVEL56 http://geoscience.wisc.edu/~chuck/MORVEL/motionframe_nnrm56.html
% AM  63.17 -122.82 0.297
% PS -46.02  -31.36 0.910

% JGR 21,20,2191-2194, 1994 NNR-NUVEL-1A 
% Plate  Lat   Lon    Omega(deg/Myr) errors
% na-pa  48.7   -78.2  0.75   1.3  1.2  -61  0.01
% ri-pa  31.0  -102.4  2.45   3.6  0.6   21  0.57 
% co-pa  36.8  -108.6  2.00   1.0  0.6  -33  0.05 
% ri-na  22.8  -109.4  1.80   1.8  0.6  -57  0.58 
% ri-co   6.8   -83.7  0.54  35.8  1.8  -56  0.52 
% co-na  27.9  -120.7  1.36   1.8  0.7  -67  0.05 
% co-nz   4.8  -124.3  0.91   2.9  1.5  -88  0.05 
% nz-pa  55.6   -90.1  1.36   1.8  0.9   -1  0.02 
% nz-an  40.5   -95.9  0.52   4.5  1.9   -9  0.02 
% nz-sa  56.0   -94.0  0.72   3.6  1.5  -10  0.02 
% an-pa  64.3   -84.0  0.87   1.2  1.0   81  0.01 
% pa-au -60.1  -178.3  1.07   1.0  0.9  -58  0.01 
% eu-pa  61.1   -85.8  0.86   1.3  1.1   90  0.02 
% co-ca  24.1  -119.4  1.31   2.5  1.2  -60  0.05 
% nz-ca  56.2  -104.6  0.55   6.5  2.2  -31  0.03 
% Atlantic  Ocean 
% eu-na  62.4   135.8  0.21   4.1  1.3  -11  0.01 
% af-na  78.8    38.3  0.24   3.7  1.0   77  0.01 
% af-eu  21.0   -20.6  0.12   6.0  0.7   -4  0.02 
% na-sa  16.3   -58.1  0.15   5.9  3.7   -9  0.01 
% af-sa  62.5   -39.4  0.31   2.6  0.8  -11  0.01 
% an-sa  86.4   -40.7  0.26   3.0  1.2  -24  0.01 
% na-ca -74.3   -26.1  0.10  24.7  2.6  -52  0.03 
% ca-sa  50.0   -65.3  0.18  14.9  4.3   -2  0.03 
% Indian  Ocean 
% au-an  13.2    38.2  0.65   1.3  1.0  -63  0.01 
% af-an   5.6   -39.2  0.13   4.4  1.3  -42  0.01 
% au-af  12.4    49.8  0.63   1.2  0.9  -39  0.01 
% au-in  -5.6    77.1  0.30   7.4  3.1  -43  0.07 
% in-af  23.6    28.5  0.41   8.8  1.5  -74  0.06 
% ar-af  24.1    24.0  0.40   4.9  1.3  -65  0.05 
% in-eu  24.4    17.7  0.51   8.8  1.8  -79  0.05 
% ar-eu  24.6    13.7  0.50   5.2  1.7  -72  0.05 
% au-eu  15.1    40.5  0.69   2.1  1.1  -45  0.01 
% in-ar   3.0    91.5  0.03  25.2  2.4  -58  0.04
% 
[NN,MM]=size(olat);
olat=olat(:);
olon=olon(:);
%
[px,py,pz]=ell2xyz(plat,plon,0);
pvec=[px, py, pz];
pvec=pomega.*pvec./norm(pvec).*(1e-6.*pi./180);
%
% local observation points : convert (lat lon) = xyz
%                             m --> mm
[ox,oy,oz]=ell2xyz(olat,olon,0);
Oxyz=[ox, oy, oz];
Oxyz = Oxyz * 1e3;
%
% V = (w E) x (r P), where E is the Euler pole unit vector
Vxyz(:,1) = -Oxyz(:,2).*pvec(3) + pvec(2).*Oxyz(:,3);
Vxyz(:,2) = -Oxyz(:,3).*pvec(1) + pvec(3).*Oxyz(:,1);
Vxyz(:,3) = -Oxyz(:,1).*pvec(2) + pvec(1).*Oxyz(:,2);
%
Vneu=zeros(3,NN*MM);
for N=1:NN*MM
  Vneu(:,N)=xyz2neu([olat(N),olon(N)],Vxyz(N,:));
end
Vneu=reshape(Vneu,3,NN,MM);
%
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
%% CONVERT TO NEU FROM XYZ
function dneu=xyz2neu(bl,dxyz)
% TRANSFORMATION FROM (DX,DY,DZ) => (DN,DE,DU)
% CODE BY T.ITO 2006/12/13 ver0.1
% BUG FIX  BY T.ITO 2007/01/13 ver0.2
% MODIFIED BY T.ITO 2015/05/09 ver0.3
%
deg2rad=pi/180;
bl=bl(1:2).*deg2rad;
dneu=zeros(size(dxyz));
R=[-sin(bl(1)).*cos(bl(2)) -sin(bl(1)).*sin(bl(2)) cos(bl(1)) ; ...
   -sin(bl(2))              cos(bl(2))             0.0        ; ...
    cos(bl(1)).*cos(bl(2))  cos(bl(1)).*sin(bl(2)) sin(bl(1))];
%
dneu(:,1)=R(1,1).*dxyz(:,1)+R(1,2).*dxyz(:,2)+R(1,3).*dxyz(:,3);
dneu(:,2)=R(2,1).*dxyz(:,1)+R(2,2).*dxyz(:,2);
dneu(:,3)=R(3,1).*dxyz(:,1)+R(3,2).*dxyz(:,2)+R(3,3).*dxyz(:,3);
end
%% CONVERT TO ELL FROM XYZ
function [lat,lon,h]=xyz2ell(X,Y,Z)
% XYZ2ELL  Converts cartesian coordinates to ellipsoidal. Vectorized.
% GRS80
a=6378137.0;                % m --> mm
b=6356752.314140356*1000;   % m --> mm
e=0.081819191042815791;
e2=0.006694380022900788;
%
lon=atan2(Y,X);
p=sqrt(X.*X+Y.*Y);
r=sqrt(p.*p+Z.*Z);
u=atan(b.*Z.*(1+e.*b./r)./(a.*p));
lat=atan((Z+e.*b.*sin(u).^3)./(p-e2.*a.*cos(u).^3));
v=a./sqrt(1-e2.*sin(lat).^2);
h=p.*cos(lat)+Z.*sin(lat)-a*a./v;
%
deg=180/pi;
lat=lat.*deg;
lon=lon.*deg;
end
%% PLATE MOTION DUE TO EULER POLE (XYZ) GPU
function Vneu=pole2velo(Pxyz,Oxyz)
% pole2velo Convert velocity from Euler pole. Vectorized.
[Nobs,~]=size(Oxyz);
[Npol,~]=size(Pxyz);
Vxyz=zeros(Npol,3,'single'); %gpuArray
Vneu=zeros(2.*Nobs,Npol,'single');%gpuArray
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
%% CONVERT TO XYZ FROM ELL AT SURFACE
function [OOxyz]=conv2ell(Olat,Olon)
Olat=Olat(:);
Olon=Olon(:);
deg2rad=pi/180;
[Oxyz(:,1),Oxyz(:,2),Oxyz(:,3)]=ell2xyz(Olat,Olon,0);
Oxyz = Oxyz*1e3;
OOxyz=[Oxyz sin(Olat*deg2rad) sin(Olon*deg2rad) cos(Olat*deg2rad) cos(Olon*deg2rad)];
end
%% CALC. DISP. DUE TO SLIP AT TRIANGULAR FAULT
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
% Implements algorithms described in the journal article:
% Meade, B. J. Algorithms for calculating displacements, 
% strains, and stresses for triangular dislocation elements
% in a uniform elastic half space
% Computers and Geosciences, submitted, 2006.
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
[n1,n2]                      = size(sz);
sz1                          = sz(:);
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
for iTri = 1:3
   % Calculate strike and dip of current leg
   strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
   segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
   [rx,ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
   dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));
   
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
   ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
   tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
   dsVec                    = cross(ssVec, tsVec);
   lss                      = dot(slipVec, ssVec);
   lts                      = dot(slipVec, tsVec);
   lds                      = dot(slipVec, dsVec);
   if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
      % First angular dislocation
      [sx1, sy1]               = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
      sx1(abs(sx1)<0.0001)=0.0001;% for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012
      sy1(abs(sy1)<0.0001)=0.0001;
      [ux1, uy1, uz1]          = adv(sx1, sy1, sz1-z(iTri), z(iTri), beta, pr, lss, lts, lds);
      % Second angular dislocation
      [sx2, sy2]               = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike); 
      sx2(abs(sx2)<0.0001)=0.0001;% for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012
      sy2(abs(sy2)<0.0001)=0.0001;
      [ux2, uy2, uz2]           = adv(sx2, sy2, sz1-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
      % Rotate vectors to correct for strike
      [uxn, uyn]                = RotateXyVec(ux1-ux2, uy1-uy2, strike);
      uzn                       = uz1-uz2;
      % Add the displacements from current leg
      U.x                       = U.x + reshape(uxn,n1,n2);
      U.y                       = U.y + reshape(uyn,n1,n2);
      U.z                       = U.z + reshape(uzn,n1,n2);
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
end
%% SUBROUTINE FOR CalcTriDisps (LinePlaneIntersect)
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
%% SUBROUTINE FOR CalcTriDisps (swap)
function [a,b] = swap(a, b)
% Swap two values
temp                            = a;
a                               = b;
b                               = temp;
end
%% SUBROUTINE FOR CalcTriDisps (RotateXyVec)
function [xp,yp] = RotateXyVec(x, y, alpha)
% Rotate a vector by an angle alpha
x                             = x(:);
y                             = y(:);
alpha                         = pi/180*alpha;
xp                            = cos(alpha).*x - sin(alpha).*y;
yp                            = sin(alpha).*x + cos(alpha).*y;
end
%% SUBROUTINE FOR CalcTriDisps (adv)
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
z1bar             =  y1.*cosbeta + y3bar.*sinbeta;
z3bar             = -y1.*sinbeta + y3bar.*cosbeta;
R2bar             = y1.*y1 + y2.*y2 + y3bar.*y3bar;
Rbar              = sqrt(R2bar);
F                 = -atan2(y2, y1) + atan2(y2, z1)    + atan2(y2.*R.*sinbeta,    y1.*z1+(y2.*y2).*cosbeta);
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
%
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
%
v1CB1             = v1CB1 ./ (4.*pi.*(1-nu));
v2CB1             = v2CB1 ./ (4.*pi.*(1-nu));
v3CB1             = v3CB1 ./ (4.*pi.*(1-nu));
%
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
%
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
%
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
%
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
%
v1B3              = v1InfB3 + v1CB3;
v2B3              = v2InfB3 + v2CB3;
v3B3              = v3InfB3 + v3CB3;
% Sum the for each slip component
v1                = B1.*v1B1 + B2.*v1B2 + B3.*v1B3;
v2                = B1.*v2B1 + B2.*v2B2 + B3.*v2B3;
v3                = B1.*v3B1 + B2.*v3B2 + B3.*v3B3;
end
%% CALC. SURFACE DISP. DUE TO SLIP AT RECT FAULT
function [U]=DISLOC_RECT(OLAT,OLON,ALAT0,ALON0,DEP,PHAI,DIP,RAK,AL,AW,NU)
% Coded    by Takeo Ito 2011/11/08 (ver 1.0)
% Modified by Takeo Ito 2012/10/26 (ver 1.1)
% Modified by Takeo Ito 2015/11/11 (ver 1.2)
%-------------------
% OKADA(1985)
%	OLAT, OLON  : Coordinates of observation points (Latitude, Longitude)
%	ALAT0,ALON0 : Coordinates of the fault centroid (Latitude, Longitude)
%	DEP         : Depth of the fault top (DEP > 0)
%	PHAI        : Strike-angle from Norh counterclock wise (in degrees)
%                 fault dips to the right side of the trace
%	DIP         : Dip-angle  (in degrees, must be scalar)
%	RAK         : Rake-angle (in degrees, must be scalar)
%	AL          : Fault length in strike direction (LEN > 0)
%	AW          : Fault width in dip direction (WIDTH > 0)
%	NU          : 0.25 isotropic Poisson's ratio
%-------------------
RAD=pi./180;
nsite=length(OLAT);
%
SR=sin(RAK.*RAD);
CR=cos(RAK.*RAD);
%
SD=sin(DIP.*RAD);
CD=cos(DIP.*RAD);
%
ST=sin(PHAI.*RAD);
CT=cos(PHAI.*RAD);
%
CD(abs(CD)<=1.0e-3)=0;
SD(CD==0 && SD > 0)=1;
SD(CD==0 && SD < 0)=-1;
%
[Xobs,Yobs]=PLTXY(OLAT(:),OLON(:),ALAT0,ALON0);
%
%	X,Y : Coordinates of observation points (relative to fault center at top)
%
Xobs=Xobs.*1000;  % Unit km --> m
Yobs=Yobs.*1000;  % Unit km --> m
DEP=DEP*1000;  % Unit km --> m
AL=AL*1000;    % Unit km --> m
AW=AW*1000;    % Unit km --> m
%
% Converts fault coordinates (E,N,DEPTH) relative to centroid
% into Okada's reference system (X,Y,DEP)
% DEP is fault's top edge
EC = Xobs + CT.*CD.*AW/2;
NC = Yobs - ST.*CD.*AW/2;
X = CT.*NC + ST.*EC + AL/2;
Y = ST.*NC - CT.*EC + CD.*AW;
% Variable substitution (independent from xi and eta)
P = Y.*CD + DEP.*SD.*ones(nsite,1);
Q = Y.*SD - DEP.*CD.*ones(nsite,1);
%
% XI,ET,Q : FAULT COORDINATE
%           ORIGIN AT CENTER OF LOWER SIDE OF FAULT
%
Uxs=zeros(nsite,1); Uys=zeros(nsite,1); Uzs=zeros(nsite,1);
Uxd=zeros(nsite,1); Uyd=zeros(nsite,1); Uzd=zeros(nsite,1);
%
[uxs,uys,uzs,uxd,uyd,uzd]=STATIC(NU,X,P,Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
%
[uxs,uys,uzs,uxd,uyd,uzd]=STATIC(NU,X,P-AW,Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
%
[uxs,uys,uzs,uxd,uyd,uzd]=STATIC(NU,X-AL,P,Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
%
[uxs,uys,uzs,uxd,uyd,uzd]=STATIC(NU,X-AL,P-AW,Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
%
U.E=(-CR.*( Uxs(:).*ST-Uys(:).*CT)-SR.*( Uxd(:).*ST-Uyd(:).*CT))./(2*pi); %E
U.N=(-CR.*( Uxs(:).*CT+Uys(:).*ST)-SR.*( Uxd(:).*CT+Uyd(:).*ST))./(2*pi); %N
U.U=(-CR.*( Uzs(:)               )-SR.*( Uzd(:)               ))./(2*pi); %D
%
end
%% SUBROUTINE FOR OKADA (1985)
function [Uxs,Uys,Uzs,Uxd,Uyd,Uzd]=STATIC(nu,xi,eta,q,sd,cd)
% Coded     by Takeo Ito 2011/11/08 (ver 1.0)
% Bug fixed by Takeo Ito 2016/02/14 (ver 2.0)
X = sqrt(xi.^2 + q.^2);
R  =sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta.*cd + q.*sd;
db = eta.*sd - q.*cd;
Rb = R + db;
cR = cd.*Rb;
Re = R + eta;
lR = log(Re);
NU12=1-2*nu;
if cd == 0
  I5 = NU12 * 2./cd .* atan((eta.*(X + q.*cd) + X.*(R + X).*sd)) ./(xi.*(R + X).*cd);
  I5(xi==0) = 0;
  I4 = NU12 * 1./cd * (log(Rb) - sd.*lR);
  I3 = NU12 * ( yb./cR) + sd./cd.*I4 - lR ;
  I1 = NU12 * (-xi./cR) - sd./cd.*I5;
else
  I5 = -NU12   * xi.*sd./Rb;
  I4 = -NU12   * q./Rb;
  I3 =  NU12/2 * (eta./Rb + yb.*q./Rb.^2 - lR);
  I1 = -NU12/2 * xi.*q./Rb.^2;
end
I2   =  NU12   * (-lR) - I3;
%
Uxs = xi.*q./(R.*Re) + I1.*sd;
Uys = yb.*q./(R.*Re) + q.*cd./Re + I2.*sd;
Uzs = db.*q./(R.*Re) + q.*sd./Re + I4.*sd;
Uxd = q./R - I3.*sd.*cd;
Uyd = yb.*q./(R.*(R + xi)) - I1.*sd.*cd;
Uzd = db.*q./(R.*(R + xi)) - I5.*sd.*cd;
%
k = find(q~=0);
ax = atan(xi(k).*eta(k)./(q(k).*R(k)));
Uxs(k) = Uxs(k) +     ax;
Uyd(k) = Uyd(k) + cd.*ax;
Uzd(k) = Uzd(k) + sd.*ax;
end