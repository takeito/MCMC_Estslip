%% Replica exchange Monte Calro
function [CHA]=PT_REMC(D,G,BLK,PRM,OBS,POL)
% Coded by H. Kimura 2018/06/02 (Original is MH_MCMC coded by T. Ito)
% Replica change Monte Calro (Parallel tempering) method
% Refering Kano et al. (2017); Earl & Deem (2005); Swendsen & Wang (1986).
logfile=fullfile(PRM.DirResult,'log.txt');
logFID=fopen(logfile,'a');
RR=(D(1).OBS./D(1).ERR)'*(D(1).OBS./D(1).ERR);
fprintf('Residual=%9.3f \n',RR);
fprintf(logFID,'Residual=%9.3f \n',RR);
% if PRM.GPU==99
%   precision='double';
% else
%   precision='single';
% end
precision='double';
RWD=PRM.RWD;
Mc.INT=1e-2;
Mp.INT=1e-10;
Mi.INT=1e-10;
La.INT=1e+1;
Mc.N=BLK(1).NB;
Mp.N=3.*BLK(1).NBlock;
Mi.N=3.*BLK(1).NBlock;
La.N=1;
Mc.STD=Mc.INT.*ones(Mc.N,1,precision);
Mp.STD=Mp.INT.*ones(Mp.N,1,precision);
Mi.STD=Mi.INT.*ones(Mi.N,1,precision);
La.STD=La.INT.*ones(La.N,1,precision);
% Mc.OLD=       randn(Mc.N,1,precision);                % Normal distribution
% Mc.OLD=(Mc.OLD-min(Mc.OLD))./max(Mc.OLD-min(Mc.OLD)); % Normal distribution
% Mc.OLD=   -0.5+rand(Mc.N,1,precision);                % Uniform distribution(-1 to 1)
Mc.OLD=  rand(Mc.N,1,precision);                % Uniform distribution( 0 to 1)
Mp.OLD= double(BLK(1).POLE);
Mi.OLD= 1e-10.*(-0.5+rand(Mi.N,1,precision));
La.OLD= zeros(La.N,1,precision);
CHA.Mc= zeros(Mc.N,PRM.KEP,precision);
CHA.Mp= zeros(Mp.N,PRM.KEP,precision);
CHA.Mi= zeros(Mi.N,PRM.KEP,precision);
CHA.La= zeros(La.N,PRM.KEP,precision);
% Set FIX POLES if POL.FIXflag=1
% MpScale=Mp.INT.*ones(Mp.N,1,precision);
Mp.OLD(POL.ID)=0; Mp.OLD=Mp.OLD+POL.FIXw;
Mp.STD(POL.ID)=0;
%
Mi.OLD=Mi.OLD.*BLK(1).IDinter;
Mi.STD=Mi.STD.*BLK(1).IDinter;
% 
RES.OLD=inf(1,1,precision);
% PRI.OLD=inf(1,1,precision);
RWDSCALE=1000*RWD/(PRM.CHA);
McScale=RWDSCALE*0.13;
MpScale=RWDSCALE*(1.3E-9).*ones(Mp.N,1,precision).*~POL.ID;
MiScale=RWDSCALE*1e-10;
% McScale=0.05;
% MpScale=3E-10.*ones(Mp.N,1,precision).*~POL.ID;
LO_Mc=0;
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
  fprintf(logFID,'USE GPU Max Chain=%4d Nitr=%2d Mc=%4d Mp=%3d res_Memory(GB)=%6.3f\n',...
           PRM.CHA,PRM.ITR,Mc.N,Mp.N,res_m./(1024.*1024.*1024));      
%   CHA.Mc=gpuArray(CHA.Mc);
%   CHA.Mp=gpuArray(CHA.Mp);
%   CHA.La=gpuArray(CHA.La);
  Mc.STD=gpuArray(Mc.STD);
  Mp.STD=gpuArray(Mp.STD);
  Mi.STD=gpuArray(Mi.STD);
  La.STD=gpuArray(La.STD);
  Mc.OLD=gpuArray(Mc.OLD);
  Mp.OLD=gpuArray(Mp.OLD);
  Mi.OLD=gpuArray(Mi.OLD);
  La.OLD=gpuArray(La.OLD);
  D(1).OBS=gpuArray(D(1).OBS);
  D(1).ERR=gpuArray(D(1).ERR);
  G.TB=gpuArray(G.TB);
  G.C=gpuArray(G.C);
  G.P=gpuArray(G.P);
  G.I=gpuArray(G.I);
  D(1).CFINV=gpuArray(D(1).CFINV);
  McScale=gpuArray(McScale);
  MpScale=gpuArray(MpScale);
  LO_Mc=gpuArray(LO_Mc);
  UP_Mc=gpuArray(UP_Mc);
  RWD=gpuArray(RWD);
else
  fprintf('USE CPU Max Chain=%4d Nitr=%2d Mc=%4d Mp=%3d \n',...
            PRM.CHA,PRM.ITR,Mc.N,Mp.N);
  fprintf(logFID,'USE CPU Max Chain=%4d Nitr=%2d Mc=%4d Mp=%3d \n',...
            PRM.CHA,PRM.ITR,Mc.N,Mp.N);
end
%
RT=0;
COUNT=0;
Burn=1;
Parallel=5;
%
while not(COUNT==PRM.THR)
  RT  =RT+1;
  NACC=0;tic
  rMc=zeros(Parallel*Mc.M,PRM.CHA);
  rMp=zeros(Parallel*Mp.M,PRM.CHA);
  rMi=zeros(Parallel*Mi.M,PRM.CHA);
  rLa=zeros(Parallel*La.M,PRM.CHA);
  if PRM.GPU~=99
    logU=log(rand(PRM.CHA,1,precision,'gpuArray'));
    for PT=1:Parallel
      rMc(Mc.N*(PT-1)+1:Mc.N*PT,:)=random('Normal',0,1,Mc.N,PRM.CHA,precision,'gpuArray');
      rMc(PT).ch=random('Normal',0,1,Mc.N,PRM.CHA,precision,'gpuArray');
      rMp(PT).ch=random('Normal',0,1,Mp.N,PRM.CHA,precision,'gpuArray');
      rMi(PT).ch=random('Normal',0,1,Mi.N,PRM.CHA,precision,'gpuArray');
      rLa(PT).ch=random('Normal',0,1,La.N,PRM.CHA,precision,'gpuArray');
    end
  else
    logU=log(rand(PRM.CHA,1,precision));
    for PT=1:Parallel
      rMc(PT).ch=random('Normal',0,1,Mc.N,PRM.CHA,precision);
      rMp(PT).ch=random('Normal',0,1,Mp.N,PRM.CHA,precision);
      rMi(PT).ch=random('Normal',0,1,Mi.N,PRM.CHA,precision);
      rLa(PT).ch=random('Normal',0,1,La.N,PRM.CHA,precision);
    end
  end
  rMp(find(POL.ID),:)=0;
  rMi(find(~BLK(1).IDinter),:)=0;
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
    Mi.SMP=Mi.OLD+RWD.*MiScale.*rMi(:,iT);
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
% CALC APRIORI AND RESIDUAL COUPLING RATE SECTION
    CAL.RIG=G.P*Mp.SMP;
    CAL.ELA=G.C*((G.TB*Mp.SMP).*D(1).CFINV.*Mc.SMPMAT);
    CAL.INE=G.I*Mi.SMP;
%     CAL.SMP=CAL.RIG+CAL.ELA;
    CAL.SMP=CAL.RIG+CAL.ELA+CAL.INE;   % including internal deformation
    if PRM.GPU~=99
      clear('CAL.RIG','CAL,ELA','CAL.INE');
    end
%   CAL.SMP=G.C*((G.TB*Mp.SMP).*Mc.SMPMAT)+G.P*Mp.SMP;       
%   CAL.SMP=G.P*Mp.SMP;
% CALC RESIDUAL SECTION
    RES.SMP=sum(((D(1).OBS-CAL.SMP)./D(1).ERR).^2,1);
% Mc is better Zero 
%     PRI.SMP=sum(abs(Mc.SMP),1);   
%% MAKE Probably Density Function
% $$ PDF_{post}=\frac{\frac{1}{\sqrt{2\pi\exp(L)}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}}{2}}\exp{\frac{-M^{2}}{2\times\exp{L}}}{\frac{1}{\sqrt{2\pi\exp(L_{old})}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}_{old}}{2}}\exp{\frac{-M^{2}_{old}}{2\times\exp{L_{old}}}} $$%%
%  log(x(x>0));
%   q1 = logproppdf(x0,y);
%   q2 = logproppdf(y,x0);
% this is a generic formula.
%   rho = (q1+logpdf(y))-(q2+logpdf(x0));  
    Pdf = -0.5.*...
         ((RES.SMP+La.SMP+exp(-La.SMP))...
         -(RES.OLD+La.OLD+exp(-La.OLD)));
%   Pdf = -0.5.*(RES.SMP-RES.OLD);
    ACC=Pdf > logU(iT);
    if ACC
      Mc.OLD  = Mc.SMP;
      Mp.OLD  = Mp.SMP;
      Mi.OLD  = Mi.SMP;
      La.OLD  = La.SMP;
      RES.OLD = RES.SMP;
%       PRI.OLD = PRI.SMP;
    end
% KEEP SECTION
    if iT > PRM.CHA-PRM.KEP
      if PRM.GPU~=99
        CHA.Mc(:,iT-(PRM.CHA-PRM.KEP))=gather(Mc.SMP);
        CHA.Mp(:,iT-(PRM.CHA-PRM.KEP))=gather(Mp.SMP);
        CHA.Mi(:,iT-(PRM.CHA-PRM.KEP))=gather(Mi.SMP);
        CHA.La(:,iT-(PRM.CHA-PRM.KEP))=gather(La.SMP);
      else
        CHA.Mc(:,iT-(PRM.CHA-PRM.KEP))=Mc.SMP;
        CHA.Mp(:,iT-(PRM.CHA-PRM.KEP))=Mp.SMP;
        CHA.Mi(:,iT-(PRM.CHA-PRM.KEP))=Mi.SMP;
        CHA.La(:,iT-(PRM.CHA-PRM.KEP))=La.SMP;
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
  Mi.STD=std(CHA.Mi,1,2); 
  La.STD=std(CHA.La,1,2);
%
  fprintf('T=%3d Res=%6.3f Accept=%5.1f RWD=%5.2f Time=%5.1fsec\n',...
           RT,1-RES.OLD./RR,100*CHA.AJR,RWD,toc)
  fprintf(logFID,'T=%3d Res=%6.3f Accept=%5.1f RWD=%5.2f Time=%5.1fsec\n',...
           RT,1-RES.OLD./RR,100*CHA.AJR,RWD,toc);
%
  for BK=1:BLK(1).NBlock
    [latp,lonp,ang]=xyzp2lla(CHA.Mp(3.*BK-2,:),CHA.Mp(3.*BK-1,:),CHA.Mp(3.*BK,:));
    fprintf('POLE OF BLOCK %2i = lat:%7.2f deg. lon:%8.2f deg. ang:%9.2e deg./m.y. \n',...
      BK,mean(latp),mean(lonp),mean(ang));
    fprintf(logFID,'POLE OF BLOCK %2i = lat:%7.2f deg. lon:%8.2f deg. ang:%9.2e deg./m.y. \n',...
      BK,mean(latp),mean(lonp),mean(ang));
  end
  fprintf('Lamda = %7.2f \n',mean(CHA.La));
  fprintf(logFID,'Lamda = %7.2f \n',mean(CHA.La));
%
  if Burn==0
    if CHA.AJR > 0.24
      RWD=RWD*1.1;
    elseif CHA.AJR < 0.22
      RWD=RWD*0.9;
    end
    COUNT=COUNT+1;
  else
    if CHA.AJR > 0.24
      RWD=RWD*1.1;
    elseif CHA.AJR < 0.22
      RWD=RWD*0.9;
    else
      COUNT=COUNT+1;
      Burn=0;
    end
  end
%   if CHA.AJR > 0.24
%     RWD=RWD*1.1;
%     COUNT=0;
%   elseif CHA.AJR < 0.22
%     RWD=RWD*0.9;
%     COUNT=0;
%   else
%     COUNT=COUNT+1;
%   end
  CHA.SMP=CAL.SMP;
  % debug-----------
  Mpmean=mean(CHA.Mp,2);
  Mcmean=mean(CHA.Mc,2);
  Mimean=mean(CHA.Mi,2);
  Mcmeanrep=repmat(Mcmean,3,D.CNT);Mcmeanrep=Mcmeanrep(D.MID);
  VEC.RIG=G.P*Mpmean;
  VEC.ELA=G.C*((G.TB*Mpmean).*D(1).CFINV.*Mcmeanrep);
  VEC.INE=G.I*Mimean;
%   VEC.SUM=VEC.RIG+VEC.ELA;
  VEC.SUM=VEC.RIG+VEC.ELA+VEC.INE;   % including internal deformation
%   vec.rel=G.C*((G.TB*poltmp).*CF);
  % debug-----------
  if PRM.GPU~=99
    cCHA.Mc=gather(CHA.Mc);
    cCHA.Mp=gather(CHA.Mp);
    cCHA.Mi=gather(CHA.Mi);
    cCHA.La=gather(CHA.La);
    cCHA.SMP=gather(CHA.SMP);
    MAKE_FIG(cCHA,BLK,OBS,RT,gather(LO_Mc),gather(UP_Mc),VEC,Mimean)
  else
    MAKE_FIG(CHA,BLK,OBS,RT,LO_Mc,UP_Mc,VEC,Mimean)
  end
  if RT > PRM.ITR; break; end;
end
if PRM.GPU~=99
  CHA.Mc=gather(CHA.Mc);
  CHA.Mp=gather(CHA.Mp);
  CHA.Mi=gather(CHA.Mi);
  CHA.La=gather(CHA.La);
  CHA.SMP=gather(CHA.SMP);
end
CHA.Res=RES.SMP;
fprintf('RMS=: %8.3f\n',CHA.Res)
fprintf('=== FINISHED MH_MCMC ===\n')
fprintf(logFID,'RMS=: %8.3f\n',CHA.Res);
fprintf(logFID,'=== FINISHED MH_MCMC ===\n');
fclose(logFID);
end
