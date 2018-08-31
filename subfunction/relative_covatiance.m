% function TCHA2 = relative_covatiance(TCHA,BLK,TRI,G)
function TCHA2 = relative_covatiance(DIR)
load(fullfile(DIR,'TCHA.mat'))
load(fullfile(DIR,'BLK.mat'))
load(fullfile(DIR,'GRN.mat'))
load(fullfile(DIR,'TRI.mat'))
slipvel=zeros(size(TCHA.SMPFLT));
MC=1;
MR=1;
relativevel=full(G.TB)*TCHA.SMPPOL;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(TRI(1).BOUND(NB1,NB2).clon,2);
    if NF~=0
      fprintf('%i and %i\n',NB1,NB2);
      velstr=relativevel(MC     :MC+  NF-1,:);
      veldip=relativevel(MC+  NF:MC+2*NF-1,:);
      veltns=relativevel(MC+2*NF:MC+3*NF-1,:);
      slip=sqrt(velstr.^2+veldip.^2+veltns.^2);
      slipvel(MR:MR+NF-1,:)=slip;
      clear velstr veldip veltns slip
      MC=MC+3*NF;
      MR=MR+NF;
    end
  end
end
clear relativevel

div=100;
MD=1;
sumsmp=0;
sumflt=0;
sumvel=0;
NS=size(TCHA.SMPFLT,2);
ns=NS/div;
for RP=1:div
    smp=[double(TCHA.SMPFLT(:,MD:MD+ns-1));slipvel(:,MD:MD+ns-1)];
    sumsmp=sumsmp+smp*smp';
    sumflt=sumflt+sum(double(TCHA.SMPFLT(:,MD:MD+ns-1)),2);
    sumvel=sumvel+sum(slipvel(:,MD:MD+ns-1),2);
    MD=MD+ns;
end
avesmp=[sumflt;sumvel]./NS;
covsmp=sumsmp./NS-avesmp*avesmp';
varsmp=diag(covsmp);
varsmp(varsmp<0)=0;
stdsmp=sqrt(varsmp);
corsmp=covsmp./(stdsmp*stdsmp');

TCHA2.slipvel=slipvel;
clear slipvel
polllw=euler_llw(TCHA,BLK);
TCHA2.polllw=polllw;
clear polllw
TCHA2.covsmp=covsmp;
TCHA2.corsmp=corsmp;
clear covsmp corsmp

save(fullfile(DIR,'TCHA2.mat'),'TCHA2','-v7.3')

end
%%
function polllw = euler_llw(TCHA,BLK)
MC=1;
polllw=zeros(size(TCHA.SMPPOL));
for NB=1:BLK(1).NBlock
  X=TCHA.SMPPOL(MC  ,:);
  Y=TCHA.SMPPOL(MC+1,:);
  Z=TCHA.SMPPOL(MC+2,:);
  lat=atan2(Z,sqrt(X.*X+Y.*Y)).*180/pi;
  lon=atan2(Y,X).*180/pi;
  ang=sqrt(X.*X+Y.*Y+Z.*Z).*(1e6.*(180./pi));
  polllw(MC:MC+2,:)=[lat;lon;ang];
  MC=MC+3;
end

end