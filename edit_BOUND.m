function edit_BOUND
BLK=READ_BLOCK_BOUND('./BLOCK/');
PAR=READ_PARAMETERS('./opt_bound_par.txt');
for nB=1:PAR(1).num
  BLK=OPT_BLOCK_BOUND(BLK,PAR(1).B1(nB),PAR(1).B2(nB),PAR(1).INT(nB));
end
SHOW_BLOCK_BOUND(BLK);
WRITE_BLOCK_BOUND(BLK,'./BLOCK_OUT/');
end
%====================================================
function WRITE_BLOCK_BOUND(BLK,oDIR)
nBLK=BLK(1).NBlock;
[NBlock,~]=size(file);
BLK(1).NBlock=NBlock;
for NB=1:nBLK
  fname=strcat(sprintf('%02i',NB),'_',BLK(NB).name,'.txt');
  fullname=fullfile(oDIR,fname);
  fileID=open(fullname); 
  fprintf(fileID,'%15.9f %15.9f/n',BLK(NB).LAT,BLK(NB).LON);
  fclose(fileID);
end
end
%====================================================
function OUT_BLK=OPT_BLOCK_BOUND(BLK,NB1,NB2,INT)
dint=0.01;
deps=dint.*INT;
%
B1.LON=BLK(NB1).LON;
B1.LAT=BLK(NB1).LAT;
B2.LON=BLK(NB2).LON;
B2.LAT=BLK(NB2).LAT;
%
B.LON=BLK(1).BOUND(NB1,NB2).LON;
B.LAT=BLK(1).BOUND(NB1,NB2).LAT;
ALAT=B.LAT(1); ALON=B.LON(1);
[B.X,B.Y]=PLTXY(B.LAT,B.LON,ALAT,ALON);
%
count=0; oBXY=[B.X B.Y]; BXY=oBXY; NBO=size(BXY,1);
fprintf('BOUNDARY: %2i %2i INT:%5.1f intPOINT:%5i \n',NB1,NB2,INT,NBO)
%
[B.LON(1) B.LAT(1) B.LON(end) B.LAT(end)]
Bind.S(1)=find((B1.LAT==B.LAT(1))  &(B1.LON==B.LON(1)));
Bind.E(1)=find((B1.LAT==B.LAT(end))&(B1.LON==B.LON(end)));
Bind.S(2)=find((B2.LAT==B.LAT(1))  &(B2.LON==B.LON(1)));
Bind.E(2)=find((B2.LAT==B.LAT(end))&(B2.LON==B.LON(end)));
%
while 1
  count=count+1;
  dXY=[diff(BXY(:,1)) diff(BXY(:,2))];
  dL=sqrt(dXY(:,1).^2+dXY(:,2).^2);
  mdL=mean(dL);
%
  if mod(count,10)==0 || mdL < 0.25*INT
    if count > 10*NBO; break; end
%
    if (std(dL-INT) < 0.1*INT) && (abs(mdL-INT) < 0.1*INT && max(dist_bo(BXY,oBXY)) < 0.01*INT )
      fprintf('Count:%5i Mean(dL):%6.1f POINT:%5i STD:%5.1f \n',count,mdL,size(BXY,1),std(dL-INT))
      OUT_BLK=update_BLK(BXY,ALAT,ALON,NB1,NB2,B1,B2,BLK,Bind,B);
      break;
    end
%
    if mdL > INT
      [~,index]=max(dL);
      BXY=[BXY(1:index,:);mean(BXY(index:index+1,:),1); BXY(index+1:end,:)];
    else
      [~,index]=min(dL);
      if index==length(dL)
        BXY=[BXY(1:index-1,:); BXY(        end,:)];
      else
        BXY=[BXY(1:index  ,:); BXY(index+2:end,:)];
      end
    end
%
    if mod(count,100)==0 
      fprintf('Count:%5i Mean(dL):%6.1f POINT:%5i STD:%5.1f \n',count,mdL,size(BXY,1),std(dL-INT))
      OUT_BLK=update_BLK(BXY,ALAT,ALON,NB1,NB2,B1,B2,BLK,Bind,B);
      SHOW_BLOCK_BOUND(OUT_BLK);
    end
    continue
  end
 
%  
  ind=dL>INT;
  ind1=[false(1);ind];
  ind2=[ind;false(1)];
  BXY(ind1,:)=BXY(ind1,:)-dint.*dXY(ind,:);
  BXY(ind2,:)=BXY(ind2,:)+dint.*dXY(ind,:);
%
  dd = dist_bo(BXY,oBXY);
  gX =(dist_bo([BXY(:,1)+deps BXY(:,2)],oBXY)-dd)/deps;
  gY =(dist_bo([BXY(:,1) BXY(:,2)+deps],oBXY)-dd)/deps;
  gXY=dd./sqrt(gX.^2+gY.^2);
  BXY=BXY-gXY.*[gX gY];
  BXY(1,:)=oBXY(1,:); BXY(end,:)=oBXY(end,:);
%
end
end
%====================================================
function d=dist_bo(po,bo)
% Distance of bo (line) and po (point)  
npo=length(po);
d=inf(npo,1);
for n=1:npo
  [~,ind]=sort(sqrt((bo(:,1)-po(n,1)).^2+(bo(:,2)-po(n,2)).^2));
  in1=ind(1); in2=ind(2);
  a=bo(in2,2)-bo(in1,2); b=bo(in2,1)-bo(in1,1);
  d(n)=abs(a.*po(n,1)-b.*po(n,2)-a.*bo(in1,1)+b.*bo(in1,2))./sqrt(a.^2+b.^2);
end
end
%====================================================
function BLK=update_BLK(BXY,ALAT,ALON,NB1,NB2,B1,B2,BLK,Bind,B)
[LAT,LON]=XYTPL(BXY(:,1),BXY(:,2),ALAT,ALON);
BLK(NB1).LON=[B1.LON(1:Bind.S(1)-1);LON;B1.LON(Bind.E(1)+1:end)];
BLK(NB1).LAT=[B1.LAT(1:Bind.S(1)-1);LAT;B1.LAT(Bind.E(1)+1:end)];
BLK(NB2).LON=[B2.LON(1:Bind.S(2)-1);LON;B2.LON(Bind.E(2)+1:end)];
BLK(NB2).LAT=[B2.LAT(1:Bind.S(2)-1);LAT;B2.LAT(Bind.E(2)+1:end)];
BLK(1).BOUND(NB1,NB2).LON=[B.LON(1);LON(2:end-1);B.LON(end)]; 
BLK(1).BOUND(NB1,NB2).LAT=[B.LAT(1);LAT(2:end-1);B.LAT(end)];
end
%====================================================
function SHOW_BLOCK_BOUND(BLK)
figure(100);
clf
%figure(h,'Name','BLOCK_AND_BOUNDARY_MAP');
for NB=1:BLK(1).NBlock
  plot(BLK(NB).LON,BLK(NB).LAT)
  hold on
  text(mean(BLK(NB).LON),mean(BLK(NB).LAT),int2str(NB))
  hold on
end
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    plot(BLK(1).BOUND(NB1,NB2).LON,BLK(1).BOUND(NB1,NB2).LAT,'o')
    hold on
  end
end
drawnow
end
%====================================================
function BLK=READ_BLOCK_BOUND(DIR)
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
%figure('Name','BLOCK_BOUNDARY_LINE')
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
      fprintf('BLOCK BOUNDARY : %2i %2i POINTS : %5i \n',NB1,NB2,size(BLK(1).BOUND(NB1,NB2).LAT,1))
%      plot(BLK(1).BOUND(NB1,NB2).LON,BLK(1).BOUND(NB1,NB2).LAT)
%      hold on
    end
  end
end
end
%====================================================
function [LAT,LON]mach_bo(BLK,NB1,NB2)
LCa=inpolygon(BLK(NB1).LON,BLK(NB1).LAT,BLK(NB2).LON,BLK(NB2).LAT);
Ca=find(LCa);
if ~isempty(Ca)
  if and(LCa(1),LCa(end))
    Ca0=find(LCa~=true,1,'last')+1:length(LCa)-1;
    Ca1=1:find(LCa~=true,1,'first')-1;
    Ca=[Ca0 Ca1];
  end
LAT=BLK(NB1).LAT(Ca);
LON=BLK(NB1).LON(Ca);
end
%====================================================
function [X,Y]=PLTXY(ALAT,ALON,ALAT0,ALON0)
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
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
function PAR=READ_PARAMETERS(infile)
tmp=load(infile);
PAR(1).num=size(tmp,1);
PAR(1).B1=tmp(:,1);
PAR(1).B2=tmp(:,2);
PAR(1).INT=tmp(:,3);
end