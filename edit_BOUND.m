function edit_BOUND
BLK=READ_BLOCK_BOUND('./BLOCK/');
SHOW_BLOCK_BOUND(BLK);
REDUCE_BLOCK_BOUND(BLK,3,5,20);
SHOW_BLOCK_BOUND(BLK);
%WRITE_BLOCK_BOUND(BLK_OUT,'./BLOCK_OUT/');
end
%====================================================
function OUT_BLK=REDUCE_BLOCK_BOUND(BLK,NB1,NB2,INT)
B1.LON=BLK(NB1).LON;
B1.LAT=BLK(NB1).LAT;
B2.LON=BLK(NB2).LON;
B2.LAT=BLK(NB2).LAT;
B.LON=BLK(1).BOUND(NB1,NB2).LON;
B.LAT=BLK(1).BOUND(NB1,NB2).LAT;
ALON=B.LON(1);
ALAT=B.LAT(1);
[B.X  ,B.Y  ]=PLTXY(B.LAT,B.LON,ALAT,ALON);
(ALAT,ALON);


[B.LAT,B.LON]=XYTPL(B.X  ,B.Y  ,ALAT,ALON);
OUT_BLK=BLK;
OUT_BLK(NB1).LON=B1.LON;
OUT_BLK(NB1).LAT=B1.LAT;
OUT_BLK(NB2).LON=B2.LON;
OUT_BLK(NB2).LAT=B2.LAT;
OUT_BLK(1).BOUND(NB1,NB2).LON=B.LON;
OUT_BLK(1).BOUND(NB1,NB2).LAT=B.LAT;
end
%====================================================
function SHOW_BLOCK_BOUND(BLK)
figure('BLOCK_AND_BOUNDARY_MAP')
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
end
%====================================================
function BLK=READ_BLOCK_BOUND(DIR)
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
fprintf('READ BLOCK FILES : %4i \n',BLK(1).NBlock)
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    BLK(1).BOUND(NB1,NB2).LAT=[];
    BLK(1).BOUND(NB1,NB2).LON=[];
    [~,ialon,~]=intersect(BLK(NB1).LON,BLK(NB2).LON);
    [~,ialat,~]=intersect(BLK(NB1).LAT,BLK(NB2).LAT);
    [Ca,~,~]=intersect(ialat,ialon);
    BLK(1).BOUND(NB1,NB2).LAT=BLK(NB1).LAT(Ca);
    BLK(1).BOUND(NB1,NB2).LON=BLK(NB1).LON(Ca);
    if sum(Ca) > 0 
      fprintf('Number of boundary points %4i between %2i and %2i \n',length(Ca),NB1,NB2)
    end
  end
end
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
