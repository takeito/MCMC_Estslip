function dist=distance_matrix(DIR)
% Compute the distance of each subfault.
load(fullfile(DIR,'BLK.mat'))
load(fullfile(DIR,'TRI.mat'))

MC=1;
MT=1;
MR=1;
clon=zeros(BLK(1).NB,1);
clat=zeros(BLK(1).NB,1);
cdep=zeros(BLK(1).NB,1);
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(TRI(1).BOUND(NB1,NB2).clon,2);
    if NF~=0
      clon(MR:MR+NF-1,1)=TRI(1).BOUND(NB1,NB2).clon;
      clat(MR:MR+NF-1,1)=TRI(1).BOUND(NB1,NB2).clat;
      cdep(MR:MR+NF-1,1)=TRI(1).BOUND(NB1,NB2).cdep;
      MC=MC+3*NF;
      MT=MT+2*NF;
      MR=MR+  NF;
    end
  end
end
[cx,cy,cz]=ell2xyz(clat,clon,1e3.*cdep);
distx2=cx*cx';
disty2=cy*cy';
distz2=cz*cz';
dist=sqrt(repmat(diag(distx2),1,BLK(1).NB) - 2.*distx2 + repmat(diag(distx2)',BLK(1).NB,1)...
         +repmat(diag(disty2),1,BLK(1).NB) - 2.*disty2 + repmat(diag(disty2)',BLK(1).NB,1)...
         +repmat(diag(distz2),1,BLK(1).NB) - 2.*distz2 + repmat(diag(distz2)',BLK(1).NB,1));
dist=dist*1e-3; % m -> km

save(fullfile(DIR,'dist.mat'),'dist','-v7.3')
end
%%
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
