function make_interface_tri
% PRM.OBS_F='./data_set.txt';
PRM.OBS_F='./geonet_jcg_nu.txt';
PRM.SUB_F='./plate_phs.txt';
% PRM.BOU_F='./bound.txt';
PRM.BOU_F='./phcont.txt';
PRM.INIP_F='./initial_point.txt';
savefolder='./figs/Mesh';
anim_savefile='./figs/tri_anim.gif';

% ------------------Parameter---------------------
PRM.NMESH=1000;
Reducerate=0.001; % Reduce rate of triangles
ARR=0.4;                                      % Reduce badangle triangle rate in each roop
ARS=0.5;                                      % Reduce badangle when roop reached to PRM.NMESH*AngRedStart
pole.lon=161.65; pole.lat=54.74; pole.omega=-1.168;   % PH plate motion relative to AM plate(REVEL)
% --------------------------------------------------

OBS=READ_OBS(PRM);
ini_size=FIX_POINT(PRM.INIP_F);
s=INIT_INTERFACE_TRI(PRM.SUB_F, PRM.BOU_F, PRM.INIP_F, PRM.NMESH.*5,ini_size);
save('./figs/plate_phs_initial','s')

s=DOWN_TRI(s,OBS,PRM.NMESH,savefolder,anim_savefile,Reducerate,ini_size,pole,ARS,ARR);

[tarea,tcenter]=calc_tri_area(s);
save(['./figs/area_tri_',num2str(PRM.NMESH)],'tarea','tcenter');
% save('plate_phs','s')
save(['./figs/plate_phs',num2str(PRM.NMESH)],'s')
% [tarea,tcenter]=calc_tri_area(s);
% save(['./figs/area_tri_',num2str(PRM.NMESH)],'tarea','tcenter');


end

%% READ OBSERVATION DATA
function [OBS]=READ_OBS(PRM)
% moditied 2015/11/11 by T.ITO
%-------------------
% format:
% site_name lon lat EW_comp. NS_comp. UD_comp. ERR_EW ERR_NS ERR_UD 
%-------------------
Fid_OBS=fopen(PRM.OBS_F,'r');
N=0;
while 1
  tline=fgetl(Fid_OBS);
  if ~ischar(tline); break; end
  str=strsplit(tline);
  N=N+1;
% TODO: FOR TEST READ TWO BLOCK REGION
%   if N<=194
%     OBS(1).LAT(N) =str2double(cellstr(str(3))); %LAT
%     OBS(1).LON(N) =str2double(cellstr(str(2))); %LON
%     OBS(1).HIG(N) =str2double(cellstr(str(4))); %HIG
%   else
%     OBS(2).LAT(N-194) =str2double(cellstr(str(3))); %LAT
%     OBS(2).LON(N-194) =str2double(cellstr(str(2))); %LON
%     OBS(2).HIG(N-194) =str2double(cellstr(str(4))); %HIG
%   end
%   OBS(1).NAME(N)=cellstr(str(1));
%   OBS(1).ALAT(N) =str2double(cellstr(str(3))); %LAT
%   OBS(1).ALON(N) =str2double(cellstr(str(2))); %LON
%   OBS(1).AHIG(N) =str2double(cellstr(str(4))); %HIG
%   OBS(1).EVEC(N) =str2double(cellstr(str(5))); %E-W
%   OBS(1).NVEC(N) =str2double(cellstr(str(6))); %N-S
%   OBS(1).HVEC(N) =str2double(cellstr(str(7))); %U-D
%   OBS(1).YYY(3*N-2) =str2double(cellstr(str(5))); %E-W
%   OBS(1).YYY(3*N-1) =str2double(cellstr(str(6))); %N-S
%   OBS(1).YYY(3*N)   =str2double(cellstr(str(7))); %U-D
%   OBS(1).EEE(3*N-2) =str2double(cellstr(str(8))); %E-W
%   OBS(1).EEE(3*N-1) =str2double(cellstr(str(9))); %N-S
%   OBS(1).EEE(3*N)   =str2double(cellstr(str(10))); %U-D
  OBS(1).ALAT(N) =single(str2double(cellstr(str(3)))); %LAT
  OBS(1).ALON(N) =single(str2double(cellstr(str(2)))); %LON
  OBS(1).AHIG(N) =single(str2double(cellstr(str(4)))); %HIG
  OBS(1).EVEC(N) =single(str2double(cellstr(str(5)))); %E-W
  OBS(1).NVEC(N) =single(str2double(cellstr(str(6)))); %N-S
  OBS(1).HVEC(N) =single(str2double(cellstr(str(7)))); %U-D
  OBS(1).YYY(3*N-2) =single(str2double(cellstr(str(5)))); %E-W
  OBS(1).YYY(3*N-1) =single(str2double(cellstr(str(6)))); %N-S
  OBS(1).YYY(3*N)   =single(str2double(cellstr(str(7)))); %U-D
  OBS(1).EEE(3*N-2) =single(str2double(cellstr(str(8)))); %E-W
  OBS(1).EEE(3*N-1) =single(str2double(cellstr(str(9)))); %N-S
  OBS(1).EEE(3*N)   =single(str2double(cellstr(str(10)))); %U-D
end
fprintf('==================\nNumber of observation site : %i \n',N)
end
%====================================================
function ininum = FIX_POINT(inip_f)

% COUNT THE NUMBER OF INITIAL FIXED POINTS

Fid=fopen(inip_f);
ini_point=textscan(Fid,'%f%f%f');
fclose(Fid);
ini_point=cell2mat(ini_point);
[ininum,~]=size(ini_point);

end
%====================================================
function [S]=DOWN_TRI(S,OBS,n_mesh,saveloc,animfile,Redu_rate,initial,POLE,ars,arr)
alat0=(mean(OBS(1).ALAT)+mean(S.lat))./2;
alon0=(mean(OBS(1).ALON)+mean(S.lon))./2;
tg.lon=mean(S.lon(S.tri),2);
tg.lat=mean(S.lat(S.tri),2);
tg.dep=mean(S.dep(S.tri),2);
[tgxyz]=enu2xyz(tg.lon,tg.lat,6371D0);
[pxyz]=enu2xyz(POLE.lon,POLE.lat,deg2rad(POLE.omega));
[PMEN]=platemotion(tgxyz.x,tgxyz.y,tgxyz.z,tg.lon,tg.lat,pxyz.x,pxyz.y,pxyz.z);
[gx,gy]=PLTXY(OBS(1).ALAT,OBS(1).ALON,alat0,alon0);
[sx,sy]=PLTXY(S.lat,S.lon,alat0,alon0);
gz=OBS(1).AHIG./1000;
sz=S.dep;
Ua=zeros(length(S.tri(:,1)),1);
minAng=zeros(length(S.tri),1);
tri_Ang=zeros(length(S.tri),3);
parfor n=1:length(S.tri)
% for n=1:length(S.tri)        % <<--------- use in debug test
    [SDT]=SDTvec(sx(S.tri(n,:)),sy(S.tri(n,:)),sz(S.tri(n,:)));
    PMcom=[PMEN.EW(n) PMEN.NS(n) 0]*SDT;
    PMcom=PMcom./norm(PMcom);
    [Ang]=triangle_angles([sx(S.tri(n,:)) sy(S.tri(n,:)) sz(S.tri(n,:))],'d');
    tri_Ang(n,:)=Ang(:);
    minAng(n,1)=min(Ang);    
    [U]=CalcTriDisps(gx',gy',gz',sx(S.tri(n,:)),sy(S.tri(n,:)),sz(S.tri(n,:)),0.25,PMcom(1),PMcom(3),PMcom(2));
    Ua(n)=sum(sqrt(U.x.^2+U.y.^2+U.z.^2));
end
clear n;
Ntri=length(S.tri);
Ntriini=Ntri;

ffanim=0;
while Ntri > n_mesh
  Redu_tri=ceil(Ntri*Redu_rate);
  r_index=ones(length(S.lat),1);
  min_tri=[];
%   [~,index]=min(Ua);  
%   min_tri=S.tri(index,:);
%   f_tri=zeros(3,1);
%   for n=1:3
%       ftrisum=find(S.tri,min_tri(n)); %%
%       f_tri(n)=sum(find(S.tri,min_tri(n)));
%   end
  [~,Uasortindex]=sort(Ua);
  f_tri=zeros(Redu_tri,3);
  
% Choice triangles apart of triangles which consisit of only initial edge point
  Ualogic = S.tri(Uasortindex,1)>initial | S.tri(Uasortindex,2)>initial | S.tri(Uasortindex,3)>initial;
  Uasortindex=Uasortindex(Ualogic);
  
  if Ntri > Ntriini*ars
      min_tri=S.tri(Uasortindex(1:Redu_tri),:);
      for mm=1:Redu_tri
          for n=1:3
              if min_tri(mm,n)<=initial
                  f_tri(mm,n)=0;              % If the vertex of triangle is initial edge point, don'n remove this point.
              else
                  f_tri(mm,n)=tri_Ang(Uasortindex(mm),n);
%                   [f_tri(mm,n),~]=size(find(S.tri==min_tri(mm,n)));
              end
          end
          clear n;
          [~,index]=max(f_tri(mm,:));
          r_index(min_tri(mm,index))=0;
      end
  else
      [~,Angsortindex]=sort(minAng);
      Anglogic=S.tri(Angsortindex,1)>initial | S.tri(Angsortindex,2)>initial | S.tri(Angsortindex,3)>initial;
      Angsortindex=Angsortindex(Anglogic);
      min_tri=S.tri(Uasortindex(1:ceil(Redu_tri*arr)),:);
%       min_tri=S.tri(Angsortindex(1:ceil(Redu_tri*angredurate)),:);
      ang_tri=S.tri(Angsortindex(1:Redu_tri-ceil(Redu_tri*arr)),:);
%       min_tri(ceil(Redu_tri*arr)+1:Redu_tri,:)=S.tri(Uasortindex(1:Redu_tri-ceil(Redu_tri*angredurate)),:);
      for mm=1:ceil(Redu_tri*arr)
          for n=1:3
              if min_tri(mm,n)<=initial
                  f_tri(mm,n)=0;              % If the vertex of triangle is initial edge point, don'n remove this point.
              else
                  f_tri(mm,n)=tri_Ang(Uasortindex(mm),n);
%                   [f_tri(mm,n),~]=size(find(S.tri==min_tri(mm,n)));
              end
          end
          clear n;
          [~,index]=max(f_tri(mm,:));
          r_index(min_tri(mm,index))=0;
      end
      for mm=1:Redu_tri-ceil(Redu_tri*arr)
          for n=1:3
              if ang_tri(mm,n)<=initial
                  f_tri(mm,n)=0;              % If the vertex of triangle is initial edge point, don'n remove this point.
              else
                  f_tri(mm,n)=tri_Ang(Angsortindex(mm),n);
              end
          end
          clear n;
          [~,index]=max(f_tri(mm,:));
          r_index(ang_tri(mm,index))=0;
      end
  end
  clear mm;
  r_index=logical(r_index);
%   min_tri=S.tri(Uasortindex(1:Redu_tri),:);
  
%   for mm=1:Redu_tri
%       for n=1:3
%           if min_tri(mm,n)<=initial
%               f_tri(mm,n)=0;              % If the vertex of triangle is initial edge point, don'n remove this point.
%           else
%               [f_tri(mm,n),~]=size(find(S.tri==min_tri(mm,n)));
%           end
%       end
%       clear n;
%       [~,index]=max(f_tri(mm,:));
%       r_index(min_tri(mm,index))=0;
%   end
%   clear mm;
%   r_index=logical(r_index);
  
%   f_tri=find(S.tri,min_tri)  
%   [~,index]=max(f_tri);
%   r_index(min_tri(index))=0;
%   r_index=logical(r_index);
  S.lat=S.lat(r_index);
  S.lon=S.lon(r_index);
  S.dep=S.dep(r_index);
  S.tri=delaunay(S.lon,S.lat);
  [sx,sy]=PLTXY(S.lat,S.lon,alat0,alon0);
  sz=S.dep;
%---------
  ntri=length(S.tri);
  Stri=[];
%   nn=0;
%   for n=1:ntri
%     glon=mean(S.lon(S.tri(n,:)));  
%     glat=mean(S.lat(S.tri(n,:)));
%     ID=inpolygon(glon,glat,S.bound(:,1),S.bound(:,2));
%     if ID==1
%       nn=nn+1;
%       Stri(nn,:)=S.tri(n,:);
%     end  
%   end

  glon=mean(S.lon(S.tri),2);
  glat=mean(S.lat(S.tri),2);
  ID=inpolygon(glon,glat,S.bound(:,1),S.bound(:,2));
  ID1ind=find(ID==1);
  [nn,~]=size(ID1ind);
  Stri=S.tri(ID1ind,:);
  tg.lon=mean(S.lon(Stri),2);
  tg.lat=mean(S.lat(Stri),2);
  [tgxyz]=enu2xyz(tg.lon,tg.lat,6371D0);
  [PMEN]=platemotion(tgxyz.x,tgxyz.y,tgxyz.z,tg.lon,tg.lat,pxyz.x,pxyz.y,pxyz.z);
  
%---------
  Ua_tmp=Ua;
  Ua=zeros(nn,1);
  
  tri_Ang_tmp=tri_Ang;
  minAng_tmp=minAng;
  minAng=zeros(length(Stri),1);
  tri_Ang=zeros(length(Stri),3);
  parfor n=1:nn
      %     for n=1:nn   %    <<----------------- use in debug test
      if Stri(n,1)~=S.tri(n,1) || Stri(n,2)~=S.tri(n,2) || Stri(n,3)~=S.tri(n,3)
          [SDT]=SDTvec(sx(Stri(n,:)),sy(Stri(n,:)),sz(Stri(n,:)));
          PMcom=[PMEN.EW(n) PMEN.NS(n) 0]*SDT;
          PMcom=PMcom./norm(PMcom);
          [Ang]=triangle_angles([sx(Stri(n,:)) sy(Stri(n,:)) sz(Stri(n,:))],'d');
          tri_Ang(n,:)=Ang(:);
          minAng(n,1)=min(Ang);
          %      [U]=CalcTriDisps(gx',gy',gz',sx(Stri(n,:)),sy(Stri(n,:)),sz(Stri(n,:)),0.25,0,0,1);
          [U]=CalcTriDisps(gx',gy',gz',sx(Stri(n,:)),sy(Stri(n,:)),sz(Stri(n,:)),0.25,PMcom(1),PMcom(3),PMcom(2));
          Ua(n)=sum(sqrt(U.x.^2+U.y.^2+U.z.^2));
      else
          Ua(n)=Ua_tmp(n);
          tri_Ang(n,:)=tri_Ang_tmp(n,:);
          minAng(n,1)=minAng_tmp(n,1);
      end
  end

  clear n;

  S.tri=Stri;
  Ntri=length([S.tri]);
  
  ffanim=ffanim+1;
% figure view and save each down_tri roop ------------
%     Fid=figure('visible','off');
%     plot(S.bound(:,1),S.bound(:,2),'r');
%     hold on;
%     plot(OBS(1).ALON,OBS(1).ALAT,'.g');
%     triplot(S.tri,S.lon,S.lat);
%     title(['Number of triangels= ',num2str(Ntri)]);
%     fprintf('Number of triangels=%4.0f \n ',Ntri)
%     print(Fid,'-depsc ',[saveloc,num2str(Ntri)]);
%     close(Fid)

% figure view only each down_tri roop ----------
%   figure(30); clf
%   plot(S.bound(:,1),S.bound(:,2),'r');
%   hold on;
%   plot(OBS(1).ALON,OBS(1).ALAT,'.g');
%   triplot(S.tri,S.lon,S.lat);
%   title(['Number of triangels= ',num2str(Ntri)]);
%   pause(.1)

% GIF animation test -----------
  if ffanim==1
      fig30=figure;
  else
      fig30=figure('visible','off');
  end
  plot(S.bound(:,1),S.bound(:,2),'r');
  hold on;
  plot(OBS(1).ALON,OBS(1).ALAT,'.g');
  triplot(S.tri,S.lon,S.lat);
  title(['Number of triangels= ',num2str(Ntri)]);
  if ffanim == 1;
      print('-depsc',[saveloc,num2str(Ntri)]);
  end  
  fprintf('Number of triangels=%4.0f \n',Ntri)
  frame=getframe(fig30);
  im=frame2im(frame);
  [A,map]=rgb2ind(im,256);
  if ffanim == 1;
      imwrite(A,map,animfile,'gif','LoopCount',Inf,'DelayTime',0.2);
  else
      imwrite(A,map,animfile,'gif','WriteMode','append','DelayTime',0.2);
  end
  hold off;clf;
end

% export when exit roop
S.x=sx;
S.y=sy;
S.z=sz;
S.U=Ua;

% save figure at the number of n-mesh
close(fig30)
clear ffanim;

Fid=figure('visible','off');
plot(S.bound(:,1),S.bound(:,2),'r');
hold on;
plot(OBS(1).ALON,OBS(1).ALAT,'.g');
triplot(S.tri,S.lon,S.lat);
title(['Number of triangels= ',num2str(Ntri)]);
fprintf('Number of triangels=%4.0f \n ',Ntri)
print(Fid,'-depsc ',[saveloc,num2str(Ntri)]);
close(Fid)

end
%====================================================
function [s]=INIT_INTERFACE_TRI(sub_f,bound_f,inip_f,int_mesh,ini_size)
%====================================================
Fid=fopen(sub_f);
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
%====================================================
F=scatteredInterpolant(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');
min_lon=min(bound(:,1)); max_lon=max(bound(:,1));
min_lat=min(bound(:,2)); max_lat=max(bound(:,2));
figure(10); clf
plot(bound(:,1),bound(:,2),'r')
hold on
% initial straint point-------
s.lon=ini_point(:,1);
s.lat=ini_point(:,2);
s.dep=ini_point(:,3);
%-----------------------------
% n=0;
n=ini_size;
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
% ntri=length(tri);
% nn=0;
% for n=1:ntri
%   glon=mean(s.lon(tri(n,:)));  
%   glat=mean(s.lat(tri(n,:)));
%   ID=inpolygon(glon,glat,bound(:,1),bound(:,2));
%   if ID==1
%     nn=nn+1;
%     s.tri(nn,:)=tri(n,:);
%   end  
% end
% clear n
glon=mean(s.lon(tri),2);
glat=mean(s.lat(tri),2);
ID=inpolygon(glon,glat,bound(:,1),bound(:,2));
ID1ind=find(ID==1);
s.tri=tri(ID1ind,:);

figure(20); clf
plot(bound(:,1),bound(:,2),'r')
hold on
triplot(s.tri,s.lon,s.lat)
pause(.1)
s.bound=bound;
s.dep_sub=dep_sub;
end
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
function [VEN] = platemotion(trigx,trigy,trigz,tlon,tlat,px,py,pz)

% Input-------------------------------------------------
%   trigx,trigy,trigz     :Center of triangle position
%   tlon,tlat             :Longitude and Latitude of center of triangle
%   px,py,pz              :Pole
% Output------------------------------------------------
%   V.EW, V.NS            :plate motion[km/yr] (Positive directions are E and N)
% ------------------------------------------------------

tlon=deg2rad(tlon);
tlat=deg2rad(tlat);

% Transformation from pole coord. to XYZ
[polexyz]=[px; py; pz];
[trigxyz]=[trigx'; trigy'; trigz'];

% Calculate plate motion
[~,trigsize]=size(trigxyz);
pm=zeros(3,trigsize);
for ii=1:trigsize
    pm(:,ii)=cross(polexyz,trigxyz(:,ii)).*1.0D-6;
end
% PM=PM';
Vx=pm(1,:);
Vy=pm(2,:);
Vz=pm(3,:);

% Export to E-N direction Velocity
VEN.EW=           -sin(tlon).*Vx' +            cos(tlon).*Vy';
VEN.NS=-sin(tlat).*cos(tlon).*Vx' - sin(tlat).*sin(tlon).*Vy' + cos(tlat).*Vz';
VEN.UD= cos(tlat).*cos(tlon).*Vx' + cos(tlat).*sin(tlon).*Vy' + sin(tlat).*Vz';    % not need

end
%====================================================
function [SDT] = SDTvec(tvx,tvy,tvz)

% Calculate Strike, Dip, Tensile vectors of triangle
% Input-----------------------------   
%  tvx  : x-coordinate(X) of triangle vertices. 
%  tvy  : y-coordinate(Y) of triangle vertices. 
%  tvz  : z-coordinate(dep) of triangle vertices. 
% Output---------------------------
%  SDT  : xyz-coordinates(XY-plane) of strike, dip, tensile vectors of triangle.
% ---------------------------------

% Strike, Dip, Tensile vectors on xyz-coordinate
normVec        = cross([tvx(2);tvy(2);tvz(2)]-[tvx(1);tvy(1);tvz(1)], [tvx(3);tvy(3);tvz(3)]-[tvx(1);tvy(1);tvz(1)]);
normVec        = normVec./norm(normVec);
if (normVec(3) < 0) % Enforce clockwise circulation
   normVec     = -normVec;
   [tvx(2),tvx(3)] = swap(tvx(2), tvx(3));
   [tvy(2),tvy(3)] = swap(tvy(2), tvy(3));
   [tvz(2),tvz(3)] = swap(tvz(2), tvz(3));
end
strikeVec      = [-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
dipVec         = cross(normVec, strikeVec);

SDT=[strikeVec(:) dipVec(:) zeros(3,1)];

end
%====================================================
function [xyz]=enu2xyz(tvlon,tvlat,abs)

% degree to radian
tvlon=deg2rad(tvlon);
tvlat=deg2rad(tvlat);

% lonlat to xyz
tvx=cos(tvlat).*cos(tvlon).*abs;
tvy=cos(tvlat).*sin(tvlon).*abs;
tvz=sin(tvlat).*abs;

xyz.x=tvx;
xyz.y=tvy;
xyz.z=tvz;
% XYZ=[tvx(:) tvy(:) tvz(:)];

end
%====================================================
function [area,center]=calc_tri_area(s)

% Export area of triangles using function triangle_area.m
% --------------------------------------------
% input
% s         : data of triangle subfault
% 
% output
% center    : center coordinate of triangle (lom lat dep)
% area      : area of triangle
% --------------------------------------------

tg.lon=mean(s.lon(s.tri),2);
tg.lat=mean(s.lat(s.tri),2);
tg.dep=mean(s.dep(s.tri),2);
% center=[tg.lon(:) tg.lat(:) tg.dep(:)];
tri.x=s.x(s.tri);
tri.y=s.y(s.tri);
tri.z=s.z(s.tri);

area=zeros(length(s.tri),1);
for nn=1:length(s.tri)
    triver=[(tri.x(nn,:))' (tri.y(nn,:))' (tri.z(nn,:))'];
    area(nn)=triangle_area(triver);
end

% sort in area order <<---- option (to know min and max of area)
[~,sortindex]=sort(area);
area=area(sortindex);
center=[tg.lon(sortindex) tg.lat(sortindex) tg.dep(sortindex)];

% save(['area_tri_',num2str(nmesh)],'center','area');

end
%====================================================
function [area]=triangle_area(P,method)
% This function gives the area of a triangle
%
% [area]=triangle_area(Points, Method)
%
% Points: The Points should be a numeric array, of size 3xn, 
%         thus the points can be 2D, 3D... nD
% Method: Can be 'h' area calculation with Heron's formula 
%         or can be 'q' Orthogonal-triangular decomposition (default)
%
% Example: 
% P1=[0 0]; P2=[1 0.5]; P3=[0.5 1];
% area = triangle_area([P1;P2;P3])
%
% Version 1.1 updated on 2007-09-21 
% Added 'Orthogonal-triangular decomposition' after a usefull review of John D'Errico 

% Default output format
if(exist('method','var')==0), method='q'; end

% Check input
if((method~='h')&&(method~='q')), error('Unknown area calculation method'); end
[k,m]=size(P); if(k~=3), error('Points are not a 3xn array'); end

if(method=='h')
    % Length of edges
    L=[sqrt(sum((P(1,:)-P(2,:)).^2)) sqrt(sum((P(2,:)-P(3,:)).^2)) sqrt(sum((P(3,:)-P(1,:)).^2))];
    
    % Area calculation with Heron's formula
    s = ((L(1)+L(2)+L(3))/2); 
    area = sqrt(s*(s-L(1))*(s-L(2))*(s-L(3)));
else
    % Area calculation with Orthogonal-triangular decomposition
    [q,r] = qr((P(2:3,:) - repmat(P(1,:),2,1))');
    area=abs(prod(diag(r)))/2;
end
    
end
%====================================================
function [angles]=triangle_angles(P,format)
% This function gives the angles of a triangle
%
% This function splits a triangle in two triangles with a corner of 90
% degrees (Heron's formula), and than uses a sine to calculate te angles
%
% [angles]=triangle_angles(Points, Format)
%
% Points: The Points should be a numeric array, of size 3xn, 
%         thus the points can be 2D, 3D... nD
% Format: Can be 'd' degrees or 'r' radians (default)
% 
%
% Example: 
% P1=[0 0]; P2=[0.88 0.5]; P3=[0 1];
% [angles]=triangle_angles([P1;P2;P3],'d')
%
% Version 1.3 updated on 2007-11-08

% Default output format
if(exist('format','var')==0), format='r'; end

% Check input
if((format~='r')&&(format~='d')), error('Unknown Output Format'); end
[k,m]=size(P); if(k~=3), error('Points are not a 3xn array'); end

% Length of edges
L=[sqrt(sum((P(1,:)-P(2,:)).^2)) sqrt(sum((P(2,:)-P(3,:)).^2)) sqrt(sum((P(3,:)-P(1,:)).^2))];
        
% Edge length when split in two 90 degrees triangles
s = ((L(1)+L(2)+L(3))/2); 
h = (2/L(3))*sqrt(s*(s-L(1))*(s-L(2))*(s-L(3)));
x = (L(1)^2-L(2)^2+L(3)^2)/(2*L(3));

% Angle calculations
if (format=='d')
    angles(1)=asind(h/L(1));
    if(x<0), angles(1)=180-angles(1); end
    angles(3)=asind(h/L(2));
    if(x>L(3)), angles(3)=180-angles(3); end
    angles(2)=180-angles(3)-angles(1);
else
    angles(1)=asin(h/L(1));
    if(x<0), angles(1)=pi-angles(1); end
    angles(3)=asin(h/L(2));
    if(x>L(3)), angles(3)=pi-angles(3); end
    angles(2)=pi-angles(3)-angles(1);
end

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