function PoleHistgram(folder,burn)
load(fullfile(folder,'TCHA.mat'));
thin = 10;

meanp.x = TCHA.AVEPOL(1:3:end);
meanp.y = TCHA.AVEPOL(2:3:end);
meanp.z = TCHA.AVEPOL(3:3:end);
[meanp.lat,meanp.lon,meanp.ang] = xyzp2lla(meanp.x,meanp.y,meanp.z);


burnid = round(size(TCHA.SMPPOL,2)*burn/100);
p.x = TCHA.SMPPOL(1:3:end,burnid:thin:end);
p.y = TCHA.SMPPOL(2:3:end,burnid:thin:end);
p.z = TCHA.SMPPOL(3:3:end,burnid:thin:end);
nsamp = size(p.x,2);
[p.lat,p.lon,p.ang] = xyzp2lla(p.x,p.y,p.z);
sor.lat = sort(p.lat);
sor.lon = sort(p.lon);
sor.ang = sort(p.ang);
cllo = round(0.5*nsamp*0.05);
clup = round(0.5*nsamp*1.95);
x95.lat = sor.lat(:,clup) - sor.lat(:,cllo);
x95.lon = sor.lon(:,clup) - sor.lon(:,cllo);
x95.ang = sor.ang(:,clup) - sor.ang(:,cllo);

% Coast line
filename   = gunzip('gshhs_c.b.gz', tempdir);
shorelines = gshhs(filename{1});
figure(100); clf(100)
geoshow([shorelines.Lat],[shorelines.Lon])
% for ns = 1:size(shorelines,1)
% plot(shorelines(ns).Lon,shorelines(ns).Lat,'k','LineWidth',1); hold on
% end

for nblk = 1:size(x95.lat,1)
  %   wbin = max(x95.lon(nblk),x95.lat(nblk));
  wbinlon = x95.lon(nblk);
  wbinlat = x95.lat(nblk);
  if wbinlon<=0; wbinlon = 1; end
  if wbinlat<=0; wbinlat = 1; end
  wbin = [wbinlon*1e-3, wbinlat*1e-3];
  [N,xedge,yedge] = histcounts2(p.lon(nblk,:),p.lat(nblk,:),'BinWidth',wbin,'Normalization','probability');
  x = (xedge(1:end-1) + xedge(2:end)) ./ 2;
  y = (yedge(1:end-1) + yedge(2:end)) ./ 2;
  [X,Y] = meshgrid(x,y);
  
  figure(100); hold on
  %   plot(p.lon(nblk,:),p.lat(nblk,:),'r.','MarkerSize',0.01); hold on
  histogram2(p.lon(nblk,:),p.lat(nblk,:),'BinWidth',wbin,'Normalization','probability','DisplayStyle','tile'); hold on
  %   plot(meanp.lon(nblk),meanp.lat(nblk),'kx'); hold on
  text(double(meanp.lon(nblk)),double(meanp.lat(nblk)),num2str(nblk)); hold on
  %   if ~(size(x,2)==1 || size(y,2)==1)
  %     contour(X,Y,N','LineWidth',1.5)
  %   end
  
  figure(nblk)
  %   plot(p.lon(nblk,:),p.lat(nblk,:),'r.','MarkerSize',0.01); hold on
  histogram2(p.lon(nblk,:),p.lat(nblk,:),'BinWidth',wbin,'Normalization','probability','DisplayStyle','tile'); hold on
  %   plot(meanp.lon(nblk),meanp.lat(nblk),'kx'); hold on
  text(double(meanp.lon(nblk)),double(meanp.lat(nblk)),num2str(nblk)); hold on
  %   if ~(size(x,2)==1 || size(y,2)==1)
  %     contour(X,Y,N','LineWidth',1.5)
  %   end
  savefig(nblk,fullfile(folder,['pole_',num2str(nblk)]))
end
savefig(100,fullfile(folder,['pole_',num2str(nblk)]))

end

function [lat,lon,ang]=xyzp2lla(X,Y,Z)
% XYZP2LLA  Converts Shpear coordinates from cartesian. Vectorized.
% GRS80
% CODE BY T.ITO 2017/03/11     ver0.1
% lat: deg, lon: deg, ang: deg/m.y.
lat=atan2(Z,sqrt(X.*X+Y.*Y)).*180/pi;
lon=atan2(Y,X).*180/pi;
ang=sqrt(X.*X+Y.*Y+Z.*Z).*(1e6.*(180./pi));
end
