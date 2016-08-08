function [PM]=platemotion(tlon,tlat,origlon,origlat)

% Site velocity relative plate momtion
% Input  : Pole position(tlon, tlat)
%        : Angular velocity(angvel)
% Output : Velocity Azimuth(PM.azi)

% Parameter
amean=6369D3;                 % radius of the Earth
psy=origlon;                  % rotarion angle around Z axis
phy=90-origlat;               % rotarion angle around Y axis
theta=0;                      % rotarion angle around X axis
psy2=90;                      % rotarion angle around Z axis again
plon=161.65;plat=54.74;pomeg=-1.168;   % PH plate subducting rate

% Transformation from degree to radian
plon=deg2rad(plon);
plat=deg2rad(plat);
pomeg=deg2rad(pomeg);
slon=deg2rad(tlon);
slat=deg2rad(tlat);
psy=deg2rad(psy);psy2=deg2rad(psy2);phy=deg2rad(phy);

% Transformation from pole coord. to XYZ
[omea]=[cos(plat).*cos(plon) cos(plat).*sin(plon) sin(plat)].*pomeg;
[tga]=[cos(slat).*cos(slon) cos(slat).*sin(slon) sin(slat)].*amean;
% ome.x=omea(:,1);ome.y=omea(:,2);ome.z=omea(:,2);
% tg.x=tga(:,1);tg.y=tga(:,2);tg.z=tga(:,3);

% Rotation matrix
rotX =[1 0 0;0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
rotY =[cos(phy) 0 -sin(phy); 0 1 0; sin(phy) 0 cos(phy)];
rotZ =[cos(psy) sin(psy) 0; -sin(psy) cos(psy) 0; 0 0 1];
rotZ2=[cos(psy2) sin(psy2) 0; -sin(psy2) cos(psy2) 0; 0 0 1];

% Transration the coordinate to origin=(alon0,alat0)
[poleXYZ]=rotZ2*rotY*rotZ*[omea]';
[trigXYZ]=rotZ2*rotY*rotZ*[tga]';

% Plate motion
[tgsize,~]=size(tga);
for ii=1:tgsize
    plm(:,ii)=cross(poleXYZ,trigXYZ(:,ii)).*1.0D-6;
    plmnorm(ii)=norm(plm(:,ii));
%     plmnorm(ii)=sqrt(plm(1,ii)^2+plm(2,ii)^2);
    unitplm(:,ii)=plm(:,ii)./plmnorm(ii);
end

% [omea]=[omea]';[tga]=[tga]';
% for ii=1:tgsize
%     plm(:,ii)=cross(omea,tga(:,ii));
%     plmnorm(ii)=norm(plm(:,ii));
%     unitplm(:,ii)=plm(:,ii)./plmnorm(ii);
% end
% 
% unitplm=rotY*rotZ*unitplm;




PM.X=unitplm(1,:);PM.Y=unitplm(2,:);PM.Z=unitplm(3,:);



end