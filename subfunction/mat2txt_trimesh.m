function mat2txt_trimesh(folder)
mat2txt_area(folder)
% mat2txt_sU(folder)
end
%%
function mat2txt_area(folder)

EXT='area*.mat';
file=dir([folder,'/',EXT]);
[NMesh,~]=size(file);
initial=[file(1).folder,'/','area_tri_initial.mat'];
for NB=1:NMesh
  if strcmp(fullfile(file(NB).folder,file(NB).name),initial)==1
    AREAINITIAL=file(NB);
  else
    AREA=file(NB);
  end
end

fname=[fullfile(AREA.folder,AREA.name)];
outareatxtl=[AREA.folder,'/area_land_HP.txt'];
outareatxts=[AREA.folder,'/area_sea_HP.txt'];

% Define the boundary of land and sea
File.land_sea='~/MasterResearch/inversion/MCMC_Estslip/land_sea_bound.txt';
Fid=fopen(File.land_sea);
bound=textscan(Fid,'%f%f%f');
fclose(Fid);
bound=cell2mat(bound);

% normalization of area
load(fname);
normtarea=tarea/max(tarea);

% separate fault to land and sea area
sID=inpolygon(tcenter(:,1),tcenter(:,2),bound(:,1),bound(:,2));
lID=ones(length(tarea),1);lID(sID)=0;
lID=logical(lID);
larea=normtarea(lID);
sarea=normtarea(sID);

Fid1=fopen(outareatxtl,'wt');
fprintf(Fid1,'%f\n',larea);
fclose(Fid1);

Fid2=fopen(outareatxts,'wt');
fprintf(Fid2,'%f\n',sarea);
fclose(Fid2);

disp(fname)
disp(['Mean area of onshore    = ',num2str(mean(larea)*max(tarea))])
disp(['Mean area of offshore   = ',num2str(mean(sarea)*max(tarea))])
disp(['Median area of onshore  = ',num2str(median(larea)*max(tarea))])
disp(['Median area of offshore = ',num2str(median(sarea)*max(tarea))])

end
%%
% function mat2txt_sU(arpha,folder)
function mat2txt_sU(folder)
    
% fname=['~/Desktop/initial_',num2str(arpha),'/Meshplate_phs_initial.mat'];
fname = fullfile(folder,'plate_phs_initial.mat');
outsUtxt = fullfile(folder,'U_initial.txt');
% outsUtxt=['~/Desktop/initial_',num2str(arpha),'/U',',num2str(arpha),','.txt'];

load(fname);

sortU=sort(s.U);
sortU=sortU/max(s.U);
fid1=fopen(outsUtxt,'wt');
fprintf(fid1,'%f\n',sortU);
fclose(fid1);

end