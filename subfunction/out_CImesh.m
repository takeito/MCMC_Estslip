% function out_CImesh(DIR,ci,BLK)
function out_CImesh(DIR)
% Export the CI for each mesh.
load(fullfile(DIR,'BLK.mat'))
load(fullfile(DIR,'CI.mat'))
pf=pwd;
folder=fullfile(pf,DIR,'coupling');
for nci=1:size(CI,2)
  subfolder=fullfile(folder,['CI',num2str(CI(nci).percentile)]);
  write_CImeshfile(subfolder,CI(nci),BLK);
end
end
%%
function write_CImeshfile(subfolder,ci,BLK)
if exist(subfolder,'dir') ~= 7
  mkdir(subfolder);
end
NN=1;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
    if NF~=0
      FIDmain = fopen([subfolder,'/CI_',num2str(NB1),'_',num2str(NB2),'.txt'],'w');
      FLTNUM = NN:NN+NF-1;
      cimin=ci.data(FLTNUM,1);
      cimax=ci.data(FLTNUM,2);
      ciint=ci.data(FLTNUM,3);
      clon = mean(BLK(1).BOUND(NB1,NB2).blon,2);
      clat = mean(BLK(1).BOUND(NB1,NB2).blat,2);
      cdep = mean(BLK(1).BOUND(NB1,NB2).bdep,2);
      outdata = [FLTNUM' ...
          BLK(1).BOUND(NB1,NB2).blon ...
          BLK(1).BOUND(NB1,NB2).blat ...
          BLK(1).BOUND(NB1,NB2).bdep ...
          clon clat cdep ...
          cimin cimax ciint];
      fprintf(FIDmain,'# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %10s %10s %10s\n',...
                      'Tri_No','Lon1','Lon2','Lon3','Lat1','Lat2','Lat3','Dep1','Dep2','Dep3',...
                      'C_Lon','C_Lat','C_Dep','Cfmin','Cfmax','CfInt');
      fprintf(FIDmain,'%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.5f %10.5f %10.5f\n',outdata');
      fclose(FIDmain);
      NN=NN+NF;
    end
  end
end
end