% function out_coupling_correlation(DIR,TCHA,BLK)
function out_coupling_correlation(DIR)
load(fullfile(DIR,'TCHA.mat'))
load(fullfile(DIR,'BLK.mat'))
NN=1;
folder=[DIR,'/coupling/correlation'];
if exist(folder,'dir')~=7; mkdir(folder); end
formatstr='# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s';
formatvar='%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f';
for NB=1:BLK(1).NB
  formatstr=[formatstr,' %6i'];
  formatvar=[formatvar,' %6.3f'];
  if NB==BLK(1).NB
    formatstr=[formatstr,'\n'];
    formatvar=[formatvar,'\n'];
  end
end
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
    if NF~=0
      FIDmain = fopen([folder,'/co_',num2str(NB1),'_',num2str(NB2),'.txt'],'w');
      FLTNUM = NN:NN+NF-1;
      cpcorr = TCHA.CORFLT(FLTNUM,:);
      clon = mean(BLK(1).BOUND(NB1,NB2).blon,2);
      clat = mean(BLK(1).BOUND(NB1,NB2).blat,2);
      cdep = mean(BLK(1).BOUND(NB1,NB2).bdep,2);
      outdata = [FLTNUM' ...
          BLK(1).BOUND(NB1,NB2).blon ...
          BLK(1).BOUND(NB1,NB2).blat ...
          BLK(1).BOUND(NB1,NB2).bdep ...
          clon clat cdep ...
          cpcorr];
      fprintf(FIDmain,formatstr,...
                      'Tri_No','Lon1','Lon2','Lon3','Lat1','Lat2','Lat3','Dep1','Dep2','Dep3',...
                      'C_Lon','C_Lat','C_Dep',(1:BLK(1).NB)');
      fprintf(FIDmain,formatvar,outdata');
      fclose(FIDmain);
      NN=NN+NF;
    end
  end
end

end