function out_slip_coupling_correlation(DIR,TCHA2,BLK,TRI)

% corsmp=SlipCC.corsmp;

varsmp=diag(TCHA2.covsmp);
varsmp(varsmp<0)=0;
stdsmp=sqrt(varsmp);
corsmp=TCHA2.covsmp./(stdsmp*stdsmp');

folder=fullfile(DIR,'coupling/slipcorrelation');
if exist(folder,'dir')~=7; mkdir(folder); end
formatstr='# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s\n';
formatvar='%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n';
NN=1;
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(TRI(1).BOUND(NB1,NB2).clon,2);
    if NF~=0
      FIDmain = fopen([folder,'/slipcc_',num2str(NB1),'_',num2str(NB2),'.txt'],'w');
      FLTNUM = NN:NN+NF-1;
      cpslipcorr = diag(corsmp(FLTNUM,FLTNUM+BLK(1).NB));
      clon = mean(BLK(1).BOUND(NB1,NB2).blon,2);
      clat = mean(BLK(1).BOUND(NB1,NB2).blat,2);
      cdep = mean(BLK(1).BOUND(NB1,NB2).bdep,2);
      outdata = [FLTNUM' ...
          BLK(1).BOUND(NB1,NB2).blon ...
          BLK(1).BOUND(NB1,NB2).blat ...
          BLK(1).BOUND(NB1,NB2).bdep ...
          clon clat cdep ...
          cpslipcorr];
      fprintf(FIDmain,formatstr,...
                      'Tri_No','Lon1','Lon2','Lon3','Lat1','Lat2','Lat3','Dep1','Dep2','Dep3',...
                      'C_Lon','C_Lat','C_Dep','SlipCC');
      fprintf(FIDmain,formatvar,outdata');
      fclose(FIDmain);
      NN=NN+NF;
    end
  end
end

end