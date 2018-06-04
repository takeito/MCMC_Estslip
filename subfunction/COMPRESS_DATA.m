%% Compress CHA sampling
function COMPRESS_DATA(CHA,PRM,ITR,NACC)
% Compressing CHA sampled parameter to int8
% 
CHA.Mc=single(CHA.Mc);
CHA.Mp=single(CHA.Mp);
CHA.Mi=single(CHA.Mi);
% if PRM.GPU==99&&gpuDeviceCount==0
if PRM.GPU==99
  MEANMc=mean(CHA.Mc,2);
  MEANMp=mean(CHA.Mp,2);
  MEANMi=mean(CHA.Mi,2);
  COVMc=cov(CHA.Mc');
  COVMp=cov(CHA.Mp');
  COVMi=cov(CHA.Mi');
else
  gCHA.Mc=gpuArray(CHA.Mc);
  gCHA.Mp=gpuArray(CHA.Mp);
  gCHA.Mi=gpuArray(CHA.Mi);
  MEANMc=mean(gCHA.Mc,2);
  MEANMp=mean(gCHA.Mp,2);
  MEANMi=mean(gCHA.Mi,2);
  COVMc=cov(gCHA.Mc');
  COVMp=cov(gCHA.Mp');
  COVMi=cov(gCHA.Mi');
  MEANMc=gather(MEANMc);
  MEANMp=gather(MEANMp);
  MEANMi=gather(MEANMi);
  COVMc=gather(COVMc);
  COVMp=gather(COVMp);
  COVMi=gather(COVMi);
end
% 
McMAX=max(CHA.Mc,[],2);
McMIN=min(CHA.Mc,[],2);
MpMAX=max(CHA.Mp,[],2);
MpMIN=min(CHA.Mp,[],2);
MiMAX=max(CHA.Mi,[],2);
MiMIN=min(CHA.Mi,[],2);
% 
Mcscale=100./(McMAX-McMIN);
McBASE=bsxfun(@times,bsxfun(@minus,CHA.Mc,McMIN),Mcscale.*2.55-128);
Mcint8=int8(McBASE);
Mpscale=100./(MpMAX-MpMIN);
MpBASE=bsxfun(@times,bsxfun(@minus,CHA.Mp,MpMIN),Mpscale.*2.55-128);
Mpint8=int8(MpBASE);
Miscale=100./(MiMAX-MiMIN);
MiBASE=bsxfun(@times,bsxfun(@minus,CHA.Mi,MiMIN),Miscale.*2.55-128);
Miint8=int8(MiBASE);
% 
binedge=int8(-128:127);
% 
for ii=1:size(Mcint8,1)
  cha.McCOMPRESS.NFLT(ii).Mcscale=Mcscale(ii);
  cha.McCOMPRESS.NFLT(ii).McMAX=McMAX(ii);
  cha.McCOMPRESS.NFLT(ii).McMIN=McMIN(ii);
  cha.McCOMPRESS.NFLT(ii).McHIST=histcounts(Mcint8(ii,:),binedge);
end
cha.McCOMPRESS.COVMc=COVMc;
cha.McCOMPRESS.MEANMc=MEANMc;
cha.McCOMPRESS.SMPMc=int8(McBASE);
% 
for ii=1:size(Mpint8,1)
  cha.MpCOMPRESS.NPOL(ii).Mpscale=Mpscale(ii);
  cha.MpCOMPRESS.NPOL(ii).MpMAX=MpMAX(ii);
  cha.MpCOMPRESS.NPOL(ii).MpMIN=MpMIN(ii);
  cha.MpCOMPRESS.NPOL(ii).MpHIST=histcounts(Mpint8(ii,:),binedge);
end
cha.MpCOMPRESS.COVMp=COVMp;
cha.MpCOMPRESS.MEANMp=MEANMp;
cha.MpCOMPRESS.SMPMp=int8(MpBASE);
% 
for ii=1:size(Miint8,1)
  cha.MiCOMPRESS.NINE(ii).Miscale=Miscale(ii);
  cha.MiCOMPRESS.NINE(ii).MiMAX=MiMAX(ii);
  cha.MiCOMPRESS.NINE(ii).MiMIN=MiMIN(ii);
  cha.MiCOMPRESS.NINE(ii).MiHIST=histcounts(Miint8(ii,:),binedge);
end
cha.MiCOMPRESS.COVMi=COVMi;
cha.MiCOMPRESS.MEANMi=MEANMi;
cha.MiCOMPRESS.SMPMi=int8(MiBASE);
% 
cha.AJR=NACC./PRM.CHA;
save(fullfile(PRM.DirResult,['CHA_test',num2str(ITR,'%03i')]),'cha','-v7.3');
% save(['./Result/CHA_test',num2str(ITR,'%03i'),'.mat'],'cha','-v7.3'); % test
% 
end
