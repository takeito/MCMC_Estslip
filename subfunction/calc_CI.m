% function [CI] = calc_CI(sortSMPFLT)
function calc_CI(DIR,burn)
load([DIR,'/TCHA.mat'])
sortSMPFLT=sort(TCHA.SMPFLT,2);
remid=size(sortSMPFLT,2)*burn/100;
tmp=remid;
remid=tmp+1;
sortSMPFLT=sortSMPFLT(:,remid:end);
% Calculate Confidence Interval (CI) for coupling estimates.

% 99% CI
l05=size(sortSMPFLT,2)*(0.5/100);
r05=size(sortSMPFLT,2)-l05-1;
l05=l05+1;
smpflt_l05=sortSMPFLT(:,l05);
smpflt_r05=sortSMPFLT(:,r05);
smpflt_int=smpflt_r05-smpflt_l05;
CI99=[smpflt_l05,smpflt_r05,smpflt_int];

% 95% CI (2sigma)
l25=size(sortSMPFLT,2)*(2.5/100);
r25=size(sortSMPFLT,2)-l25-1;
l25=l25+1;
smpflt_l25=sortSMPFLT(:,l25);
smpflt_r25=sortSMPFLT(:,r25);
smpflt_int=smpflt_r25-smpflt_l25;
CI95=[smpflt_l25,smpflt_r25,smpflt_int];

% 90% CI
l50=size(sortSMPFLT,2)*(5.0/100);
r50=size(sortSMPFLT,2)-l50-1;
l50=l50+1;
smpflt_l50=sortSMPFLT(:,l50);
smpflt_r50=sortSMPFLT(:,r50);
smpflt_int=smpflt_r50-smpflt_l50;
CI90=[smpflt_l50,smpflt_r50,smpflt_int];

% 70% CI
l15=size(sortSMPFLT,2)*(15.0/100);
r15=size(sortSMPFLT,2)-l15-1;
l15=l15+1;
smpflt_l15=sortSMPFLT(:,l15);
smpflt_r15=sortSMPFLT(:,r15);
smpflt_int=smpflt_r15-smpflt_l15;
CI70=[smpflt_l15,smpflt_r15,smpflt_int];

% 68% CI (1sigma)
l16=size(sortSMPFLT,2)*(16/100);
r16=size(sortSMPFLT,2)-l16-1;
l16=l16+1;
smpflt_l16=sortSMPFLT(:,l16);
smpflt_r16=sortSMPFLT(:,r16);
smpflt_int=smpflt_r16-smpflt_l16;
CI68=[smpflt_l16,smpflt_r16,smpflt_int];

CI(1).data=CI99;
CI(1).percentile=99;
CI(2).data=CI95;
CI(2).percentile=95;
CI(3).data=CI90;
CI(3).percentile=90;
CI(4).data=CI70;
CI(4).percentile=70;
CI(5).data=CI68;
CI(5).percentile=68;

save(fullfile(DIR,'/CI.mat'),'ci','-v7.3');

end