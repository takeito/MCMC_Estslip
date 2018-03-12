function show_histogram_C(DIR)
% -------------------------------------------------------------------------
% This program draws and saves 1-D histogram and time series sampling graph
% of coupling at each mesh.
% -------------------------------------------------------------------------
fprintf('Now loading %s ...',[DIR,'/TCHA.mat'])
load([DIR,'/TCHA.mat']);fprintf('load\n')
Mcedges=TCHA.Mcbin;
TSint=TCHA.Smpint;
% meshcode=0;
% prompt=['Enter the MESH number [1-',num2str(size(TCHA.AVEFLT,1)),']: '];
meshnum=[702 703 699 701 700 1190 1170 1149 1153 1162 785 786 695 696 787 789 697 698 1174 582 834 611 632 631 581 580 578 579 577 576 639 643 644 692 575 574 689 690 694,...
    927 929 931 1203 1204 1205 1018 1191 1193 1216 1145 1146 1132 1144 944 939 941 943 1119 1138 1218 846 1185 851 849 1184 1221,...
    1230 805 1088 991 986 728 730 987 720 901 1077 918 953 726 922 951 727 964 854 962 857 855 856 782 780 958 821 824 747 751 826 823 825 822 808 812 814 1208 1013 1101 1014 683 686,...
    711 741 739 740 1073 1187 1188 1186 1183 1189 685 616 617 623 618 1045 679 672 673 674 831 667 664 638 832 668 634 636 957 614 613 1210 1067 872 873 612 798 869 678 1211 799 667 1031 891 585 586 661 603 955 956 680 681 670 671 663 660 572 573 652 602 647 646 597]; %OGnew
% meshnum=1:size(TCHA.AVEFLT,1);
for ii=1:size(meshnum,2)
  meshcode=meshnum(ii);
% while meshcode>=0
%   meshcode=input(prompt);
  if meshcode>=1&&meshcode<=size(TCHA.AVEFLT,1)
%     Graph of Histogram
    fprintf('Drawing PDF of subfault %i ...',meshcode)
    fig1=figure('visible','off');
    IDmin=find(TCHA.HISTFLT(meshcode,:),1,'first');
    IDmax=find(TCHA.HISTFLT(meshcode,:),1,'last');
    if IDmin~=1;IDmin=IDmin-1;end
    if IDmax~=size(TCHA.HISTFLT,2); IDmax=IDmax+2;else;IDmax=IDmax+1;end
    h=histogram('BinEdges',Mcedges,'BinCounts',TCHA.HISTFLT(meshcode,:),'Normalization','probability');
    ax1=gca;
    wdirhisto=[DIR,'/histo'];
    EXID=exist(wdirhisto);
    if EXID~=7; mkdir(wdirhisto); end
    wfilehisto=fullfile(wdirhisto,['histo_coupling_mesh_',num2str(meshcode)]);
    savefig([wfilehisto,'.fig'])
    print(fig1,'-dpng',wfilehisto)
    ax1.XLim=[Mcedges(IDmin) Mcedges(IDmax)];
    xlabel('Coupling');
    ylabel('Probability density');
    ax1.FontSize=30;
    ax1.PlotBoxAspectRatio=[1 0.8 0.8];
    wfilehisto=fullfile(wdirhisto,['histo_coupling_mesh_',num2str(meshcode),'_zoom']);
    savefig([wfilehisto,'.fig'])
    print(fig1,'-dpng',wfilehisto)
    close(fig1);
%     Graph of time series
    fig2=figure('visible','off');
    TSx=1:TSint:TSint*size(TCHA.SMPFLT(meshcode,:),2);
    BNINID=TSx>0.01*TCHA.Burnin*TSx(end);
    plot(TSx,TCHA.SMPFLT(meshcode,:),'-r','LineWidth',1)
    hold on
    plot(TSx(BNINID),TCHA.SMPFLT(meshcode,BNINID),'-k','LineWidth',2)
    ax2=gca;
%     ax2.YLim=[-1,1];  % Coupling ratio raging -1 ~ 1
    ax2.YLim=[ 0,1];  % Coupling ratio raging  0 ~ 1
    ax2.FontSize=20;
    ax2.PlotBoxAspectRatio=[1 0.3 1];
    wdirts=[DIR,'/timeseries'];
    EXID=exist(wdirts);
    if EXID~=7; mkdir(wdirts); end
    wfilets=fullfile(wdirts,['timeseries_mesh_',num2str(meshcode)]);
    savefig([wfilets,'.fig'])
    print(fig2,'-dpng',wfilets)
    close(fig2);
    fprintf('Done. \n')
    continue
  elseif meshcode>size(TCHA.AVEFLT,1)
    fprintf('Enter the number smaller than %d\n',size(TCHA.AVEFLT,1))
  else
    return
  end
end  
end