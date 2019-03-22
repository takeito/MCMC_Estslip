%% output triangle vertices
function write4gmt(dir,meshN)

% read file-------------------------------------
infile1=[dir,'/plate_phs_initial.mat']    % initial tringles
infile2=[dir,'/plate_phs',num2str(meshN),'.mat']        % optimized triangles
% export name of file---------------------------
exfile1=[dir,'/triangle_initial.txt']
exfile2=[dir,'/triangle_',num2str(meshN),'.txt']
% -------------------------------------------

mat2gmt(infile1,exfile1);
mat2gmt(infile2,exfile2);

end

%% read .mat file and export triangle vertices for GMT format
function mat2gmt(file_in,file_out)

load(file_in)
[r_index,~]=size(s.tri);
indexs=0;
fidt=fopen(file_out,'wt');
for ii=1:r_index
    fprintf(fidt,'%s\n','>');
    indexs=indexs+1;indexe=indexs+2;
    lon(1:3)=s.lon(s.tri(ii,:));
    lon(4)=s.lon(s.tri(ii,1));
    lat(1:3)=s.lat(s.tri(ii,:));
    lat(4)=s.lat(s.tri(ii,1));
    dep(1:3)=s.dep(s.tri(ii,:));
    dep(4)=s.dep(s.tri(ii,1));
    indexs=indexe+1;
    for mm=1:4
        fprintf(fidt,'%f %f %f\n',lon(mm),lat(mm),dep(mm));
    end
end

fclose(fidt);

% ****** original (use fortran program output4gmt.f) ********
% for ii=1:r_index
%     indexs=indexs+1;indexe=indexs+2;
%     lon(indexs:1:indexe,1)=s.lon(s.tri(ii,:));
%     lon(indexe+1,1)=s.lon(s.tri(ii,1));
%     lat(indexs:1:indexe,1)=s.lat(s.tri(ii,:));
%     lat(indexe+1,1)=s.lat(s.tri(ii,1));
%     dep(indexs:1:indexe,1)=s.dep(s.tri(ii,:));
%     dep(indexe+1,1)=s.dep(s.tri(ii,1));
%     indexs=indexe+1;
% end
% triangle=[lon(:) lat(:) dep(:)]';
% fid=fopen(file_out,'wt');
% fprintf(fid,'%f %f %f\n',triangle);
% fclose(fid);
% ***************************************************

end