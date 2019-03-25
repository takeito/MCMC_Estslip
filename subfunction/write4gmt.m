%% output triangle vertices
function write4gmt(folder)

file = dir(fullfile(folder,'area*.mat'));
namesprt = strsplit(file.name,{'_','.mat'});
meshN = namesprt{3};
% read file-------------------------------------
infile1 = fullfile(folder,'plate_phs_initial.mat');    % initial tringles
infile2 = fullfile(folder,['plate_phs',num2str(meshN),'.mat']);        % optimized triangles
% export name of file---------------------------
exfile1 = fullfile(folder,'triangle_initial.txt');
exfile2 = fullfile(folder,['triangle_',num2str(meshN),'.txt']);
% -------------------------------------------

mat2gmt(infile1,exfile1);
mat2gmt(infile2,exfile2);

end

%% read .mat file and export triangle vertices for GMT format
function mat2gmt(file_in,file_out)

load(file_in)
[r_index,~]=size(s.tri);
fidt=fopen(file_out,'wt');
for ii=1:r_index
    fprintf(fidt,'%s\n','>');
    lon(1:3)=s.lon(s.tri(ii,:));
    lon(4)=s.lon(s.tri(ii,1));
    lat(1:3)=s.lat(s.tri(ii,:));
    lat(4)=s.lat(s.tri(ii,1));
    dep(1:3)=s.dep(s.tri(ii,:));
    dep(4)=s.dep(s.tri(ii,1));
    fprintf(fidt,'%f %f %f\n',[lon;lat;dep]);
end
fclose(fidt);

end