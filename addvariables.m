condition='Carbon Pole';
d=strcat('./videosselectedtextures/',condition,'/horizontal/');
d2='C:\Users\mqbpwac5\Desktop\videosforcontact2\';
ff = dir([d '*.tr4']);
tic
nr_files=size(ff,1);
string_list=cell(nr_files,1);

filename=strcat('C:\Users\mqbpwac5\Documents\MATLAB\videosselectedtextures\',condition,'\horizontal\result',condition,'.mat');
v=load(filename,'-mat');
result=v.result;

for i=1:nr_files
    string_list{i}=strcat(d,ff(i).name);
    
fname=char(string_list{i})
tmp=load(fname,'-mat');
tmp.whisker.theta_all=result(:,1,i);
tmp.whisker.kappa_all=result(:,4,i);
save(strcat(d2,fname(end-22:end)),'-struct','tmp')
clear tmp
    
end




