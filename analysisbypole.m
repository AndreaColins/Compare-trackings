d='./tr4s/tr4s/air/';
condition='black sandpaper';
ff = dir([d '*.tr4']);
tic
nr_files=size(ff,1);
string_list=cell(nr_files,1);
for i=1:nr_files
    string_list{i}=strcat(d,ff(i).name);
end
[result,v,power]=trackcorr2(char(string_list),1);
%save(strcat(d,strcat('autocov',condition,'.mat')),'-struct','v')
%save(strcat(d,strcat('power',condition,'.mat')),'-struct','power')
%save(strcat(d,strcat('result',condition,'.mat')),'result')

toc