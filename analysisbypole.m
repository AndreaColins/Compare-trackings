d='./tr4s/tr4s/open coil/';
ff = dir([d '*.tr4']);
tic
nr_files=size(ff,1);
string_list=cell(nr_files,1);
for i=1:nr_files
    string_list{i}=strcat(d,ff(i).name);
end
[result,v,power]=trackcorr2(char(string_list),1);
%save(strcat(d,'autocovopen coil.mat'),'-struct','v')
%save(strcat(d,'powercarbon pole.mat'),'-struct','power')
%save(strcat(d,'resultcarbon pole.mat'),'result')

toc