%filterH=bandpass
% condition='air';
% d=strcat('./tr4s/tr4s/',condition,'/');
% d2='./tr4s/processed/bandpass/upto5/';
condition='Wood';
d=strcat('./videosselectedtextures/',condition,'/horizontal/');
d2='./videosselectedtextures/processed/FrameCorrected/';
ff = dir([d '*.tr4']);
tic
nr_files=size(ff,1);
string_list=cell(nr_files,1);
for i=1:nr_files
    string_list{i}=strcat(d,ff(i).name);
end
[result,v,power]=trackcorr2(char(string_list),0);
  save(strcat(d2,strcat('autocov',condition,'.mat')),'-struct','v')
  save(strcat(d2,strcat('power',condition,'.mat')),'-struct','power')
  save(strcat(d2,strcat('result',condition,'.mat')),'result')

toc