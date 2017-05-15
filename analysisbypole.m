%filterH=bandpass
% condition='air';
% d=strcat('./tr4s/tr4s/',condition,'/');
% d2='./tr4s/processed/bandpass/upto5/';
conditions=char('Smooth pole','Closed coil', 'Open coil','Black sandpaper','Carbon pole','Toothpick','Bamboo','Wood','Cardboard','air');
for j=1:size(conditions,1)
    condition=conditions(j,:)
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
  save(strcat(d2,strcat('autocovacc',condition,'.mat')),'-struct','v')
  save(strcat(d2,strcat('poweracc',condition,'.mat')),'-struct','power')
  save(strcat(d2,strcat('resultacc',condition,'.mat')),'result')
end
toc