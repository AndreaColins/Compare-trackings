d='./tr4s/tr4s/closed coil/';
ff = dir([d '*.tr4']);
tic
nr_files=size(ff,1);
string_list=cell(nr_files,1);
for i=1:nr_files
    string_list{i}=strcat(d,ff(i).name);
end
v=trackcorr2(char(string_list),0)
save(strcat(d,'autocovclosed coil.mat'),'-struct','v')
toc