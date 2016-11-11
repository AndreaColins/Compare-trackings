d='./tr4s/tr4s/smooth pole/';
ff = dir([d '*.tr4']);
tic
nr_files=size(ff,1);
string_list=cell(nr_files,1);
for i=1:nr_files
    string_list{i}=strcat(d,ff(i).name);
end
timelag=trackcorr(char(string_list))
toc