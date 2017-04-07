function goodtrials
%%%%%%define data
B=zeros(100,206);
condition=char('Smooth pole', 'Open coil','Closed coil', 'Bamboo','Black Sandpaper','Carbon Pole','Cardboard','Toothpick','Wood');
radius=[14,22,24,15,30,3,15,18,15];
ncondition=size(condition,1);
fileidx=1;

%%%%get vertical video
[xls_info,txt] = xlsread('Selectedvideos.xlsx');
filesv=char(txt(2:end,2));
C=str2num(filesv(:,end-10:end-4));

for con=1:ncondition
d=strcat('./videosselectedtextures/',condition(con,:),'/horizontal/')
ff = dir([d '*.dat'])
nr_files=size(ff,1)
string_list=cell(nr_files,1);

for i=2:nr_files
    string_list{i}=strcat(d,ff(i).name);
    
fname=char(string_list{i});
 video.fid = fopen(fname,'r')
 video.header = read_mikrotron_datfile_header(video.fid);

B(fileidx,1)=str2num(fname(end-9:end-4));
B(fileidx,2)=fileidx;
B(fileidx,3)=video.header.startframe;
B(fileidx,4)=video.header.triggerframe;
B(fileidx,5)=video.header.nframes;
B(fileidx,7)=1;%firsttouch
B(fileidx,8)=5;%firsttouch
B(fileidx,9)=radius(con);%radius
B(fileidx,10)=C(fileidx);%vertical video
B(fileidx,204)=1;%trial type/Left by default/go
B(fileidx,205)=1;%response/Left by default/go
B(fileidx,206)=1;%response/Left by default/go
fileidx=fileidx+1;
fclose(video.fid);
end
end
% %%%%%%save in xlsfile
 filename = 'good_trials.xlsx';
 A = {'movie name',	'movie_id',	'start_frame',	'trigger',	'nframes','trigger_movie'};
 sheet = 1;
 xlRange = 'A1';
 xlswrite(filename,A,sheet,xlRange)
 xlRange = 'J2';
 xlswrite(filename,C,sheet,xlRange)
 xlRange = 'A2';
 xlswrite(filename,B,sheet,xlRange)
 close all
end