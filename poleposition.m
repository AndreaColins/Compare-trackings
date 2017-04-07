% function [barpos,despos,radialD,radialdes,string_list]=poleposition
% 
% d='.Z:\Pole\Video data\2016_08_12\TT3\';
% firstv=0;
% lastv=inf;
% ff = dir([d '*.dat']);
% tic
% nr_files=size(ff,1);
% string_list=cell(nr_files,1);
% 
%     
% for i=1:nr_files
%     string_list{i}=strcat(d,ff(i).name);
%    %string_list{i}(11:12)
%     %char(string_list{i}(34:39))
% if strcmp(string_list{i}(21:22),'TT')&& str2num(string_list{i}(34:39))>=firstv && str2num(string_list{i}(34:39))<=lastv 
% char(string_list{i})
% [posvideo,~]=poletracker(char(string_list{i}),18,0,10);
% posvideo=posvideo';
% radial=sqrt(posvideo(:,1).^2+posvideo(:,2).^2);
% radialD(i,:)=mean(radial);
% radialdes(i,:)=std(radial);
% barpos(i,:)=mean(posvideo)';
% despos(i,:)=std(posvideo)';
% 
% %pause(1)
% end
% end
% 
% %barpos=barpos';
% toc
% end
function poleposition

d='Z:\Pole\Video data\2016_09_08\TT3\';
txt=importdata('videostomeasure.txt');
textdata=char(txt);
textdata=textdata(:,2:end-1);
nvideos=size(textdata,1);
radius=15;
fileID=fopen('measurements.txt','w');

%setGlobalID(fileID);
    
for i=1:nvideos
    
    strcat(d,textdata(i,:))
[posvideo,~]=poletracker(strcat(d,textdata(i,:)),radius,1,10);
posvideo=posvideo';
%radial=sqrt(posvideo(:,1).^2+posvideo(:,2).^2);
%radialD(i,:)=mean(radial);
%radialdes(i,:)=std(radial);
barpos(i,:)=mean(posvideo)'
despos(i,:)=std(posvideo)'
fileID=fopen('measurements.txt','a');
fprintf(fileID,'%s %f %f %f %f\r\n',textdata(i,:),barpos(i,:),despos(i,:));
%pause(1)
end

fclose(fileID);
%barpos=barpos';
toc
end