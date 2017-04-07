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
function poleposition2
%%%%%%%%%%%%selected
xlsfile='Selectedvideos.xlsx';
[xls_info,txt] = xlsread(xlsfile);
textdata=char(txt(2:end,1));
dir1='Z:\Pole\Video data\';
fileID=fopen('measurements2.txt','w');
 nvideos=size(textdata,1);
%%%%%%preselected
xlsfile2='preselectedvideos2.xlsx';
[xls_info2,txt2] = xlsread(xlsfile2);
texture=2;
day=1;
I=find((xls_info2(:,4)==texture)&(xls_info2(:,5)==day))
subplot(1,2,1)
scatter(xls_info2(I,7),470-xls_info2(I,8),'b');
hold on   
ifiles=[11:1:15];
for i=ifiles
    date=strcat('2016_',textdata(i,[9,10]),'_',textdata(i,[11,12]),'\TT3\');
    file=textdata(i,1:end-1);
    videoh=strcat(dir1,date,file);
    radius=xls_info(i,1);
    
[posvideo,posvideo2,~]=contourvideos(videoh,radius,1,200);
posvideo=posvideo';
posvideo2=posvideo2';
subplot(1,2,1)
scatter(xls_info(i,2),470-xls_info(i,3),'r')
hold on
%radial=sqrt(posvideo(:,1).^2+posvideo(:,2).^2);
%radialD(i,:)=mean(radial);
%radialdes(i,:)=std(radial);
barpos(i,:)=mean(posvideo)';
despos(i,:)=std(posvideo)';
radiald(i,:)=mean(sqrt((posvideo(:,1)-posvideo2(:,1)).^2+(posvideo(:,2)-posvideo2(:,2)).^2));
radialstd(i,:)=std(sqrt((posvideo(:,1)-posvideo2(:,1)).^2+(posvideo(:,2)-posvideo2(:,2)).^2));
anteriord(i,:)=mean(posvideo2(:,1)-posvideo(:,1));
laterald(i,:)=mean(posvideo2(:,2)-posvideo(:,2));
fileID=fopen('measurements2.txt','a');
fprintf(fileID,'%f %f %f %f\r\n',radiald(i,:),radialstd(i,:),anteriord(i,:),laterald(i,:));
%pause(1)
end
ylabel('y position [px]')
xlabel('x position [px]')
title('Distribution of positions')
subplot(1,2,2)
scatter(xls_info2(I,7),470-xls_info2(I,8),'b');
hold on
scatter(xls_info(ifiles,2),470-xls_info(ifiles,3),'r');
title('Close up distribution of positions')
legend('Preselected','Selected')
ylabel('y position [px]')
xlabel('x position [px]')
fclose(fileID);
%barpos=barpos';
toc
end