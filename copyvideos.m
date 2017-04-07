function copyvideos
xlsfile='Selectedvideos.xlsx';
%condition=char('smooth pole','Open coil','Closed coil','Bamboo','Black Sandpaper','Carbon Pole','Cardboard','Toothpick','Wood');
[xls_info,txt] = xlsread(xlsfile);

filesh=char(txt(2:end,1));
filesv=char(txt(2:end,2));
condition=char(txt(2:end,3));
for i=91:size(filesh,1)
    
    %%%%%%%%%%horizontal

dir1='Z:\Pole\Video data\';
date=strcat('2016_',filesh(i,[9,10]),'_',filesh(i,[11,12]),'\TT3\');
i
file=filesh(i,1:end-1)
videoh=strcat(dir1,date,file)
destinationh=strcat('.\videosselectedtextures\',condition(i,:),'\horizontal\')
copyfile(videoh,destinationh)


%%%%%vertical
date=strcat('2016_',filesh(i,[9,10]),'_',filesh(i,[11,12]),'_cam2\TT3\');
file=filesv(i,1:end-1);
videov=strcat(dir1,date,file)
destinationv=strcat('.\videosselectedtextures\',condition(i,:),'\vertical\')
copyfile(videov,destinationv)

end

end