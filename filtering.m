function filtering
tic
condition=char('air','Smooth pole', 'Carbon Pole', 'Black Sandpaper', 'Closed coil','Open coil','Cardboard','Bamboo','Toothpick','Wood');
ncondition=size(condition,1);
for con=1:ncondition

d=strcat('./videosselectedtextures/processed/FrameCorrected/result',condition(con,:),'.mat')
d2='./videosselectedtextures/processed/bandpass/from20/';
bandpass=[1,20];
% ff = dir([d '*.tr4']);
% tic
% nr_files=size(ff,1);
% string_list=cell(nr_files,1);
% for i=1:nr_files
%     string_list{i}=strcat(d,ff(i).name);
% end
% char(string_list)
% [idx,result]=comparetr4(char(string_list),0);
% %result=result(200:end,:,:);

v=load(d,'-mat');
    result=v.result;
    idx=[1:1:size(result,1)];
    
filteredresult=zeros(size(result));
for i=1:size(result,3)
    %meanr=mean(result);
    azimuth=timeseries(result(:,1,i),[0:1:numel(result(:,1,i))-1]./1000);
    elevation=timeseries(result(:,2,i),[0:1:numel(result(:,2,i))-1]./1000);
    kcoronal=timeseries(result(:,3,i),[0:1:numel(result(:,3,i))-1]./1000);
    khorizontal=timeseries(result(:,4,i),[0:1:numel(result(:,4,i))-1]./1000);
    
    filteraz=idealfilter(azimuth,bandpass,'pass');
    filterel=idealfilter(elevation,bandpass,'pass');
    filterkc=idealfilter(kcoronal,bandpass,'pass');
    filterkh=idealfilter(khorizontal,bandpass,'pass');
    
    filteredresult(:,1,i)=filteraz.Data;
    filteredresult(:,2,i)=filterel.Data;
    filteredresult(:,3,i)=filterkc.Data;
    filteredresult(:,4,i)=filterkh.Data;
     plot(1:1:3488,[filteredresult(:,1,i),result(:,1,i)])
     pause(2)

end
clear result
result=filteredresult;
size(result)
save(strcat(d2,strcat('result',condition(con,:),'1to20.mat')),'result')
toc
end
end