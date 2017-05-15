function compareautocov
conditions=char('Smooth pole','Closed coil', 'Open coil','Black sandpaper','Carbon pole','Toothpick','Bamboo','Wood','Cardboard','air');
npole=size(conditions,1);
colors=distinguishable_colors(npole);
for j=1:size(conditions,1)
    condition=conditions(j,:)

d2='./videosselectedtextures/processed/FrameCorrected/';

% autocov=load(strcat(d2,strcat('autocovacc',condition,'.mat')));
% azimuth(:,:,j)=autocov.autoazimuth;
% elevation(:,:,j)=autocov.autoelevation;
% kcoronal(:,:,j)=autocov.autokcoronal;
% khorizontal(:,:,j)=autocov.autokhorizontal;
autocov=load(strcat(d2,strcat('poweracc',condition,'.mat')));
freq=autocov.freq;
azimuth(:,:,j)=autocov.azimuth;
elevation(:,:,j)=autocov.elevation;
kcoronal(:,:,j)=autocov.kcoronal;
khorizontal(:,:,j)=autocov.khorizontal;
end
mid=floor(size(azimuth,1)/2+1);
figure
subplot(2,2,1)
for j=1:size(conditions,1)
    
plot(mean(azimuth(:,:,j),2),'Color',colors(j,:,:))
hold on
end
title('Azimuth acceleration autocovariance')
xlabel('Time [ms]')
hold off

subplot(2,2,2)
for j=1:size(conditions,1)
    
plot(mean(elevation(:,:,j),2),'Color',colors(j,:,:))
hold on
end
title('Elevation acceleration autocovariance')
xlabel('Time [ms]')
hold off

subplot(2,2,3)
for j=1:size(conditions,1)
    
plot(mean(kcoronal(:,:,j),2),'Color',colors(j,:,:))
hold on
end
title('K coronal acceleration autocovariance')
xlabel('Time [ms]')
hold off
subplot(2,2,4)
for j=1:size(conditions,1)
    
plot(mean(khorizontal(:,:,j),2),'Color',colors(j,:,:))
hold on
end
title('K horizontal acceleration autocovariance')
xlabel('Time [ms]')
hold off
end