function compareslope
conditions=char('air','smooth pole','closed coil', 'open coil','black sandpaper','carbon pole');
for i=1:size(conditions,1)
    [slope,intercept]=powerslope2(conditions(i,:));
    meanslope(i,:)=mean(slope);
    meanintercept(i,:)=mean(intercept);
    stdslope(i,:)=std(slope);
    stdintercept(i,:)=std(intercept);
end
%%%%%%%%%%%%%%%%%%slopes
figure

subplot(2,2,1)
errorbar(meanslope(:,1),stdslope(:,1),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Azimuth')
hold off

subplot(2,2,2)
errorbar(meanslope(:,2),stdslope(:,2),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Elevation')
hold off

subplot(2,2,3)
errorbar(meanslope(:,3),stdslope(:,3),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Kappa Coronal')
hold off

subplot(2,2,4)
errorbar(meanslope(:,4),stdslope(:,4),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Kappa Horizontal')
hold off

%%%%%%%%%%%%%%%%%%%%intercept
figure
subplot(2,2,1)
errorbar(meanintercept(:,1),stdintercept(:,1),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Azimuth')
hold off

subplot(2,2,2)
errorbar(meanintercept(:,2),stdintercept(:,2),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Elevation')
hold off

subplot(2,2,3)
errorbar(meanintercept(:,3),stdintercept(:,3),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Kappa Coronal')
hold off

subplot(2,2,4)
errorbar(meanintercept(:,4),stdintercept(:,4),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title('Kappa Horizontal')
hold off
end

