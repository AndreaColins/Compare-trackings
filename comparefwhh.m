function comparefwhh
conditions=char('air','smooth pole','closed coil', 'open coil');
conditions(1,:)
size(conditions,1)
for i=1:size(conditions,1)

hw=fwhh(conditions(i,:));
meanhw(i,:)=mean(hw);
stdhw(i,:)=std(hw);
end
%%%%%%%%%%%%%%plot con 4 subplot, uno por cada variable.
subplot(2,2,1)
errorbar(1:4,meanhw(:,1),stdhw(:,1),'s')
hold on
axis([0 5 0 700])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:4])
title('Azimuth')
hold off


subplot(2,2,2)
errorbar(1:4,meanhw(:,2),stdhw(:,2),'s')
hold on
axis([0 5 0 700])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:4])
title('Kappa Coronal')

hold off
subplot(2,2,3)
errorbar(1:4,meanhw(:,3),stdhw(:,3),'s')
hold on
axis([0 5 0 700])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:4])
title('Elevation')

subplot(2,2,4)
errorbar(1:4,meanhw(:,4),stdhw(:,4),'s')
hold on
axis([0 5 0 700])
xlab={'Air','Smooth pole','Closed coil', 'Open Coil'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:4])
title('Kappa Horizontal')
hold off
hold off
end