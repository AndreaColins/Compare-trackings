function comparefwhh
conditions=char('air','smooth pole','closed coil', 'open coil','black sandpaper','carbon pole');

ncon=size(conditions,1);
for i=1:size(conditions,1)
%%%%%%to compare fwhh with gaussing fit use
%hw=fwhh(conditions(i,:));
hw=fwhh(conditions(i,:));
meanhw(i,:)=mean(hw);
stdhw(i,:)=std(hw);
end
%%%%%%%%%%%%%%plot con 4 subplot, uno por cada variable.
subplot(2,2,1)
errorbar(meanhw(:,1),stdhw(:,1),'s')
hold on
axis([0 ncon+1 0 700])
xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon'};

 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title('Azimuth')
hold off


subplot(2,2,2)
errorbar(meanhw(:,2),stdhw(:,2),'s')
hold on
axis([0 ncon+1 0 700])
xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title('Elevation')

hold off
subplot(2,2,3)
errorbar(meanhw(:,3),stdhw(:,3),'s')
hold on
axis([0 ncon+1 0 700])
xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon'};

 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title('Kappa Coronal')

subplot(2,2,4)
errorbar(meanhw(:,4),stdhw(:,4),'s')
hold on
axis([0 ncon+1 0 400])
xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon'};

 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title('Kappa Horizontal')
hold off
hold off
end