function comparefwhh
conditions=char('air','smooth pole','closed coil', 'open coil','black sandpaper','carbon pole','Bamboo','Toothpick');
azimuth=[];
elevation=[];
kcoronal=[];
khorizontal=[];
group=[];
ncon=size(conditions,1);
for i=1:size(conditions,1)
%%%%%%to compare fwhh with gaussing fit use
%hw=fwhh(conditions(i,:));
hw=fwhh(conditions(i,:));
%%%%%%crear un vector con todos los datos de hw en la misma columna 
azimuth=[azimuth; hw(:,1)];
elevation=[elevation;hw(:,2)];
kcoronal=[kcoronal;hw(:,3)];
khorizontal=[khorizontal;hw(:,4)];
%%%%%crear un vector group que indique que los primeros nfiles sonde un
%%%%%grupo, de ahi hasta nfiles1+nfiles2 es otro, y asi hasta n
nfiles=size(hw,1);
group=[group;ones(nfiles,1)*i];
meanhw(i,:)=mean(hw);
stdhw(i,:)=std(hw);
end


% [p,tbl,stats] = anova1(kcoronal,group,'off')
% multcompare(stats)
% [p,tbl,stats] = anova1(khorizontal,group,'off')
% multcompare(stats)


%%%%%%%%%%%%%%plot con 4 subplot, uno por cada variable.

%%%%%%%%%%%%%%%%%%are they statistical different?
%%%%%%%%%%%%%%%%%p<0.05 indicates they are different
[paz,tbl,stats] = anova1(azimuth,group,'off');

m=multcompare(stats);
sigidx=find(m(:,end)<0.05);
sigazimuth=m(sigidx,1:2)


[pel,tbl,stats] = anova1(elevation,group,'off');

m=multcompare(stats);
sigidx=find(m(:,end)<0.05);
sigelevation=m(sigidx,1:2)


[pkc,tbl,stats] = anova1(kcoronal,group,'off');

m=multcompare(stats);
sigidx=find(m(:,end)<0.05);
sigkcoronal=m(sigidx,1:2)


[pkh,tbl,stats] = anova1(khorizontal,group,'off');

m=multcompare(stats);
sigidx=find(m(:,end)<0.05);
sighorizontal=m(sigidx,1:2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%azimuth
subplot(2,2,1)
errorbar(meanhw(:,1),stdhw(:,1),'s')
hold on
axis([0 ncon+1 0 400])
xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon','Bamboo','Toothpick'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title(strcat('Azimuth (p=',num2str(paz,2),')'))
B=mat2cell(sigazimuth,ones(1,size(sigazimuth,1)),[2]);
sigstar(B,[]);

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%elevation


subplot(2,2,2)
errorbar(meanhw(:,2),stdhw(:,2),'s')
hold on
axis([0 ncon+1 0 400])
%xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon',};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title(strcat('Elevation (p=',num2str(pel,2),')'))
B=mat2cell(sigelevation,ones(1,size(sigelevation,1)),[2]);
sigstar(B,[]);
hold off


%%%%%%%%%%%%%%%%%%%%%%%kcoronal
subplot(2,2,3)
errorbar(meanhw(:,3),stdhw(:,3),'s')
hold on
axis([0 ncon+1 0 400])
%xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon'};

 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title(strcat('Kappa Coronal (p=',num2str(pkc,2),')'))
B=mat2cell(sigkcoronal,ones(1,size(sigkcoronal,1)),[2]);
sigstar(B,[]);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%khorizontal
subplot(2,2,4)
errorbar(meanhw(:,4),stdhw(:,4),'s')
hold on
axis([0 ncon+1 0 400])

%xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon'};

 set(gca,...
 'XTickLabel',xlab,'XTick',[1:ncon])
title(strcat('Kappa Horizontal (p=',num2str(pkh,2),')'))
B=mat2cell(sighorizontal,ones(1,size(sighorizontal,1)),[2]);
sigstar(B,[]);
hold off

end