function compareslope
conditions=char('air','smooth pole','closed coil', 'open coil','black sandpaper','carbon pole');

azimuth=[];
elevation=[];
kcoronal=[];
khorizontal=[];
group=[];
azimuth2=[];
elevation2=[];
kcoronal2=[];
khorizontal2=[];


for i=1:size(conditions,1)
    [slope,intercept]=powerslope2(conditions(i,:));
    
    %%%%%grupo, de ahi hasta nfiles1+nfiles2 es otro, y asi hasta n
nfiles=size(slope,1);
group=[group;ones(nfiles,1)*i];
%%%%%%crear un vector con todos los datos de hw en la misma columna 
azimuth=[azimuth; slope(:,1)];
elevation=[elevation;slope(:,2)];
kcoronal=[kcoronal;slope(:,3)];
khorizontal=[khorizontal;slope(:,4)];

%%%%%%%%%%%%%%%%%%%%%%same but intecept 
azimuth2=[azimuth2; intercept(:,1)];
elevation2=[elevation2;intercept(:,2)];
kcoronal2=[kcoronal2;intercept(:,3)];
khorizontal2=[khorizontal2;intercept(:,4)];

    meanslope(i,:)=mean(slope);
    meanintercept(i,:)=mean(intercept);
    stdslope(i,:)=std(slope);
    stdintercept(i,:)=std(intercept);
end


%%%%%%%%%%%%%%%%%%are slopes statistical different?
%%%%%%%%%%%%%%%%%p<0.05 indicates they are different
[paz,tbl,stats] = anova1(azimuth,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigazimuth=m(sigidx,1:2)

[pel,tbl,stats] = anova1(elevation,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigelevation=m(sigidx,1:2)


[pkc,tbl,stats] = anova1(kcoronal,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigkcoronal=m(sigidx,1:2)


[pkh,tbl,stats] = anova1(khorizontal,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sighorizontal=m(sigidx,1:2)


%%%%%%%%%%%%%%%%%%slopes
%%%%%%%%%%%%%%%%%%%%%%%%
figure
%%%%%%%%%%%%%%%%%azimuth
subplot(2,2,1)
errorbar(meanslope(:,1),stdslope(:,1),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Azimuth (p=',num2str(paz,2),')'))
B=mat2cell(sigazimuth,ones(1,size(sigazimuth,1)),[2]);
sigstar(B,[]);
hold off
%%%%%%%%%%%%%%%elevation
subplot(2,2,2)
errorbar(meanslope(:,2),stdslope(:,2),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Elevation (p=',num2str(pel,2),')'))
B=mat2cell(sigelevation,ones(1,size(sigelevation,1)),[2]);
sigstar(B,[]);
hold off

%%%%%%%%%%%%kcoronal
subplot(2,2,3)
errorbar(meanslope(:,3),stdslope(:,3),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Kappa Coronal (p=',num2str(pkc,2),')'))
B=mat2cell(sigkcoronal,ones(1,size(sigkcoronal,1)),[2]);
sigstar(B,[]);
hold off

%%%%%%%%%%%%%%khorizontal
subplot(2,2,4)
errorbar(meanslope(:,4),stdslope(:,4),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Kappa Horizontal (p=',num2str(pkh,2),')'))
B=mat2cell(sighorizontal,ones(1,size(sighorizontal,1)),[2]);
sigstar(B,[]);
hold off



%%%%%%%%%%%%%%%%%%are intercepts statistical different?
%%%%%%%%%%%%%%%%%p<0.05 indicates they are different
[paz,tbl,stats] = anova1(azimuth2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigazimuth=m(sigidx,1:2)

[pel,tbl,stats] = anova1(elevation2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigelevation=m(sigidx,1:2)


[pkc,tbl,stats] = anova1(kcoronal2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigkcoronal=m(sigidx,1:2)


[pkh,tbl,stats] = anova1(khorizontal2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sighorizontal=m(sigidx,1:2)

%%%%%%%%%%%%%%%%%%%%intercept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,2,1)
errorbar(meanintercept(:,1),stdintercept(:,1),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Azimuth (p=',num2str(paz,2),')'))
B=mat2cell(sigazimuth,ones(1,size(sigazimuth,1)),[2]);
sigstar(B,[]);

hold off

subplot(2,2,2)
errorbar(meanintercept(:,2),stdintercept(:,2),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Elevation (p=',num2str(pel,2),')'))
B=mat2cell(sigelevation,ones(1,size(sigelevation,1)),[2])
sigstar(B,[]);
hold off

subplot(2,2,3)
errorbar(meanintercept(:,3),stdintercept(:,3),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Kappa Coronal (p=',num2str(pkc,2),')'))
B=mat2cell(sigkcoronal,ones(1,size(sigkcoronal,1)),[2]);
sigstar(B,[]);
hold off

subplot(2,2,4)
errorbar(meanintercept(:,4),stdintercept(:,4),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Kappa Horizontal (p=',num2str(pkh,2),')'))
B=mat2cell(sighorizontal,ones(1,size(sighorizontal,1)),[2]);
sigstar(B,[]);
hold off

clear all

end

