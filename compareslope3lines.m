function compareslope3lines
conditions=char('Smooth pole','Closed coil', 'Open coil','Black sandpaper','Carbon pole','Toothpick','Bamboo','Wood','Cardboard','air');
colors=distinguishable_colors(size(conditions,1));

azimuth=[];
elevation=[];
kcoronal=[];
khorizontal=[];
group=[];
azimuth2=[];
elevation2=[];
kcoronal2=[];
khorizontal2=[];
azimuth3=[];
elevation3=[];
kcoronal3=[];
khorizontal3=[];
for i=1:size(conditions,1)
    if i==10
        [slope,slope2,slope3]=powerslopeexploratory2('air',colors(i,:));
    else
[slope,slope2,slope3]=powerslopeexploratory2(conditions(i,:),colors(i,:));
        
    end

    %%%%%grupo, de ahi hasta nfiles1+nfiles2 es otro, y asi hasta n
nfiles=size(slope,1);
group=[group;ones(nfiles,1)*i];
%%%%%%crear un vector con todos los datos de hw en la misma columna 
azimuth=[azimuth; slope(:,1)];
elevation=[elevation;slope(:,2)];
kcoronal=[kcoronal;slope(:,3)];
khorizontal=[khorizontal;slope(:,4)];

%%%%%%%%%%%%%%%%%%%%%%same but slope2
azimuth2=[azimuth2; slope2(:,1)];
elevation2=[elevation2;slope2(:,2)];
kcoronal2=[kcoronal2;slope2(:,3)];
khorizontal2=[khorizontal2;slope2(:,4)];

    
    %%%%%%%%%%%%%%%%%%%%%%same but slope2
azimuth3=[azimuth3; slope3(:,1)];
elevation3=[elevation3;slope3(:,2)];
kcoronal3=[kcoronal3;slope3(:,3)];
khorizontal3=[khorizontal3;slope3(:,4)];

    meanslope(i,:)=mean(slope);
    meanslope2(i,:)=mean(slope2);
    meanslope3(i,:)=mean(slope3);
    stdslope(i,:)=std(slope);
    stdslope2(i,:)=std(slope2);
    stdslope3(i,:)=std(slope3);
end

% d=strcat('./videosselectedtextures/processed/FrameCorrected/power','air','.mat');
%     v=load(d,'-mat');
%     logazimuth=log(v.azimuth);
%         plot(log(v.freq),logazimuth,'r')
%%%%%%%%%%%%%%%%%are slopes statistically different?
%%%%%%%%%%%%%%%%p<0.05 indicates they are different
[paz,tbl,stats] = anova1(azimuth,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigazimuth=m(sigidx,1:2);
sigazimuth(find(sigazimuth(:,1)~=1),:)=[];

[pel,tbl,stats] = anova1(elevation,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigelevation=m(sigidx,1:2);
sigelevation(find(sigelevation(:,1)~=1),:)=[];

[pkc,tbl,stats] = anova1(kcoronal,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigkcoronal=m(sigidx,1:2);
sigkcoronal(find(sigkcoronal(:,1)~=1),:)=[]

[pkh,tbl,stats] = anova1(khorizontal,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sighorizontal=m(sigidx,1:2);
sighorizontal(find(sighorizontal(:,1)~=1),:)=[]

%%%%%%%%%%%%%%%%%%slopes
%%%%%%%%%%%%%%%%%%%%%%%%
figure
%%%%%%%%%%%%%%%%%azimuth
subplot(1,2,1)
errorbar(meanslope(:,1),stdslope(:,1),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
xlab={'Air','Smooth','Closed', 'Open','Sandpaper','Carbon','Toothpick','Bamboo','Wood','Cardboard'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Azimuth (p=',num2str(paz,2),')'))
B=mat2cell(sigazimuth,ones(1,size(sigazimuth,1)),[2]);
sigstar(B,[]);
hold off
%%%%%%%%%%%%%%%elevation
subplot(1,2,2)
errorbar(meanslope(:,2),stdslope(:,2),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
%xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Elevation (p=',num2str(pel,2),')'))
B=mat2cell(sigelevation,ones(1,size(sigelevation,1)),[2]);
sigstar(B,[]);
hold off
% 
% %%%%%%%%%%%%kcoronal
% subplot(2,2,3)
% errorbar(meanslope(:,3),stdslope(:,3),'s')
% hold on
% %axis([0 size(conditions,1)+1 0 10])
% %xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon',''};
%  set(gca,...
%  'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
% title(strcat('Kappa Coronal (p=',num2str(pkc,2),')'))
% B=mat2cell(sigkcoronal,ones(1,size(sigkcoronal,1)),[2]);
% sigstar(B,[]);
% hold off
% 
% %%%%%%%%%%%%%%khorizontal
% subplot(2,2,4)
% errorbar(meanslope(:,4),stdslope(:,4),'s')
% hold on
% %axis([0 size(conditions,1)+1 0 10])
% %xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
%  set(gca,...
%  'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
% title(strcat('Kappa Horizontal (p=',num2str(pkh,2),')'))
% B=mat2cell(sighorizontal,ones(1,size(sighorizontal,1)),[2]);
% sigstar(B,[]);
% hold off



%%%%%%%%%%%%%%%%%%are slope2 statistical different?
%%%%%%%%%%%%%%%%%p<0.05 indicates they are different
[paz,tbl,stats] = anova1(azimuth2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigazimuth=m(sigidx,1:2);
sigazimuth(find(sigazimuth(:,1)~=1),:)=[];

[pel,tbl,stats] = anova1(elevation2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigelevation=m(sigidx,1:2);
sigelevation(find(sigelevation(:,1)~=1),:)=[];

[pkc,tbl,stats] = anova1(kcoronal2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigkcoronal=m(sigidx,1:2);
sigkcoronal(find(sigkcoronal(:,1)~=1),:)=[];

[pkh,tbl,stats] = anova1(khorizontal2,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sighorizontal=m(sigidx,1:2);
sighorizontal(find(sighorizontal(:,1)~=1),:)=[]
%%%%%%%%%%%%%%%%%%%%slope2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1)
errorbar(meanslope2(:,1),stdslope2(:,1),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
%xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon',''};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Azimuth (p=',num2str(paz,2),')'))
B=mat2cell(sigazimuth,ones(1,size(sigazimuth,1)),[2]);
sigstar(B,[]);

hold off

subplot(1,2,2)
errorbar(meanslope2(:,2),stdslope2(:,2),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
%xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Elevation (p=',num2str(pel,2),')'))
B=mat2cell(sigelevation,ones(1,size(sigelevation,1)),[2]);
sigstar(B,[]);
hold off

% subplot(2,2,3)
% errorbar(meanslope2(:,3),stdslope2(:,3),'s')
% hold on
% %axis([0 size(conditions,1)+1 0 10])
% %xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
%  set(gca,...
%  'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
% title(strcat('Kappa Coronal (p=',num2str(pkc,2),')'))
% B=mat2cell(sigkcoronal,ones(1,size(sigkcoronal,1)),[2]);
% sigstar(B,[]);
% hold off
% 
% subplot(2,2,4)
% errorbar(meanslope2(:,4),stdslope2(:,4),'s')
% hold on
% %axis([0 size(conditions,1)+1 0 10])
% %xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
%  set(gca,...
%  'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
% title(strcat('Kappa Horizontal (p=',num2str(pkh,2),')'))
% B=mat2cell(sighorizontal,ones(1,size(sighorizontal,1)),[2]);
% sigstar(B,[]);
% hold off


%%%%%%%%%%%%%%%%%%are slope3 statistical different?
%%%%%%%%%%%%%%%%%p<0.05 indicates they are different
[paz,tbl,stats] = anova1(azimuth3,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigazimuth=m(sigidx,1:2);
sigazimuth(find(sigazimuth(:,1)~=1),:)=[];

[pel,tbl,stats] = anova1(elevation3,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigelevation=m(sigidx,1:2);
sigelevation(find(sigelevation(:,1)~=1),:)=[];

[pkc,tbl,stats] = anova1(kcoronal3,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sigkcoronal=m(sigidx,1:2);
sigkcoronal(find(sigkcoronal(:,1)~=1),:)=[];

[pkh,tbl,stats] = anova1(khorizontal3,group,'off');

m=multcompare(stats,'display','off');
sigidx=find(m(:,end)<0.05);
sighorizontal=m(sigidx,1:2);
sighorizontal(find(sighorizontal(:,1)~=1),:)=[]

%%%%%%%%%%%%%%%%%%%%slope3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1)
errorbar(meanslope3(:,1),stdslope3(:,1),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
%xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Azimuth (p=',num2str(paz,2),')'))
B=mat2cell(sigazimuth,ones(1,size(sigazimuth,1)),[2]);
sigstar(B,[]);

hold off

subplot(1,2,2)
errorbar(meanslope3(:,2),stdslope3(:,2),'s')
hold on
%axis([0 size(conditions,1)+1 0 10])
%xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
 set(gca,...
 'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
title(strcat('Elevation (p=',num2str(pel,2),')'))
B=mat2cell(sigelevation,ones(1,size(sigelevation,1)),[2])
sigstar(B,[]);
hold off

% subplot(2,2,3)
% errorbar(meanslope3(:,3),stdslope3(:,3),'s')
% hold on
% %axis([0 size(conditions,1)+1 0 10])
% %xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
%  set(gca,...
%  'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
% title(strcat('Kappa Coronal (p=',num2str(pkc,2),')'))
% B=mat2cell(sigkcoronal,ones(1,size(sigkcoronal,1)),[2]);
% sigstar(B,[]);
% hold off
% 
% subplot(2,2,4)
% errorbar(meanslope3(:,4),stdslope3(:,4),'s')
% hold on
% %axis([0 size(conditions,1)+1 0 10])
% %xlab={'Air','Smooth pole','Closed', 'Open','Sandpaper','Carbon'};
%  set(gca,...
%  'XTickLabel',xlab,'XTick',[1:size(conditions,1)])
% title(strcat('Kappa Horizontal (p=',num2str(pkh,2),')'))
% B=mat2cell(sighorizontal,ones(1,size(sighorizontal,1)),[2]);
% sigstar(B,[]);
% hold off


clear all

end

