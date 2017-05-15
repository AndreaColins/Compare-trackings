function comparepower
tic
condition=char('Smooth pole', 'Carbon Pole', 'Black Sandpaper', 'Closed coil','Open coil','Cardboard','Bamboo','Toothpick','Wood','air');
ncondition=size(condition,1);
colors=distinguishable_colors(size(condition,1));

logazimuth=[];
varazimuth=[];
logelevation=[];
varelevation=[];
logkcor=[];
varkcor=[];
logkhor=[];
varkhor=[];
for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/poweracc',condition(con,:),'.mat');
    v=load(d,'-mat');
    logazimuth=v.azimuth;
    varazimuth=std(log(v.azimuth),0,2);
%     idx=find(v.freq>100,1,'First');
%     idx2=find(v.freq>300,1,'First');
%     power(con,:)=trapz(v.freq(idx:idx2),v.azimuth(idx:idx2,:,:),1);

    freq=log(v.freq);
    subplot(2,2,1)
    if con==ncondition
        plot(v.freq,mean(logazimuth,2),'Color',colors(con,:))
    else
        plot(v.freq,mean(logazimuth,2),'Color',colors(con,:))
    end
    hold on
   
end
title('Azimuth')
xlabel('Frequency [Hz]')
ylabel('Power frecuency')
hold off
legend(condition)

for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/poweracc',condition(con,:),'.mat');
    v=load(d,'-mat');
    logelevation=v.elevation;
    varelevation=std(log(v.elevation),0,2);
%     idx=find(v.freq>100,1,'First');
%     idx2=find(v.freq>300,1,'First');
%     power(con,:)=trapz(v.freq(idx:idx2),v.elevation(idx:idx2,:,:),1);


    subplot(2,2,2)
    if con==ncondition
        plot(v.freq,mean(logelevation,2),'Color',colors(con,:))
    else
        plot(v.freq,mean(logelevation,2),'Color',colors(con,:))
    end
    hold on
end
title('Elevation')
xlabel('Frequency [Hz]')
ylabel('Power frecuency')
hold off


for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/poweracc',condition(con,:),'.mat');
    v=load(d,'-mat');
    logkcoronal=v.kcoronal;
    varkcoronal=std(log(v.kcoronal),0,2);
%     idx=find(v.freq>100,1,'First');
%     idx2=find(v.freq>300,1,'First');
%     power(con,:)=trapz(v.freq(idx:idx2),v.kcoronal(idx:idx2,:,:),1);
    subplot(2,2,3)
     if con==ncondition
        plot(v.freq,mean(logkcoronal,2),'Color',colors(con,:))
    else
    plot(v.freq,mean(logkcoronal,2),'Color',colors(con,:))
     end
    hold on
end
title('Kappa coronal')
xlabel('Frequency [Hz]')
ylabel('Power frecuency')
hold off

for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/poweracc',condition(con,:),'.mat');
    v=load(d,'-mat');
    logkhorizontal=v.khorizontal;
    varkhorizontal=std(log(v.khorizontal),0,2);
    idx=find(v.freq>5,1,'First');
    idx2=find(v.freq>50,1,'First');
    power(con,:)=trapz(v.freq(idx:idx2),v.khorizontal(idx:idx2,:,:),1);
    subplot(2,2,4)
     if con==ncondition
         
     p1=plot(v.freq,mean(logkhorizontal,2),'Color',colors(con,:));
     else
        
    p2=plot(v.freq,mean(logkhorizontal,2),'Color',colors(con,:));
     end
    hold on
    
    
%     ax = get(fig,'CurrentAxes');
% set(ax,'XScale','log','YScale','log')
end
title('Kappa Horizontal')
xlabel('Frequency [Hz]')

ylabel('Power frecuency')
hold off
figure
plot(mean(power,2))
hold on 
errorbar([1:size(power,2)],mean(power,2),std(power,[],2))
set(gca,...
 'XTickLabel',condition,'XTick',[1:size(condition,1)])
end