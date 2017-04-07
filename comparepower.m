function comparepower
tic
condition=char('Smooth pole', 'Carbon Pole', 'Black Sandpaper', 'Closed coil','Open coil','Cardboard','Bamboo','Toothpick','Wood','air');
ncondition=size(condition,1);
logazimuth=[];
varazimuth=[];
logelevation=[];
varelevation=[];
logkcor=[];
varkcor=[];
logkhor=[];
varkhor=[];
for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/power',condition(con,:),'.mat');
    v=load(d,'-mat');
    logazimuth=log(v.azimuth);
    varazimuth=std(log(v.azimuth),0,2);
    freq=log(v.freq);
    subplot(2,2,1)
    if con==ncondition
        plot(log(v.freq),logazimuth,'r')
    else
        plot(log(v.freq),logazimuth,'b')
    end
    hold on
   
end
title('Azimuth')
xlabel('Log(Freq)')
ylabel('Power frecuency')
hold off

for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/power',condition(con,:),'.mat');
    v=load(d,'-mat');
    logelevation=log(v.elevation);
    varelevation=std(log(v.elevation),0,2);
    subplot(2,2,2)
    if con==ncondition
        plot(log(v.freq),logelevation,'r')
    else
        plot(log(v.freq),logelevation,'b')
    end
    hold on
end
title('Elevation')
xlabel('Log(Freq)')
ylabel('Power frecuency')
hold off

for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/power',condition(con,:),'.mat');
    v=load(d,'-mat');
    logkcoronal=log(v.kcoronal);
    varkcoronal=std(log(v.kcoronal),0,2);
    subplot(2,2,3)
     if con==ncondition
        plot(log(v.freq),logkcoronal,'r')
    else
    plot(log(v.freq),logkcoronal,'b')
     end
    hold on
end
title('Kappa coronal')
xlabel('Log(Freq)')
ylabel('Power frecuency')
hold off

for con=1:ncondition
    
    d=strcat('./videosselectedtextures/processed/FrameCorrected/power',condition(con,:),'.mat');
    v=load(d,'-mat');
    logkhorizontal=log(v.khorizontal);
    varkhorizontal=std(log(v.khorizontal),0,2);
   
    subplot(2,2,4)
     if con==ncondition
         
     p1=plot(freq,logkhorizontal,'r');
     else
        
    p2=plot(freq,logkhorizontal,'b');
     end
    hold on
    
%     ax = get(fig,'CurrentAxes');
% set(ax,'XScale','log','YScale','log')
end
title('Kappa Horizontal')
xlabel('Log(Freq)')
ylabel('Power frecuency')
hold off
end