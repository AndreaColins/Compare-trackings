%%%setpoint of <5Hz
function setpointtextures 
pole=char('Smooth pole','Closed coil', 'Open coil','Black sandpaper','Carbon pole','Toothpick','Bamboo','Wood','Cardboard');
pole2=pole;
colors={'b','r','g','k','c','m','y'};
npole=size(pole,1);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
resultair=[];
azimuth=[];
elevation=[];
kcoronal=[];
khorizontal=[];
amplitude=[];
amplitudeair=[];
for i=1:npole
   
%%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%data between 0 to 5 for setpoint
    %filename=strcat('./videosselectedtextures/processed/bandpass/upto5/result',pole(i,:),'.mat');
    
    %%%file result(pole)1to5.mat actually has the results for filtered
    %%%data between 0.01 to 5 for hilbert tranform
    %filename=strcat('./videosselectedtextures/processed/bandpass/upto5/result',pole(i,:),'1to5.mat');
   d=strcat('./videosselectedtextures/processed/touchdetection/',pole(i,:),'/');
   ff = dir([d '*.mat']);
        nfiles(i)=size(ff,1);
        string_list=cell(nfiles(i),1);
        for j=1:nfiles(i)
            string_list{j}=strcat(d,ff(j).name);
            touchfile=strcat(d,ff(j).name)
            v=load(filename,'-mat');
            result=v.result;
            vtouch=load(touchfile,'-mat');
            idx=[1:3488];
            idx2=circshift(idx,size(idx,1)-vtouch.start_frame-1,2);
            %idx2=circshift(idx,-(size(idx,1)-vtouch.start_frame-1));
            result=[[0,0,0,0];diff(result(idx,:,j),1,1)];
            vtouch.touches=vtouch.touches(idx2);
            startp=find(vtouch.touches,1,'first');
            endp=find(vtouch.touches,1,'last');
            
%             %%%%%%%%%%%%%%%for exploratory periods
% %             [result,touchperiods]=centredk(endp,vtouch.touches,result);
% %             if (3488-endp)>=50
% %                 endp=endp+50;
% %             else
% %                 endp=3488;
% %             end
% %             result=result(startp:endp,:);
%             
  %%%%%%%%%%%%for touch periods
  
            I=find(vtouch.touches);
            result=result(I,:);

  
            h1=hilbert(result);
            %[h1,h2]=envelope(result);
             plot(result(:,1))
%             hold on 
%             envelope(result(:,1),100)
%             pause(2)
%             hold off
            amplitude=[amplitude;abs(h1)];
            azimuth=[azimuth;result(:,1)];
            histcounts([azimuth],80)
            pause(2)
            elevation=[elevation;result(:,2)];
            kcoronal=[kcoronal;result(:,3)];
            khorizontal=[khorizontal;result(:,4)];
        end
  
end
subplot(2,2,1)
[~,edges]=histcounts([azimuth],80);
hold on
histogram(azimuth,edges,'Normalization','pdf');
legend(pole2)
xlabel('Angle (\circ)')
title('Azimuth')
hold off
subplot(2,2,2)
[~,edges]=histcounts([elevation],80);
hold on

histogram(elevation,edges,'Normalization','pdf');
legend(pole2)
xlabel('Angle (\circ)')
title('Elevation')
hold off
subplot(2,2,3)
[~,edges]=histcounts([kcoronal],80);
hold on

histogram(kcoronal,edges,'Normalization','pdf');
title('Kappa Coronal')
legend(pole2)
hold off
subplot(2,2,4)
[~,edges]=histcounts([khorizontal],80);
hold on

histogram(khorizontal,edges,'Normalization','pdf');
legend(pole2)
title('Kappa Horizontal')
hold off

figure
subplot(2,2,1)
[~,edges]=histcounts([squeeze(amplitude(:,1,:))],80);
hold on

histogram(squeeze(amplitude(:,1,:)),edges,'Normalization','pdf');
legend(pole2)
xlabel('Angle (\circ)')
title('Azimuth')
subplot(2,2,2)
[~,edges]=histcounts([squeeze(amplitude(:,2,:))],40);
hold on

histogram(squeeze(amplitude(:,2,:)),edges,'Normalization','pdf');
legend(pole2)
xlabel('Angle (\circ)')
title('Elevation')
subplot(2,2,3)
[~,edges]=histcounts([squeeze(amplitude(:,3,:))],40);
hold on

histogram(squeeze(amplitude(:,3,:)),edges,'Normalization','pdf');
title('Kappa Coronal')
legend(pole2)
subplot(2,2,4)
[~,edges]=histcounts([squeeze(amplitude(:,4,:))],40);
hold on

histogram(squeeze(amplitude(:,4,:)),edges,'Normalization','pdf');
legend(pole2)
title('Kappa Horizontal')
end
function [result,touchperiods]=centredk(endp,touches,result)
i=1;
touchescopy=touches;
touchcount=0;
while i<endp
    ini=find(touchescopy,1,'first');
    touchescopy(1:ini)=1;
    fin=find(~touchescopy,1,'first')-1;
     if isempty(fin)
        fin=3488;
    end
    touchescopy(1:fin)=0;
    hold on
    touchperiods(touchcount+1,1:2)=[ini,fin];
    touchcount=touchcount+1;
    i=fin;
end
%touchperiods

kappacoronal=result(:,3,1);
kappahorizontal=result(:,4,1);
kappacoronal(~touches)=0;
kappahorizontal(~touches)=0;
for i=1:touchcount
    if i==1
        if touchperiods(1,1)>10
            k0coronal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,4));
        else
            k0coronal=median(result(1:touchperiods(i,1)-1,3));
            k0horizontal=median(result(1:touchperiods(i,1)-1,4));
        end
    else
        if touchperiods(i,1)-touchperiods(i-1,2)>10
            k0coronal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,4));
        else
            k0coronal=median(result(touchperiods(i-1,2)+1:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i-1,2)+1:touchperiods(i,1)-1,4));
        end
    end
        kappacoronal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),3)-k0coronal;
    kappahorizontal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),4)-k0horizontal;
   
end
result(:,3)=kappacoronal;
result(:,4)=kappahorizontal;
clear touchescopy ini fin
end