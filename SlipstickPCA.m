function SlipstickPCA
%%%%this version of PCA analyse data in exploratory/touch periods, i.e,
%%%%times series may have different lenghts

pole=char('air','Smooth pole','Black Sandpaper','Bamboo','Wood','Closed coil','Open coil','Cardboard','Toothpick','Carbon pole');
%pole=char('Smooth pole','Open coil');
%pole2=char('Smooth', 'Carbon', 'Sandpaper', 'Closed','Open','Cardboard','Bamboo','Toothpick','Wood');
pole2=pole;
npole=size(pole,1);
colors=distinguishable_colors(npole);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
ms1=50;%number of miliseconds
overlap=25;%miliseconds
poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
tazimuth=[];
televation=[];
tkcoronal=[];
tkhorizontal=[];
%%%%%%for measuring length of touch and exploratory periods
lengthtouch=[];
lengthexpl=[];
fig=figure
for i=1:npole
    i
    ms=ms1;%number of miliseconds
    %%%%%%%%%%%%%%%raw data
    filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    %filename=strcat('./videosselectedtextures/processed/bandpass/from20/result',pole(i,:),'.mat');
    if strcmp(pole(i,1:3),'air')
        v=load(filename,'-mat');
        result=v.result;
        nfiles(i)=size(result,3);
        ms=ms-overlap;
        nresult(i)=floor(size(result,1)/ms);
        result=diff(diff(result,1,1));
        result((nresult(i)*ms+1):end,:,:)=[];
        
        
        poleaz=reshape(squeeze(result(:,1,:)),nfiles(i).*nresult(i).*ms,1);
        poleaz=reshape(poleaz,ms,nresult(i).*nfiles(i));

        poleel=reshape(squeeze(result(:,2,:)),nfiles(i).*nresult(i).*ms,1);
        poleel=reshape(poleel,ms,nresult(i).*nfiles(i));

           thresholdh=poleaz>1.3343 | poleaz<-1.3342;

%             thresholdh=thresholdh|[ zeros(3,size(thresholdh,2)); thresholdh(1:end-3,:)];
%             poleaz=thresholdh.*poleel;
%             thresholdv=poleel>1.2211 | poleel<-1.2213;
%            thresholdv=thresholdv|[ zeros(3,size(poleaz,2)); thresholdv(1:end-3,:)];
%             poleel=thresholdv.*poleel;
        
        
        
        polekcor=reshape(squeeze(result(:,3,:)),nfiles(i).*nresult(i).*ms,1)*0;
        polekcor=reshape(polekcor,ms,nresult(i).*nfiles(i))*0;
   
        polekhor=reshape(squeeze(result(:,4,:)),nfiles(i).*nresult(i).*ms,1)*0;
        polekhor=reshape(polekhor,ms,nresult(i).*nfiles(i))*0;
         
        if overlap~=0
            for j=1:nresult(i)*nfiles(i)-1
                poleaz2(:,j)=[poleaz(:,j);poleaz(1:overlap,j+1)];
                poleel2(:,j)=[poleel(:,j);poleel(1:overlap,j+1)];
                polekcor2(:,j)=[polekcor(:,j);polekcor(1:overlap,j+1)]*0;
                polekhor2(:,j)=[polekhor(:,j);polekhor(1:overlap,j+1)]*0;
            end
           
            poleaz= poleaz2;
            poleel=poleel2;
            polekcor=polekcor2*0;
            polekhor=polekhor2*0;
        end
        nresult(i)=size(poleaz,2);
    else 
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
            result=result(idx,:,j);
            vtouch.touches=vtouch.touches(idx2);
            %%%%%%%select touch periods
            startp=find(vtouch.touches,1,'first');
            endp=find(vtouch.touches,1,'last');

            %%%%%%%%%%%%%%%for exploratory periods
%             [result2,touchp]=centredk(endp,vtouch.touches,result);
%             if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
%             touchp=[startp endp];
            %result=result2(startp:endp,:);
%             
%%%%%%%%%for debugging
           
          
%             [result2,touchp]=centredk(endp,vtouch.touches,result);
%              if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
% %             result=notouchselection(vtouch.touches,result);
%             result=touchselection(vtouch.touches,result);
%                          hold on
%                          %result=result(startp:endp,:);
%                          plot(result2(:,4),'r')
%                        
%                          hold off
%                          pause(0.1)
            %%%%%%%%%%%%for touch periods
             
            [result2,touchp]=centredk(endp,vtouch.touches,result);
            
%            hold on
%            plot(result2(:,4),'r')
           diffkhor=diff(diff(result2,1,1));
           
%            plot(diff(diffkhor),'k')
%            pause(2)
%            hold off

            diffkhor=[[0,0,0,0];diffkhor];
%             thresholdh=diffkhor(:,1)>1.3343 | diffkhor(:,1)<-1.3342;
%             thresholdh=thresholdh|[ zeros(3,1); thresholdh(1:end-3)];
%             diffkhor(:,1)=thresholdh.*diffkhor(:,1);
%             thresholdv=diffkhor(:,2)>1.2211 | diffkhor(:,2)<-1.2213;
%            thresholdv=thresholdv|[ zeros(3,1); thresholdv(1:end-3)];
%             diffkhor(:,2)=thresholdv.*diffkhor(:,2);
%             thresholdv=diffkhor(:,3)>1.0e-03*0.5076 | diffkhor(:,3)<-1.0e-03*0.5077;
%             thresholdv=thresholdv|[ zeros(3,1); thresholdv(1:end-3)];
%             diffkhor(:,3)=thresholdv.*diffkhor(:,4);
%             thresholdk=diffkhor(:,4)<1.0e-03*0.4437 |diffkhor(:,4)<-1.0e-03*0.4437;
%             thresholdv=thresholdv|[ zeros(3,1); thresholdv(1:end-3)];
%             diffkhor(:,4)=thresholdv;
             result2=diffkhor;
%             subplot(4,1,1)
%             plot(result2(:,1))
%             subplot(4,1,2)
%             plot(result2(:,2))
%             subplot(4,1,3)
%             plot(result2(:,3)*10000)
%             subplot(4,1,4)
%             plot(result2(:,4)*10000)
%       
%               clf(fig)
            %%%%%%for measuring length of touch and exploratory periods
            lengthtouch=[lengthtouch; touchp(:,2)-touchp(:,1)];
            lengthexpl=[lengthexpl;touchp(end,2)-touchp(1,1)];
            %discard small touch periods
            counter=1;
            touchperiods=[];
             for t=1:size(touchp,1)
                 if (touchp(t,2)-touchp(t,1))>=ms
                     touchperiods(counter,:)=touchp(t,:);
                     counter=counter+1;
                 end
             end
           
             clear touchp 
             if ~isempty(touchperiods)
            for t=1:size(touchperiods,1)
                if touchperiods(t,2)>=3477
                    touchperiods(t,2)=3466;
                end
                result=result2(touchperiods(t,1):touchperiods(t,2),:);

            nsamples=floor(size(result,1)/(ms-overlap));
            if mod(size(result,1),(ms-overlap))<overlap
                nsamples=nsamples-1;
            end
            nresult(i)=nsamples+nresult(i);
            result((nsamples*(ms-overlap)+overlap+1):end,:)=[];
            
             if overlap==0

                poleaz2=reshape(result(:,1),ms,nsamples);
                poleel2=reshape(result(:,2),ms,nsamples);
                polekcor2=reshape(result(:,3),ms,nsamples);
                polekhor2=reshape(result(:,4),ms,nsamples);
             else
                
                for j=0:nsamples-1
                    poleaz2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,1);

                    poleel2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,2);
                    polekcor2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,3);
                    polekhor2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,4);
                end
                
            end
    
            poleaz=[poleaz, poleaz2];

            poleel=[poleel,poleel2];
            
            polekcor=[polekcor,polekcor2];
            
            polekhor=[polekhor, polekhor2];
            
            clear poleaz2 poleel2 polekcor2 polekhor2 result 
            end
             end
            clear touchperiods
        end 
        
    end
   
   
    
    %%%%%%%%%%%%%%%%%%Add new condition

        tazimuth=[tazimuth,poleaz];
        televation=[televation,poleel];
        tkcoronal=[tkcoronal,polekcor];
        tkhorizontal=[tkhorizontal,polekhor];

    
    poleaz=[];
    poleel=[];
    polekcor=[];
    polekhor=[];
    clear result poleaz2 poleel2 polekcor2 polekhor2 result2
end
%%
%%%%%%for measuring length of touch and exploratory periods
% figure
% subplot(1,2,1)
% histogram(lengthtouch,60)
% xlabel('Length of touch [ms]')
% ylabel('Number of touches')
% title('Touch periods')
% subplot(1,2,2)
% histogram(lengthexpl,20)
% title('Exploratory periods')
% xlabel('Length of exploratory periods [ms]')
% ylabel('Number of exploratory periods')
%%
% curvatures=[tkhorizontal;tkcoronal]';
% curvatures=(curvatures-mean(curvatures(:)))/std(curvatures(:));
% 
 azel=[tazimuth;televation]';
 azel=(azel-mean(azel(:)))/std(azel(:));
% azel=[azel,curvatures];
% azel=[[tkhorizontal',tkcoronal',tazimuth',televation']];
%%%%%%%%%Raw Data

plot(azel(1,:)','b')
hold on 
plot(azel(nresult(1)+1,:)','r')
%plot(azel(nresult(1)+nresult(2)+1,:)','g')
plot(azel(1:nresult(1),:)','b')
plot(azel(nresult(1)+1:nresult(1)+nresult(2),:)','r')
%plot(azel(nresult(1)+nresult(2)+1:end,:)','g')
hold off

%%
% %%%%Apply PCA angles
 [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(azel);

size(SCORE)
%%%%%%%%Plot results
figure
%%%graficar en los nuevos vectores
subplot(2,2,1)
if strcmp(pole(1,1:3),'air')
    totalfiles=nresult(1)
else
    totalfiles=nresult(1)
end
scatter3(SCORE(1:totalfiles,1),SCORE(1:totalfiles,2),SCORE(1:totalfiles,3),[],'MarkerEdgeColor',colors(1,:))
%means(1)=mean(mean(azel(1:totalfiles,1:100)));
for i=2:npole
hold on
totalfiles+nresult(i)
scatter3(SCORE(totalfiles+1:(totalfiles+nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nresult(i)),2),SCORE(totalfiles+1:(totalfiles+nresult(i)),3),[],'MarkerEdgeColor',colors(2,:))
%means(i)=mean(mean(azel(totalfiles+1:(totalfiles+nresult(i)),1:100)));
title('[Azimuth,Elevation]')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nresult(i);
end
legend(pole2,'Location','best','Orientation','horizontal') 
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot Firsts components
subplot(2,2,2)
plot(COEFF(:,1),'r')
hold on
plot(COEFF(:,2),'b')
legend('First Component','Second Component')
xlabel('Time [miliseconds]')
vexplained=EXPLAINED(1)+EXPLAINED(2);
title(strcat('Variance explained=',num2str(vexplained)))


%%%%%%%%%%%%%%%%%%%%%%%%explained variance
subplot(2,2,3)
plot(cumsum(EXPLAINED),'*-')
xlabel('Number of principal component')
ylabel('Variance explained [%]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components

subplot(2,2,4)
scatter(COEFF(:,1),COEFF(:,2))
xlabel('First Component')
ylabel('Second Component')
% %xlab={pole2(3,:),pole2(1,:),pole2(5,:),pole2(2,:),pole2(6,:),pole2(7,:),pole2(4,:)};
% %hold off
% %bar([means(3),means(1),means(5),means(2),means(6),means(7),means(4)])
% % set(gca,...
% % 'XTickLabel',xlab,'XTick',[1:size(pole2,1)])
% %hold off
% %%
% %%%%%%%%%%%%%%Apply PCA azimuth and kappa horizontal
% [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(zscore([tazimuth',tkhorizontal']));

end
function slips=findslips(pole,threshold)
outliers=(pole>=threshold(1,1))|(pole<=threshold(2,1));
i=1;
outlierscopy=outliers;
countslips=0;
slips(1,1:2)=[0,0];
while i<size(outliers,1)-1
    ini=find(outlierscopy,1,'first');
        if isempty(ini) || ini==size(outliers,1)
        break;
        end
    outlierscopy(1:ini)=1;
    fin=find(~outlierscopy,1,'first')-1;
     if isempty(fin)
        fin=size(outliers,1);
    end
    outlierscopy(1:fin)=0;
    hold on
    slips(countslips+1,1:2)=[ini,fin];
    if countslips>1 && slips(countslips+1,1)-slips(countslips,2)<=2
        countslips=countslips;
    else
    countslips=countslips+1;
    end
    i=fin;
end
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
        elseif touchperiods(1,1)==1
            k0coronal=median(result(touchperiods(i+1,1)-10:touchperiods(i+1,1)-1,3));
            k0horizontal=median(result(touchperiods(i+1,1)-10:touchperiods(i+1,1)-1,4));
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
function [result]=notouchselection(touches,result)
% plot(result(:,4))
% hold on
I=find(touches,1,'first');
I2=find(touches,1,'last');
if I2==3488
result=result(1:I,:);
else
result=[result(1:I,:);result(I2:end,:)];
end

clear I I2
end
function [result]=touchselection(touches,result)
% plot(result(:,4))
% hold on
I=find(touches);
result=result(I,:);
% plot(result(:,4),'.r')
% pause(5)
% hold off
clear I
end