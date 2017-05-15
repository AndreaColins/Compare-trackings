function Slipstick(slipsver)
%%%%this version of PCA analyse data in exploratory/touch periods, i.e,
%%%%times series may have different lenghts

pole=char('Toothpick','Bamboo','Carbon pole','Wood','Smooth pole','Closed coil','air','Carbon pole','Open coil','Cardboard','Black Sandpaper');
%pole=char('Smooth pole','Open coil');
pole2=char('air','Wood','Closed coil','Sandpaper','Carbon pole','Open coil','Smooth','Bamboo','Toothpick','Cardboard');
npole=size(pole,1);
colors=distinguishable_colors(npole);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
tazimuth=[];
televation=[];
tkcoronal=[];
tkhorizontal=[];
slipstexture=zeros(1,npole);
%%%%%%for measuring length of touch and exploratory periods
lengthtouch=[];
lengthexpl=[];
figprot1=figure;
figprot2=figure;
figprot3=figure;
%figprot4=figure;
figthresh=figure;
n=3486*10;
fromn=3486*0+1;
for i=1:npole
    %%%%%%%%%%%%%%%raw data
    filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
        v=load(filename,'-mat');
        result=v.result;
        nfiles(i)=size(result,3);
        result=diff(diff(result,1,1),1,1);
        
        %%%%%calculate threshold for slipstick as 4 standard deviations
        
        poleaz=[poleaz;reshape(squeeze(result(:,1,:)),nfiles(i)*size(result,1),1)];
        thresholdaz=[mean(poleaz)+4*std(poleaz);mean(poleaz)-4*std(poleaz)];
        
        poleel=[poleel;reshape(squeeze(result(:,2,:)),nfiles(i)*size(result,1),1)];
        thresholdel=[mean(poleel)+4*std(poleel);mean(poleel)-4*std(poleel)];
        
        polekcor=[polekcor;reshape(squeeze(result(:,3,:)),nfiles(i)*size(result,1),1)];
        thresholdkcor=[mean(polekcor)+4*std(polekcor);mean(polekcor)-4*std(polekcor)];
        
        polekhor=[polekhor;reshape(squeeze(result(:,4,:)),nfiles(i)*size(result,1),1)];
        thresholdkhor=[mean(polekhor)+4*std(polekhor);mean(polekhor)-4*std(polekhor)]
        
        
     figure(figthresh)
     ax=subplot(4,1,1);
     plot(poleaz(fromn:n))
        hold on
        plot(thresholdaz(1,:).*ones(1,size(poleaz(fromn:n),1)),'r')
        plot(thresholdaz(2,:).*ones(1,size(poleaz(fromn:n),1)),'g')
        plot(slipsver(:,1),slipsver(:,2),'*g')
        title('Azimuth acceleration')
        xlabel('Time [ms]')
        ylabel('Acceleration [\circ/ms^2]')
        hold off
        ax2=subplot(4,1,2);
        plot(poleel(fromn:n))
        hold on
        plot(thresholdel(1,:).*ones(1,size(poleaz(fromn:n),1)),'r')
        plot(thresholdel(2,:).*ones(1,size(polekhor(fromn:n),1)),'g')
        plot(slipsver(:,1),slipsver(:,2),'*g')
        title('Elevation acceleration')
        xlabel('Time [ms]')
        ylabel('Acceleration [\circ/ms^2]')
        hold off
        ax3=subplot(4,1,3);
        plot(polekcor(fromn:n))
        hold on
        plot(thresholdkcor(1,:).*ones(1,size(polekcor(fromn:n),1)),'r')
        plot(thresholdkcor(2,:).*ones(1,size(polekcor(fromn:n),1)),'g')
        plot(slipsver(:,1),slipsver(:,2)*5*10^-4,'*g')
        title('Vertical curvature acceleration')
        xlabel('Time [ms]')
        ylabel('Acceleration')
        hold off
        ax4=subplot(4,1,4);
        plot(polekhor(fromn:n))
        hold on
        plot(thresholdkhor(1,:).*ones(1,size(polekhor(fromn:n),1)),'r')
        plot(thresholdkhor(2,:).*ones(1,size(polekhor(fromn:n),1)),'g')
        plot(slipsver(:,1),slipsver(:,2)*5*10^-4,'*g')
        title('Horizontal curvature acceleration')
        xlabel('Time [ms]')
        ylabel('Acceleration')
        hold off
        linkaxes([ax,ax2,ax3,ax4],'x') 
end
        

poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
%%
       %%%%find slipsticks in azimuth
   for i=1:npole
        filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
        v=load(filename,'-mat');
        result=v.result;
        nfiles(i)=size(result,3);
        result=diff(diff(result,1,1),1,1);
        slipsazprot=[];
        slipsazret=[];
        poleaz=reshape(squeeze(result(:,1,:)),nfiles(i)*size(result,1),1);
        slips=findslips(poleaz,thresholdaz);
        nprot=1;
        nret=1;
        slipstextures(i)=size(slips,1);
        for j=1:size(slips,1)
            if (slips(j,1)-2)>=1&& (slips(j,1)+10)<size(poleaz,1)
                if poleaz(slips(j,1))>0
                    slipsazprot(:,nprot)=poleaz(slips(j,1)-2:slips(j,1)+10);
                    
                    nprot=nprot+1;
                else
                    slipsazret(:,nret)=poleaz(slips(j,1)-2:slips(j,1)+10);
                    nret=nret+1;
                end      
            end
                  
        end
                    figure(figprot1);
                    subplot(2,1,1)
                   
                    plot([-2:1:10],mean(slipsazprot,2),'Color',colors(i,:))
                    hold on
                
                    plot([-2:1:10],mean(slipsazret,2),'Color',colors(i,:))
                    hold on
   end
    
 for i=1:1
        filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
        v=load(filename,'-mat');
        result=v.result;
        nfiles(i)=size(result,3);
        result=diff(diff(result,1,1),1,1);
        slipsazprot=[];
        slipsazret=[];
        poleaz=reshape(squeeze(result(:,1,:)),nfiles(i)*size(result,1),1);
        slipsbin(:,1)=1:1:size(poleaz);
        slipsbin(:,2)=poleaz>thresholdaz(1)|poleaz<thresholdaz(2);
        poleel=reshape(squeeze(result(:,2,:)),nfiles(i)*size(result,1),1);
        slipsbin(:,3)=poleel>thresholdel(1,1)|poleel<thresholdel(2,1);
        polekcor=reshape(squeeze(result(:,3,:)),nfiles(i)*size(result,1),1);
        slipsbin(:,4)=polekcor>thresholdkcor(1,1)|polekcor<thresholdkcor(2,1);
        polekhor=reshape(squeeze(result(:,4,:)),nfiles(i)*size(result,1),1);
        slipsbin(:,5)=polekhor>thresholdkhor(1,1)|polekhor<thresholdkhor(2,1);
        
        %%%plot
        figure(figprot2)
        ax=subplot(4,1,1);
        plot(slipsbin(fromn:n,2))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title('Azimuth acceleration')
        xlabel('Time [ms]')
        hold off
        ax2=subplot(4,1,2);
        plot(slipsbin(fromn:n,3))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title('Elevation acceleration')
        xlabel('Time [ms]')
        hold off
        ax3=subplot(4,1,3);
        plot(slipsbin(fromn:n,4))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title('Vertical curvature acceleration')
        xlabel('Time [ms]')
        hold off
        ax4=subplot(4,1,4);
        hold on
        plot(slipsbin(fromn:n,5))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title('Horizontal acceleration')
        xlabel('Time [ms]')
        hold off
        linkaxes([ax,ax2,ax3,ax4],'x') 
        %%%%%%%
        %%%%%%plot combination of variables
        %%%%allvariables
        figure(figprot3)
        allvariables=slipsbin(:,2)&slipsbin(:,3)&slipsbin(:,4)&slipsbin(:,5);
        ax1=subplot(4,1,1);
        plot(allvariables(fromn:n,1))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title('All variables')
        xlabel('Time [ms]')
        hold off
       
        
        %%%%atleast 3
        threevar=(slipsbin(:,2)+slipsbin(:,3)+slipsbin(:,4)+slipsbin(:,5))>=3;
        ax2=subplot(4,1,2);
        plot(threevar(fromn:n,1))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title(' At least three variables')
        xlabel('Time [ms]')
        hold off
        %%%%angles+kappa horizontal
        %threevar=(slipsbin(:,2)&slipsbin(:,5))|(slipsbin(:,2)&circshift(slipsbin(:,5),1))|(slipsbin(:,2)&circshift(slipsbin(:,5),-1));
        threevar=(slipsbin(:,2)&slipsbin(:,5));
        ax3=subplot(4,1,3);
        plot(threevar(fromn:n,1))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title('Azimuth & Kappa horizontal acceleration')
        xlabel('Time [ms]')
        hold off
        %%%%atleast 3
        threevar=slipsbin(:,3)&slipsbin(:,4);
        ax4=subplot(4,1,4);
        plot(threevar(fromn:n,1))
        hold on
        plot(slipsver(:,1),slipsver(:,2)/3,'*g')
        title(' Elevation & Kappa coronal acceleration')
        xlabel('Time [ms]')
        hold off
        
        
        delete1=slipsbin(:,2)|slipsbin(:,3)|slipsbin(:,4)|slipsbin(:,5);
        delete=find(delete1==0);
        slipsbin(delete,:)=[]
        k=1;
        size(slipsbin,1)
        for j=1:10
                
         h=find(slipsbin(1:end,1)<=j*3486,1,'last')
         j
         if ~isempty(h)&(h>=k)
         xlswrite('Slipsticksbin.xlsx',slipsbin(k:h,:),j,'A1');
         k=h+1
         end  
        end
        clear slipsbin;
        
        
        
   end  
   
   
   
%    figure(figprot1);
%    subplot(2,1,1)
%    legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on protraction Slip/stick [Azimuth angle]')
%    subplot(2,1,2)
%   legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on retraction Slip/stick [Azimuth angle]')
 
 %%
%%%%find slipsticks in elevation
%    for i=1:npole
%         filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
%         v=load(filename,'-mat');
%         result=v.result;
%         nfiles(i)=size(result,3);
%         result=diff(diff(result,1,1),1,1);
%         slipsazprot=[];
%         slipsazret=[];
%         poleel=reshape(squeeze(result(:,2,:)),nfiles(i)*size(result,1),1);
%         slips=findslips(poleel,thresholdel);
%         nprot=1;
%         nret=1;
%         slipstextures(i)=size(slips,1);
%         for j=1:size(slips,1)
%             if (slips(j,1)-2)>=1&& (slips(j,1)+10)<size(poleel,1)
%                 if poleel(slips(j,1))>0
%                     slipselprot(:,nprot)=poleaz(slips(j,1)-2:slips(j,1)+10);
%                     
%                     nprot=nprot+1;
%                 else
%                     slipselret(:,nret)=poleaz(slips(j,1)-2:slips(j,1)+10);
%                     nret=nret+1;
%                 end      
%             end
%                   
%         end
%                     figure(figprot2);
%                     subplot(2,1,1)
%                     size(mean(slipselprot,2))
%                     plot([-2:1:10],mean(slipselprot,2),'Color',colors(i,:))
%                     hold on
%                 
%                     subplot(2,1,2)
%                     plot([-2:1:10],mean(slipselret,2),'Color',colors(i,:))
%                     hold on
%    end
%      
%    figure(figprot2);
%    subplot(2,1,1)
%    legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on elevation Slip/stick [Elevation angle]')
%    subplot(2,1,2)
%   legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on  Slip/stick [Elevation angle]')
% figure
% bar(slipstextures)
% set(gca,...
% 'XTickLabel',pole2,'XTick',[1:npole])   
% ylabel('Number of slips/sticks')

%%
       %%%%find slipsticks in kappa coronal
%    for i=1:npole
%         filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
%         v=load(filename,'-mat');
%         result=v.result;
%         nfiles(i)=size(result,3);
%         result=diff(diff(result,1,1),1,1);
%         slipsazprot=[];
%         slipsazret=[];
%         polekcor=reshape(squeeze(result(:,3,:)),nfiles(i)*size(result,1),1);
%         slips=findslips(polekcor,thresholdkcor);
%         nprot=1;
%         nret=1;
%         slipstextures(i)=size(slips,1);
%         for j=1:size(slips,1)
%             if (slips(j,1)-2)>=1&& (slips(j,1)+10)<size(polekcor,1)
%                 if polekcor(slips(j,1))>0
%                     slipskcorprot(:,nprot)=polekcor(slips(j,1)-2:slips(j,1)+10);
%                     
%                     nprot=nprot+1;
%                 else
%                     slipskcorret(:,nret)=polekcor(slips(j,1)-2:slips(j,1)+10);
%                     nret=nret+1;
%                 end      
%             end
%                   
%         end
%                     figure(figprot3);
%                     subplot(2,1,1)
% 
%                     plot([-2:1:10],mean(slipskcorprot,2),'Color',colors(i,:))
%                     hold on
%                 
%                     subplot(2,1,2)
%                     plot([-2:1:10],mean(slipskcorret,2),'Color',colors(i,:))
%                     hold on
%    end
     
%    figure(figprot3);
%    subplot(2,1,1)
%    legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on coronal curvature Slip/stick [Kappa coronal]')
%    subplot(2,1,2)
%   legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on coronal curvature Slip/stick [Kappa coronal]')
%    
   %%
       %%%%find slipsticks in khorizontal
%    for i=1:npole
%         filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
%         v=load(filename,'-mat');
%         result=v.result;
%         nfiles(i)=size(result,3);
%         result=diff(diff(result,1,1),1,1);
%         slipsazprot=[];
%         slipsazret=[];
%         polekhor=reshape(squeeze(result(:,4,:)),nfiles(i)*size(result,1),1);
%         slips=findslips(polekhor,thresholdkhor);
%         %counters for protraction and retraction slipssticks 
%         nprot=1;
%         nret=1;
%         slipstextures(i)=size(slips,1);
%         for j=1:size(slips,1)
%             if (slips(j,1)-2)>=1&& (slips(j,1)+10)<size(polekhor,1)
%                 if polekhor(slips(j,1))>0
%                     slipskhorprot(:,nprot)=polekhor(slips(j,1)-2:slips(j,1)+10);
% %                     figure(figprot4);
% %                     subplot(2,1,1)
% % 
% %                     plot([-2:1:10],slipskhorprot(:,nprot),'Color',colors(i,:))
% %                     hold on
%        
%                     nprot=nprot+1;
%                 else
%                     slipskhorret(:,nret)=polekhor(slips(j,1)-2:slips(j,1)+10);
%                     subplot(2,1,2)
%                     plot([-2:1:10],slipskhorret(:,nret),'Color',colors(i,:))
%                     hold on
%                     nret=nret+1;
%                 end      
%             end
%                   
%         end
%         
%         
% % %                     figure(figprot4);
% % %                     subplot(2,1,1)
% % % 
% % %                     plot([-2:1:10],mean(slipskhorprot,2),'Color',colors(i,:))
% % %                     hold on
% % %                 
% % %                     subplot(2,1,2)
% % %                     plot([-2:1:10],mean(slipskhorret,2),'Color',colors(i,:))
% % %                     hold on
%     end
%      figure
%  bar(slipstextures)
%  set(gca,...
%  'XTickLabel',pole2,'XTick',[1:npole])   
%  ylabel('Number of slips/sticks')
%    figure(figprot4);
%    subplot(2,1,1)
%    legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on protraction Slip/stick [Curvature horizontal]')
%    subplot(2,1,2)
%   legend(pole2)
%    xlabel('Time [ms]')
%    ylabel('Acceleration [\circ/ms^2]')
%    title('Mean acceleration on retraction Slip/stick [Curvature horizontal]')
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
        fin=n;
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
if I2==3486
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