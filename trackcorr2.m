function [result,v,power]=trackcorr2(files,meanplot)
%meanplot=0;
[idx,result]=comparetr4(files,0);
%result=[idx,sin(idx*(1/100)+pi/4),sin(idx.*(idx*0.00001+1/100))+rand(size(idx)),sin(idx.*(idx*0.00001+1/100))];
nfiles=size(files,1);
acor=zeros(size(xcov(result(1:end,:,1)),1),16,nfiles);


for i=1:nfiles
[acor(:,:,i),lag]=xcov(result(1:end,:,i),'coeff');
%[~,idlag(i,:)]=max(abs(acor(:,:,i)));
end
v.autoazimuth=squeeze(acor(:,[1],:));
v.autoelevation=squeeze(acor(:,[6],:));
v.autokcoronal=squeeze(acor(:,[11],:));
v.autokhorizontal=squeeze(acor(:,[16],:));
if meanplot==0
figure(1)
        
        hold on 
         subplot(4,2,1)
       
         plot(idx,squeeze(result(:,1,:)),'b')
        
          title('Azimuth')
          xlabel('Time [msec]')
          ylabel('Angle [degrees]')
         
         
       h1(1)=subplot(4,2,2);
        
        plot(lag, v.autoazimuth,'b')
       
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance azimuth')
        
        subplot(4,2,3)
         plot(idx,squeeze(result(:,2,:)),'g')
        
          title('Elevation')
          xlabel('Time [msec]')
          ylabel('Angle [degrees]')
          
      h1(2)=subplot(4,2,4);
       
        plot(lag,v.autoelevation,'g')
        
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance elevation')
        
        subplot(4,2,5)
         plot(idx,squeeze(result(:,3,:)),'r')
        
          title('Kappa Coronal')
          xlabel('Time [msec]')
          %ylabel('')
          
          
      h1(3)=subplot(4,2,6);
     
        plot(lag,v.autokcoronal,'r')
       
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance kappa coronal')
     
        subplot(4,2,7)
         plot(idx,squeeze(result(:,4,:)),'c')
          title('Kappa Horizontal')
          xlabel('Time [msec]')
          %ylabel('')
          
        h1(4)=subplot(4,2,8);
       
        plot(lag,v.autokhorizontal,'c')
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance kappa horizontal')
      
        hold off
       
        power=powerspectra3(result,0)
        powerspectra2(result,0)
        
else
  figure(1)
        
        hold on 
         subplot(4,2,1)
         
         %plot(idx,mean(squeeze(result(:,1,:)),2),'b')
        shadedErrorBar(idx,mean(squeeze(result(:,1,:)),2),std(squeeze(result(:,1,:)),0,2),'-b')
          title('Azimuth')
          xlabel('Time [msec]')
          ylabel('Angle [degrees]')
         
         
       subplot(4,2,2)
       % plot(lag,mean(squeeze(acor(:,[1],:)),2),'b')
        shadedErrorBar(lag,mean(squeeze(acor(:,[1],:)),2),std(squeeze(acor(:,[1],:)),0,2),'-b')
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance azimuth')
        
        subplot(4,2,3)
         %plot(idx,squeeze(result(:,2,:)),'g')
        shadedErrorBar(idx,mean(squeeze(result(:,2,:)),2),std(squeeze(result(:,2,:)),0,2),'-g')
          title('Elevation')
          xlabel('Time [msec]')
          ylabel('Angle [degrees]')
          
      subplot(4,2,4)
        shadedErrorBar(lag,mean(squeeze(acor(:,[6],:)),2),std(squeeze(acor(:,[6],:)),0,2),'-g')
        
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance elevation')
        
        subplot(4,2,5)
        % plot(idx,squeeze(result(:,3,:)),'r')
        shadedErrorBar(idx,mean(squeeze(result(:,3,:)),2),std(squeeze(result(:,3,:)),0,2),'-r')
          title('Kappa Coronal')
          xlabel('Time [msec]')
          %ylabel('')
          
          
      subplot(4,2,6)
       % plot(lag,squeeze(acor(:,[11],:)),'r')
       shadedErrorBar(lag,mean(squeeze(acor(:,[11],:)),2),std(squeeze(acor(:,[11],:)),0,2),'-r')
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance kappa coronal')
     
        subplot(4,2,7)
         %plot(idx,squeeze(result(:,4,:)),'c')
         shadedErrorBar(idx,mean(squeeze(result(:,4,:)),2),std(squeeze(result(:,4,:)),0,2),'-c')
          title('Kappa Horizontal')
          xlabel('Time [msec]')
          %ylabel('')
          
        subplot(4,2,8)
        %plot(lag,squeeze(acor(:,[16],:)),'c')
        shadedErrorBar(lag,mean(squeeze(acor(:,[16],:)),2),std(squeeze(acor(:,[16],:)),0,2),'-c')
         xlabel('lag')
        ylabel('autocovariance')
        title('Autocovariance kappa horizontal')
      
        hold off      
        %%%%%%%%%%%%%%%%%%%%%%%%%power spectrum plot

        powerspectra2(result,1)
        power=powerspectra3(result,1)
end

        
      
        
end