function powerspectra2(result,meanplot)
%%%%%%%%%%%%using fft
%%%%%%%%%%%%meanplot==0 shows every curve in result
%%%%%%%%%%%meanplot/=0 shows the mean and standar deviation curve for every variable 
N = size(result,1); 
T = N/1000; %% define time of interval, 3.4 seconds
t = [0:N-1]/N; 
t = t*T; 
f = result; 
for j=1:size(f,3)
for i=1:size(f,2)
    tff=fft(f(:,i,j));
%     tff2=tff;
%     threshold = max(abs(tff))/1000; %tolerance threshold
% tff2(abs(tff)<threshold) = 0;
power(:,i,j) = abs(tff)/(N/2); %% absolute value of the fft
phase(:,i,j) = mod(unwrap(angle(tff))*180/pi,360);
end
end
phase = phase(1:N/2,:,:);
power = power(1:N/2,:,:); %% take the power of positve freq. half
freq = [0:N/2-1]/T; 

if meanplot==0
    figure
h1(1)=subplot(4,2,1);
semilogy(freq,squeeze(power(:,1,:)),'b'); 
h2(1)=subplot(4,2,2);
plot(freq,squeeze(phase(:,1,:)),'b'); 
h1(2)=subplot(4,2,3);
semilogy(freq,squeeze(power(:,2,:)),'g');
h2(2)=subplot(4,2,4);
plot(freq,squeeze(phase(:,2,:)),'g'); 
h1(3)=subplot(4,2,5);
semilogy(freq,squeeze(power(:,3,:)),'r'); 
h2(3)=subplot(4,2,6);
plot(freq,squeeze(phase(:,3,:)),'r'); 
h1(4)=subplot(4,2,7);
semilogy(freq,squeeze(power(:,4,:)),'c'); 
h2(4)=subplot(4,2,8);
plot(freq,squeeze(phase(:,4,:)),'c'); 
else
    figure
    h1(1)=subplot(4,2,1);
shadedErrorBar(freq,mean(log(squeeze(power(:,1,:))),2),std(log(squeeze(power(:,1,:))),0,2),'b'); 
h2(1)= subplot(4,2,2);
shadedErrorBar(freq,mean(squeeze(phase(:,1,:)),2),std(squeeze(phase(:,1,:)),0,2),'b'); 
h1(2)=subplot(4,2,3);
shadedErrorBar(freq,mean(log(squeeze(power(:,2,:))),2),std(log(squeeze(power(:,2,:))),0,2),'g'); 
h2(2)=subplot(4,2,4);
shadedErrorBar(freq,mean(squeeze(phase(:,2,:)),2),std(squeeze(phase(:,2,:)),0,2),'g'); 
h1(3)=subplot(4,2,5);
shadedErrorBar(freq,mean(log(squeeze(power(:,3,:))),2),std(log(squeeze(power(:,3,:))),0,2),'r'); 
h2(3)=subplot(4,2,6);
shadedErrorBar(freq,mean(squeeze(phase(:,3,:)),2),std(squeeze(phase(:,3,:)),0,2),'r'); 
h1(4)=subplot(4,2,7);
shadedErrorBar(freq,mean(log(squeeze(power(:,4,:))),2),std(log(squeeze(power(:,4,:))),0,2),'c'); 
h2(4)=subplot(4,2,8);
shadedErrorBar(freq,mean(squeeze(phase(:,4,:)),2),std(squeeze(phase(:,4,:)),0,2),'c'); 
end

 linkaxes(h2,'x')
    linkaxes(h1,'x')
%legend('Azimuth','Elevation','Kappa Coronal','Kappa Horizontal')
%axis([0 20 0 1]); %% zoom in
% figure
% semilogy(freq',p./[freq',freq',freq',freq'])
% legend('Azimuth','Elevation','Kappa Coronal','Kappa Horizontal')
end