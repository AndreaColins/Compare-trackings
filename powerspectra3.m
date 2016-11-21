function powerspectra3(result,meanplot)
%%%%%%%%%using multitapper method
%%%%%%%%%%%%meanplot==0 shows every curve in result
%%%%%%%%%%%meanplot/=0 shows the mean and standar deviation curve for every variable 
N = 3488; 
T = 3.487; %% define time of interval, 3.4 seconds
t = [0:N-1]/N; 
t = t*T; 
f = result; 
for j=1:size(f,3)
for i=1:size(f,2)
    p(:,i,j)=pmtm(f(:,i,j),4);
%p(:,i,j) = abs(fft(f(:,i,j)))/(N/2); %% absolute value of the fft
end
end
p = p(1:N/2,:,:); %% take the power of positve freq. half
freq = [0:N/2-1]/T; 
if meanplot==0
    figure
h1(1)=subplot(4,1,1);
semilogy(freq,squeeze(p(:,1,:)),'b'); 
%pmtm(squeeze(f(:,1,:)))
h1(2)=subplot(4,1,2);
semilogy(freq,squeeze(p(:,2,:)),'g'); 
h1(3)=subplot(4,1,3);
semilogy(freq,squeeze(p(:,3,:)),'r'); 
h1(4)=subplot(4,1,4);
semilogy(freq,squeeze(p(:,4,:)),'c'); 
else
    figure
    h1(1)=subplot(4,1,1)
shadedErrorBar(freq,mean(log(squeeze(p(:,1,:))),2),std(log(squeeze(p(:,1,:))),0,2),'b'); 
h1(2)=subplot(4,1,2);
shadedErrorBar(freq,mean(log(squeeze(p(:,2,:))),2),std(log(squeeze(p(:,2,:))),0,2),'g'); 
h1(3)=subplot(4,1,3);
shadedErrorBar(freq,mean(log(squeeze(p(:,3,:))),2),std(log(squeeze(p(:,3,:))),0,2),'r'); 
h1(4)=subplot(4,1,4);
shadedErrorBar(freq,mean(log(squeeze(p(:,4,:))),2),std(log(squeeze(p(:,4,:))),0,2),'c');  
end

    linkaxes(h1,'x')
%legend('Azimuth','Elevation','Kappa Coronal','Kappa Horizontal')
%axis([0 20 0 1]); %% zoom in
% figure
% semilogy(freq',p./[freq',freq',freq',freq'])
% legend('Azimuth','Elevation','Kappa Coronal','Kappa Horizontal')
end