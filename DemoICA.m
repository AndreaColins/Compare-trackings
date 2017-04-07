function DemoICA
filename=strcat('./videosselectedtextures/processed/FrameCorrected/result','air','.mat');
v=load(filename);
result=v.result;
%mixedsig=squeeze(result(:,1,:));
mixedsig=reshape(squeeze(result(1:3400,1,:)),3400/10,10*10);
%[sig,mixedsig]=demosig();
figure
subplot(4,1,1)
plot(mixedsig(1:3,:)')
hold on 
title('Mixed signals')
hold off
subplot(4,1,2)
plot(mixedsig(4:6,:)')
subplot(4,1,3)
plot(mixedsig(6:8,:)')
subplot(4,1,4)
plot(mixedsig(9:10,:)')

% figure
% subplot(4,1,1)
% plot(sig(1,:))
% hold on 
% title('Original signals')
% hold off
% subplot(4,1,2)
% plot(sig(2,:))
% subplot(4,1,3)
% plot(sig(3,:))
% subplot(4,1,4)
% plot(sig(4,:))

[ic, Out2, Out3] = fastica(mixedsig, 'numOfIC', 5);
figure
subplot(5,1,1)
plot(ic(1,:))
hold on 
title('Estimated independent components')
hold off
subplot(5,1,2)
plot(ic(2,:))
% subplot(5,1,3)
% plot(ic(3,:))
% subplot(5,1,4)
% plot(ic(4,:))
% subplot(5,1,5)
% plot(ic(5,:))
figure
scatter(ic(1,:)',ic(2,:)')

