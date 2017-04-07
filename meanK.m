[idx,result]=comparetr4('C:\Users\mqbpwac5\Documents\MATLAB\videos\2016_08_24\TT3_20160824_120039.tr4',0);
meank=mean(result(1:500,3:4))
coronal=result(:,3)-meank(1);
horizontal=result(:,4)-meank(2);
plot(horizontal)

limith1=-(-0.001-meank(2));
limith2=(-0.001-meank(2));
hold on 
plot(limith1,'r')
plot(limith2,'g')
hold off

I=find(horizontal<limith1 | horizontal>limith2);

 figure
 subplot(2,1,1)
 plot(horizontal(I))
 subplot(2,1,2)
 plot(coronal(I))