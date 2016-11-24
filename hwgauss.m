function hw=hwgauss(pole)
filename=strcat('./tr4s/tr4s/',pole,'/autocov',pole,'.mat')
%to analize from frame 500
%%filename=strcat('./tr4s/tr4s/',pole,'/autocov',pole,'500.mat')
v=load(filename,'-mat');
nfiles=size(v.autoazimuth,2);
N0=round(size(v.autoazimuth(:,1),1)/2);
hw=zeros(nfiles,4);



xdata=1:1:size(v.autoazimuth(:,1),1);
p0 = [N0,100];
for i=1:nfiles
ydata=v.autoazimuth(:,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);
[p]= lsqcurvefit(fun,p0,xdata,ydata);

hw(i,1)=2.35482*p(2)/sqrt(2);
% ysim=exp(-((xdata-p(1))/p(2)).^2);
% plot(ydata,'bo')
% hold on 
% plot(ysim,'b')
% pause(1)
% hold off

ydata=v.autoelevation(:,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);
[p]= lsqcurvefit(fun,p0,xdata,ydata);
hw(i,2)=2.35482*p(2)/sqrt(2);

ydata=v.autokcoronal(:,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);
[p]= lsqcurvefit(fun,p0,xdata,ydata);
hw(i,3)=2.35482*p(2)/sqrt(2);

ydata=v.autokhorizontal(:,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);
[p]= lsqcurvefit(fun,p0,xdata,ydata);
hw(i,4)=2.35482*p(2)/sqrt(2);

end
end