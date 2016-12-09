function [hw,rsquared]=hwgauss(pole)
filename=strcat('./tr4s/tr4s/',pole,'/autocov',pole,'.mat')
%to analize from frame 500
%%filename=strcat('./tr4s/tr4s/',pole,'/autocov',pole,'500.mat')
v=load(filename,'-mat');
nfiles=size(v.autoazimuth,2);
hw=zeros(nfiles,4);

tin=3000;
%tin=1;
tfin=4000;
%tfin=2*3488-1;
N0=(tfin-tin)/2+tin;
xdata=tin:1:tfin;
p0 = [N0,30];
for i=1:nfiles
ydata=v.autoazimuth(tin:tfin,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);
[p]= lsqcurvefit(fun,p0,xdata,ydata);

hw(i,1)=2.35482*p(2)/sqrt(2);
ysim=exp(-((xdata-p(1))/p(2)).^2);

rsquared(i,1)=corr(ysim',ydata')^2;
 plot(xdata,ydata,'bo')
 hold on 
 plot(xdata,ysim,'b')
 hold off
 pause(0.5)


ydata=v.autoelevation(tin:tfin,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);

[p]= lsqcurvefit(fun,p0,xdata,ydata);
hw(i,2)=2.35482*p(2)/sqrt(2);
ysim=exp(-((xdata-p(1))/p(2)).^2);
rsquared(i,2)=corr(ysim',ydata')^2;

%  plot(xdata,ydata,'bo')
%  hold on 
%  plot(xdata,ysim,'b')
%  hold off
%  pause(0.5)

ydata=v.autokcoronal(tin:tfin,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);
[p]= lsqcurvefit(fun,p0,xdata,ydata);
hw(i,3)=2.35482*p(2)/sqrt(2);
ysim=exp(-((xdata-p(1))/p(2)).^2);

rsquared(i,3)=corr(ysim',ydata')^2;
%  plot(xdata,ydata,'bo')
%  hold on 
%  plot(xdata,ysim,'b')
%  hold off
%  pause(0.5)

ydata=v.autokhorizontal(tin:tfin,i)';
fun = @(p,xdata)exp(-((xdata-p(1))/p(2)).^2);
[p]= lsqcurvefit(fun,p0,xdata,ydata);
hw(i,4)=2.35482*p(2)/sqrt(2);
ysim=exp(-((xdata-p(1))/p(2)).^2);
rsquared(i,4)=corr(ysim',ydata')^2;

%  plot(xdata,ydata,'bo')
%  hold on 
%  plot(xdata,ysim,'b')
%  hold off
%  pause(0.5)
end
end