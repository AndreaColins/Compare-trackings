function [slope,intercept,rsquared]=powersquared(pole)
%%%sacar los datos de la estructura 
filename=strcat('./tr4s/tr4s/',pole,'/power',pole,'.mat')
power=load(filename,'-mat');
nfiles=size(power.azimuth,2);
%%%hacer un vector que contenga las repeticiones 
xdata=power.freq(5:end);
slope=zeros(nfiles,4);
intercept=slope;
rsquared=slope;
p0=[10 5 0.5]

for i=1:nfiles
%%%%%%%%%%%%%%%%%%%%azimuth fit    
ydata=log(power.azimuth(5:end,i));
fun = @(p,xdata)(p(1)./xdata.^p(3)-p(2));
[p]= lsqcurvefit(fun,p0,xdata,ydata');
slope(i,1)=p(1);
intercept(i,1)=p(3);
ysim=p(1)./xdata.^p(3)-p(2);
rsquared(i,1)=corr(ysim',ydata)^2;

% plot(xdata,ydata,'*')
% hold on 
% plot(xdata,ysim,'r')
% hold off
% pause(1)
%
% %%%%%%%%%%%%%%%%%%%%elevation fit   
 ydata=log(power.elevation(5:end,i));
 fun = @(p,xdata)(p(1)./xdata.^p(3)-p(2));
 [p]= lsqcurvefit(fun,p0,xdata,ydata');
 slope(i,2)=p(1);
 intercept(i,2)=p(3);
 ysim=p(1)./xdata.^p(3)-p(2);
 rsquared(i,2)=corr(ysim',ydata)^2;
% 
% 
% %%%%%%%%%%%%%%%%%%%%kappa coronal
 ydata=log(power.kcoronal(5:end,i));
fun = @(p,xdata)(p(1)./xdata.^p(3)-p(2));
 [p]= lsqcurvefit(fun,p0,xdata,ydata');
 slope(i,3)=p(1);
 intercept(i,3)=p(3);
 ysim=p(1)./xdata.^p(3)-p(2);
 rsquared(i,3)=corr(ysim',ydata)^2;
 
 %%%%%%%%%%%%%%%%%%%%kappa horizontal
 ydata=log(power.khorizontal(5:end,i));
fun = @(p,xdata)(p(1)./xdata.^p(3)-p(2));

 [p]= lsqcurvefit(fun,p0,xdata,ydata');
 slope(i,4)=p(1);
intercept(i,4)=p(3);
ysim=p(1)./xdata.^p(3)-p(2);
rsquared(i,4)=corr(ysim',ydata)^2;

% plot(ydata,ysim)
% hold on 
% ysim=p(1)*xdata+p(2);
% plot(xdata,ysim)
% hold off
% pause(1)
% 
 end
 rsquared=mean(rsquared);
end