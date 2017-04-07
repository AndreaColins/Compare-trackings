function [slope1,slope2,rsquared]=powerslope2(pole)
%%%%%%%%%%divide power spectrum in 2 sections and fit to 2 straight lines
%%%%%%%%%%%%First section: 0-20 Hz
%%%%%%%%%%%%Second section:20-500 Hz
%%%sacar los datos de la estructura 
filename=strcat('./tr4s/tr4s/',pole,'/power',pole,'.mat')
power=load(filename,'-mat');
nfiles=size(power.azimuth,2);
%%%hacer un vector que contenga las repeticiones 

xdata=log(power.freq(2:end));
xdata1=xdata(1:67);
xdata2=xdata(68:end);
slope1=zeros(nfiles,4);
slope2=slope1;
rsquared=slope1;
p0=[-1 0];

for i=1:nfiles
%%%%%%%%%%%%%%%%%%%%azimuth fit    
ydata=log(power.azimuth(2:end,i));
ydata1=ydata(1:67);
ydata2=ydata(68:end);

fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun1 = @(p,xdata)(xdata1*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
slope1(i,1)=p1(1);
slope2(i,1)=p2(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
rsquared(i,1)=corr([ysim1';ysim2'],ydata)^2;

plot(xdata,ydata)
hold on 
plot(xdata1,ysim1)
plot(xdata2,ysim2)
hold off
pause(1)


%%%%%%%%%%%%%%%%%%%%elevation fit   
ydata=log(power.elevation(2:end,i));
ydata1=ydata(1:67);
ydata2=ydata(68:end);

fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun1 = @(p,xdata)(xdata1*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
slope1(i,2)=p1(1);
slope2(i,2)=p2(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
rsquared(i,2)=corr([ysim1';ysim2'],ydata)^2;

% plot(xdata,ydata)
% hold on 
% plot(xdata1,ysim1)
% plot(xdata2,ysim2)
% hold off
% pause(0.5)


%%%%%%%%%%%%%%%%%%%%kappa coronal
ydata=log(power.kcoronal(2:end,i));
ydata1=ydata(1:67);
ydata2=ydata(68:end);

fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun1 = @(p,xdata)(xdata1*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
slope1(i,3)=p1(1);
slope2(i,3)=p2(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
rsquared(i,3)=corr([ysim1';ysim2'],ydata)^2;

% plot(xdata,ydata)
% hold on 
% plot(xdata1,ysim1)
% plot(xdata2,ysim2)
% hold off
% pause(1)

%%%%%%%%%%%%%%%%%%%%kappa horizontal
ydata=log(power.khorizontal(2:end,i));
ydata1=ydata(1:67);
ydata2=ydata(68:end);

fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun1 = @(p,xdata)(xdata1*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
slope1(i,4)=p1(1);
slope2(i,4)=p2(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
rsquared(i,4)=corr([ysim1';ysim2'],ydata)^2;

% plot(xdata,ydata)
% hold on 
% plot(xdata1,ysim1)
% plot(xdata2,ysim2)
% hold off
% pause(1)
end
rsquared=mean(rsquared);

end