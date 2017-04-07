function [slope1,slope2,slope3,rsquared]=powerslope3(pole)
%%%%%%%%%%divide power spectrum in 3 sections and fit to 3 straight lines
%%%%%%%%%%%%First section: 0-5 Hz
%%%%%%%%%%%%Second section:5-20 Hz
%%%%%%%%%%%%Third section:20-500 Hz
%%%sacar los datos de la estructura 

%%%%%tr4s
%filename=strcat('./tr4s/tr4s/',pole,'/power',pole,'.mat');

%%%%selected videos
filename=strcat('./videosselectedtextures/processed/FrameCorrected/power',pole,'.mat')
power=load(filename,'-mat');
nfiles=size(power.azimuth,2);
%%%hacer un vector que contenga las repeticiones 

xdata=log(power.freq(2:end))

xdata1=xdata(1:16);
xdata2=xdata(17:70);
xdata3=xdata(71:end);
slope1=zeros(nfiles,4);
slope2=slope1;
rsquared=slope1;
p0=[-1 0 1];

for i=1:nfiles
%%%%%%%%%%%%%%%%%%%%azimuth fit    
ydata=log(power.azimuth(2:end,i));
ydata1=ydata(1:16);
ydata2=ydata(17:70);
ydata3=ydata(71:end);


fun1 = @(p,xdata)(xdata1*p(1)+p(2));
fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun3 = @(p,xdata)(xdata3*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
[p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
slope1(i,1)=p1(1);
slope2(i,1)=p2(1);
slope3(i,1)=p3(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
ysim3=p3(1)*xdata3+p3(2);

rsquared(i,1)=corr([ysim1';ysim2';ysim3'],ydata)^2;

% plot(xdata,ydata)
% hold on 
% plot(xdata1,ysim1)
% plot(xdata2,ysim2)
% plot(xdata3,ysim3)
% hold off
% pause(1)


%%%%%%%%%%%%%%%%%%%%elevation fit   
ydata=log(power.elevation(2:end,i));
ydata1=ydata(1:16);
ydata2=ydata(17:70);
ydata3=ydata(71:end);


fun1 = @(p,xdata)(xdata1*p(1)+p(2));
fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun3 = @(p,xdata)(xdata3*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
[p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
slope1(i,2)=p1(1);
slope2(i,2)=p2(1);
slope3(i,2)=p3(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
ysim3=p3(1)*xdata3+p3(2);
rsquared(i,2)=corr([ysim1';ysim2';ysim3'],ydata)^2;

% plot(xdata,ydata)
% hold on 
% plot(xdata1,ysim1)
% plot(xdata2,ysim2)
% plot(xdata2,ysim2)
% hold off
% pause(1)

%%%%%%%%%%%%%%%%%%%%kappa coronal
ydata=log(power.kcoronal(2:end,i));
ydata1=ydata(1:16);
ydata2=ydata(17:70);
ydata3=ydata(71:end);


fun1 = @(p,xdata)(xdata1*p(1)+p(2));
fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun3 = @(p,xdata)(xdata3*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
[p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
slope1(i,3)=p1(1);
slope2(i,3)=p2(1);
slope3(i,3)=p3(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
ysim3=p3(1)*xdata3+p3(2);
rsquared(i,3)=corr([ysim1';ysim2';ysim3'],ydata)^2;

% plot(xdata,ydata)
% hold on 
% plot(xdata1,ysim1)
% plot(xdata2,ysim2)
% plot(xdata2,ysim2)
% hold off
% pause(1)

%%%%%%%%%%%%%%%%%%%%kappa horizontal
ydata=log(power.khorizontal(2:end,i));
ydata1=ydata(1:16);
ydata2=ydata(17:70);
ydata3=ydata(71:end);


fun1 = @(p,xdata)(xdata1*p(1)+p(2));
fun2 = @(p,xdata)(xdata2*p(1)+p(2));
fun3 = @(p,xdata)(xdata3*p(1)+p(2));

[p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
[p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
[p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
slope1(i,4)=p1(1);
slope2(i,4)=p2(1);
slope3(i,4)=p3(1);
ysim1=p1(1)*xdata1+p1(2);
ysim2=p2(1)*xdata2+p2(2);
ysim3=p3(1)*xdata3+p3(2);
rsquared(i,4)=corr([ysim1';ysim2';ysim3'],ydata)^2;

% plot(xdata,ydata)
% hold on 
% plot(xdata1,ysim1)
% plot(xdata2,ysim2)
% plot(xdata2,ysim2)
% hold off
% pause(1)
end
rsquared=mean(rsquared);

end