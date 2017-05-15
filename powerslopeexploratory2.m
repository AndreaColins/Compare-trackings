function [slope1,slope2,slope3,ydata5]=powerslopeexploratory2(pole,color)
%%%%%%%%%%Analyse exploratory periods and divide power spectrum in 3 sections and fit to 3 straight lines
%%%%%%%%%%%%First section: 0-5 Hz
%%%%%%%%%%%%Second section:5-20 Hz
%%%%%%%%%%%%Third section:20-500 Hz
%%%sacar los datos de la estructura
ms=500;
ydata5=[];
nresult=0;
counter=1;
%power is the matrix with power spectrum information where every columns is
%a sample
%filename=strcat('./tr4s/processed/FrameCorrected/result',pole,'.mat');
filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole,'.mat');
if strcmp(pole,'air')
    
    v=load(filename,'-mat');
    result=v.result;
    idx=[1:3488];
    nfiles=size(result,3);
    startp=1;
    endp=3488;
    slope1=zeros(nfiles,4);
    slope2=slope1;
    slope3=slope1
    rsquared=slope1;

    for j=1:nfiles
        result=v.result;
        result=diff(diff(result(startp:endp,:,j),1,1),1,1);
        nresult=floor(size(result,1)/ms)
        result((nresult*ms+1):end,:,:)=[];
        for i=1:nresult
            power=powerspectra3(result(ms*(i-1)+1:ms*i,:));
            %%%find the points of x=5 and 20
            xdata=log(power.freq(2:end));
            x5=find(xdata<=log(5),1,'last');
            x20=find(xdata>=log(20),1,'first');
            
            xdata1=xdata(1:x5);
            xdata2=xdata(x5+1:x20);
            xdata3=xdata(x20+1:end);
            
            p0=[-1 -1 -1];
            
            
            %%%%%%%%%%%%%%%%%%%%azimuth fit
            ydata=log(power.azimuth(2:end));
            ydata1=ydata(1:x5);
            ydata2=ydata(x5+1:x20);
            ydata3=ydata(x20+1:end);
            
            fun1 = @(p,xdata)(xdata1*p(1)+p(2));
            fun2 = @(p,xdata)(xdata2*p(1)+p(2));
            fun3 = @(p,xdata)(xdata3*p(1)+p(2));
            
            [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
            [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
            [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
            
            slope1(counter,1)=p1(1);
            slope2(counter,1)=p2(1);
            slope3(counter,1)=p3(1);
            ysim1=p1(1)*xdata1+p1(2);
            ysim2=p2(1)*xdata2+p2(2);
            ysim3=p3(1)*xdata3+p3(2);
            
            rsquared(counter,1)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%elevation fit
            ydata=log(power.elevation(2:end));
            ydata1=ydata(1:x5);
            ydata2=ydata(x5+1:x20);
            ydata3=ydata(x20+1:end);
            
            fun1 = @(p,xdata)(xdata1*p(1)+p(2));
            fun2 = @(p,xdata)(xdata2*p(1)+p(2));
            fun3 = @(p,xdata)(xdata3*p(1)+p(2));
            
            [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
            [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
            [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
            
            slope1(counter,2)=p1(1);
            slope2(counter,2)=p2(1);
            slope3(counter,2)=p3(1);
            ysim1=p1(1)*xdata1+p1(2);
            ysim2=p2(1)*xdata2+p2(2);
            ysim3=p3(1)*xdata3+p3(2);
            
            rsquared(counter,2)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
            
            
            %%%%%%%%%%%%%%%%%%%%kcoronal fit
            ydata=log(power.kcoronal(2:end));
            ydata1=ydata(1:x5);
            ydata2=ydata(x5+1:x20);
            ydata3=ydata(x20+1:end);
            
            fun1 = @(p,xdata)(xdata1*p(1)+p(2));
            fun2 = @(p,xdata)(xdata2*p(1)+p(2));
            fun3 = @(p,xdata)(xdata3*p(1)+p(2));
            
            [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
            [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
            [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
            
            slope1(counter,3)=p1(1);
            slope2(counter,3)=p2(1);
            slope3(counter,3)=p3(1);
            ysim1=p1(1)*xdata1+p1(2);
            ysim2=p2(1)*xdata2+p2(2);
            ysim3=p3(1)*xdata3+p3(2);
            
            rsquared(counter,3)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
            
            
            %%%%%%%%%%%%%%%%%%%%khorizontal fit
            ydata=log(power.khorizontal(2:end));
            ydata1=ydata(1:x5);
            ydata2=ydata(x5+1:x20);
            ydata3=ydata(x20+1:end);
            
            fun1 = @(p,xdata)(xdata1*p(1)+p(2));
            fun2 = @(p,xdata)(xdata2*p(1)+p(2));
            fun3 = @(p,xdata)(xdata3*p(1)+p(2));
            
            [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
            [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
            [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
            
            slope1(counter,4)=p1(1);
            slope2(counter,4)=p2(1);
            slope3(counter,4)=p3(1);
            ysim1=p1(1)*xdata1+p1(2);
            ysim2=p2(1)*xdata2+p2(2);
            ysim3=p3(1)*xdata3+p3(2);
            
            rsquared(counter,4)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
            plot(exp(xdata),exp(ydata),'Color',color)
            hold on
            counter=counter+1;
        end
        
    end
    rsquared=mean(rsquared);
else
    
    %%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    
    
    d=strcat('./videosselectedtextures/processed/touchdetection/',pole,'/');
    ff = dir([d '*.mat']);
    nfiles=size(ff,1);
    string_list=cell(nfiles,1);
    slope1=zeros(nfiles,4);
    slope2=slope1;
    rsquared=slope1;
    for j=1:nfiles
        string_list{j}=strcat(d,ff(j).name);
        touchfile=strcat(d,ff(j).name)
        v=load(filename,'-mat');
        result=v.result;
        vtouch=load(touchfile,'-mat');
        idx=[1:3488];
        idx2=circshift(idx,size(idx,1)-vtouch.start_frame-1,2);
        result=diff(diff(result(idx,:,j),1,1),1,1);
        vtouch.touches=vtouch.touches(idx2);
        %%%%%%%select touch periods
        startp=find(vtouch.touches,1,'first');
        endp=find(vtouch.touches,1,'last');
        
        %%%%%%%%%%%%%%%for exploratory periods
        %[result,touchperiods]=centredk(endp,vtouch.touches,result);
        if (3488-endp)>=50
            endp=endp+50;
        else
            endp=3486;
        end
        
        
        if endp-startp>=ms
            result=result(startp:endp,:);
            nresult=floor(size(result,1)/ms);
            result((nresult*ms+1):end,:,:)=[];
            for i=1:nresult
                power=powerspectra3(result(ms*(i-1)+1:ms*i,:));
                %%%find the points of x=5 and 20
                xdata=log(power.freq(2:end));
                x5=find(xdata<=log(5),1,'last');
                x20=find(xdata>=log(20),1,'first');
                
                xdata1=xdata(1:x5);
                xdata2=xdata(x5+1:x20);
                xdata3=xdata(x20+1:end);
                
                p0=[-1 -1 -1];
                
                
                %%%%%%%%%%%%%%%%%%%%azimuth fit
                ydata=log(power.azimuth(2:end));
                ydata1=ydata(1:x5);
                ydata2=ydata(x5+1:x20);
                ydata3=ydata(x20+1:end);
                
                fun1 = @(p,xdata)(xdata1*p(1)+p(2));
                fun2 = @(p,xdata)(xdata2*p(1)+p(2));
                fun3 = @(p,xdata)(xdata3*p(1)+p(2));
                
                [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
                [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
                [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
                
                slope1(counter,1)=p1(1);
                slope2(counter,1)=p2(1);
                slope3(counter,1)=p3(1);
                ysim1=p1(1)*xdata1+p1(2);
                ysim2=p2(1)*xdata2+p2(2);
                ysim3=p3(1)*xdata3+p3(2);
                
                rsquared(counter,1)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
                
%                 plot(xdata,ydata,'Color',color)
%                 hold on
                %             plot(xdata1,ysim1)
                %             plot(xdata2,ysim2)
                %             plot(xdata3,ysim3)
                %             hold off
                %             pause(5)
                
                %%%%%%%%%%%%%%%%%%%%elevation fit
                ydata=log(power.elevation(2:end));
                ydata1=ydata(1:x5);
                ydata2=ydata(x5+1:x20);
                ydata3=ydata(x20+1:end);
                
                fun1 = @(p,xdata)(xdata1*p(1)+p(2));
                fun2 = @(p,xdata)(xdata2*p(1)+p(2));
                fun3 = @(p,xdata)(xdata3*p(1)+p(2));
                
                [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
                [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
                [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
                
                slope1(counter,2)=p1(1);
                slope2(counter,2)=p2(1);
                slope3(counter,2)=p3(1);
                ysim1=p1(1)*xdata1+p1(2);
                ysim2=p2(1)*xdata2+p2(2);
                ysim3=p3(1)*xdata3+p3(2);
                
                rsquared(counter,2)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
                
                %             plot(xdata,ydata,'b')
                %             hold on
                %             plot(xdata1,ysim1)
                %             plot(xdata2,ysim2)
                %             plot(xdata3,ysim3)
                %             hold off
                %             pause(5)
                %%%%%%%%%%%%%%%%%%%%kcoronal fit
                ydata=log(power.kcoronal(2:end));
                ydata1=ydata(1:x5);
                ydata2=ydata(x5+1:x20);
                ydata3=ydata(x20+1:end);
                
                fun1 = @(p,xdata)(xdata1*p(1)+p(2));
                fun2 = @(p,xdata)(xdata2*p(1)+p(2));
                fun3 = @(p,xdata)(xdata3*p(1)+p(2));
                
                [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
                [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
                [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
                
                slope1(counter,3)=p1(1);
                slope2(counter,3)=p2(1);
                slope3(counter,3)=p3(1);
                ysim1=p1(1)*xdata1+p1(2);
                ysim2=p2(1)*xdata2+p2(2);
                ysim3=p3(1)*xdata3+p3(2);
                
                rsquared(counter,3)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
                
                %             plot(xdata,ydata,'b')
                %             hold on
                %             plot(xdata1,ysim1)
                %             plot(xdata2,ysim2)
                %             plot(xdata3,ysim3)
                %             hold off
                %             pause(5)
                %%%%%%%%%%%%%%%%%%%%khorizontal fit
                ydata=log(power.khorizontal(2:end));
                ydata1=ydata(1:x5);
                ydata2=ydata(x5+1:x20);
                ydata3=ydata(x20+1:end);
                
                fun1 = @(p,xdata)(xdata1*p(1)+p(2));
                fun2 = @(p,xdata)(xdata2*p(1)+p(2));
                fun3 = @(p,xdata)(xdata3*p(1)+p(2));
                
                [p1]= lsqcurvefit(fun1,p0,xdata1,ydata1');
                [p2]= lsqcurvefit(fun2,p0,xdata2,ydata2');
                [p3]= lsqcurvefit(fun3,p0,xdata3,ydata3');
                
                slope1(counter,4)=p1(1);
                slope2(counter,4)=p2(1);
                slope3(counter,4)=p3(1);
                ysim1=p1(1)*xdata1+p1(2);
                ysim2=p2(1)*xdata2+p2(2);
                ysim3=p3(1)*xdata3+p3(2);
                
                rsquared(counter,4)=corr([ysim1';ysim2';ysim3'],[ydata1;ydata2;ydata3])^2;
                 plot(exp(xdata),exp(ydata),'Color',color)
                hold on
                %             plot(xdata,ydata,'b')
                %             hold on
                %             plot(xdata1,ysim1)
                %             plot(xdata2,ysim2)
                %             plot(xdata3,ysim3)
                %             hold off
                %             pause(5)
                counter=counter+1;
            end
        end
    end
    rsquared=mean(rsquared);
end
end


function power=powerspectra3(result)
%%%%%%%%%using multitapper method
%%%%%%%%%%%%meanplot==0 shows every curve in result
%%%%%%%%%%%meanplot/=0 shows the mean and standar deviation curve for every variable
N = size(result,1);
T = N/1000; %% define time of interval, 3.4 seconds
t = [0:N-1]/N;
t = t*T;
f = result;
freq =[0:N/2-1]/T;

for i=1:size(f,2)
    p(:,i)=pmtm(f(:,i),4,[],1000);
end
p = p(1:N/2,:); %% take the power of positve freq. half
power.freq=freq;
power.azimuth=p(:,1);
power.elevation=p(:,2);
power.kcoronal=p(:,3);
power.khorizontal=p(:,4);
% if meanplot==0
%     figure
%  h1(1)=subplot(4,1,1);
%
% %
% semilogy(freq,power.azimuth,'b');
% %pmtm(squeeze(f(:,1,:)))
% h1(2)=subplot(4,1,2);
% semilogy(freq,power.elevation,'g');
% h1(3)=subplot(4,1,3);
% semilogy(freq,power.kcoronal,'r');
% h1(4)=subplot(4,1,4);
% semilogy(freq,power.khorizontal,'c');
% else
%     figure
%     h1(1)=subplot(4,1,1)
% shadedErrorBar(freq,mean(log(squeeze(p(:,1,:))),2),std(log(squeeze(p(:,1,:))),0,2),'b');
% hold on
% xlabel('Frequency [Hz]')
% title('Azimuth')
% hold off
% h1(2)=subplot(4,1,2);
% shadedErrorBar(freq,mean(log(squeeze(p(:,2,:))),2),std(log(squeeze(p(:,2,:))),0,2),'g');
% hold on
% xlabel('Frequency [Hz]')
% title('Elevation')
% hold off
% h1(3)=subplot(4,1,3);
% shadedErrorBar(freq,mean(log(squeeze(p(:,3,:))),2),std(log(squeeze(p(:,3,:))),0,2),'r');
% hold on
% xlabel('Frequency [Hz]')
% title('Kappa Coronal')
% hold off
% h1(4)=subplot(4,1,4);
% shadedErrorBar(freq,mean(log(squeeze(p(:,4,:))),2),std(log(squeeze(p(:,4,:))),0,2),'c');
% hold on
% xlabel('Frequency [Hz]')
% title('Kappa Horizontal')
% hold off
% end
%
%     linkaxes(h1,'x')

end
function [result,touchperiods]=centredk(endp,touches,result)
i=1;
touchescopy=touches;
touchcount=0;
while i<endp
    ini=find(touchescopy,1,'first');
    touchescopy(1:ini)=1;
    fin=find(~touchescopy,1,'first')-1;
    if isempty(fin)
        fin=3488;
    end
    touchescopy(1:fin)=0;
    hold on
    touchperiods(touchcount+1,1:2)=[ini,fin];
    touchcount=touchcount+1;
    i=fin;
end
%touchperiods

kappacoronal=result(:,3,1);
kappahorizontal=result(:,4,1);
kappacoronal(~touches)=0;
kappahorizontal(~touches)=0;
for i=1:touchcount
    if i==1
        if touchperiods(1,1)>10
            k0coronal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,4));
        else
            k0coronal=median(result(1:touchperiods(i,1)-1,3));
            k0horizontal=median(result(1:touchperiods(i,1)-1,4));
        end
    else
        if touchperiods(i,1)-touchperiods(i-1,2)>10
            k0coronal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,4));
        else
            k0coronal=median(result(touchperiods(i-1,2)+1:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i-1,2)+1:touchperiods(i,1)-1,4));
        end
    end
    kappacoronal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),3)-k0coronal;
    kappahorizontal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),4)-k0horizontal;
    
end
result(:,3)=kappacoronal;
result(:,4)=kappahorizontal;
clear touchescopy ini fin
end
function [result]=touchselection(touches,result)
% plot(result(:,4))
% hold on
I=find(touches);
result=result(I,:);
% plot(result(:,4),'.r')
% pause(5)
% hold off
clear I
end