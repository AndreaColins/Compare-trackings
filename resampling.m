function resampling
pole=char('air');
pole2=pole;
colors={'b','r','g','k','c','m','y'};
npole=size(pole,1);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
resultair=[];
azimuth=[];
elevation=[];
kcoronal=[];
khorizontal=[];
amplitude=[];
amplitudeair=[];
fs=1;
if strcmp(pole(1,1:3),'air')
        filename=strcat('./tr4s/processed/FrameCorrected/result',pole(1,:),'.mat');
        v=load(filename,'-mat');
        result=v.result;
        figure
        ax1=subplot(3,1,1);
        plot(result(:,1,1),'+')
        hold on 
        y = resample(result(:,1,1),500,1000);
        size(y)
        t2=(0:(size(y,1)-1))*1000/(500*fs);
        stem(t2,y,'o')
        legend('Raw Azimuth','Resampling 500 [Hz]')
        xlabel('Time [ms]')
        ylabel('Angle (\circ)')
        hold off
        
        ax2=subplot(3,1,2);
        plot(result(:,1,1),'+')
        hold on 
        y = resample(result(:,1,1),250,1000);
        size(y)
        t2=(0:(size(y,1)-1))*1000/(250*fs);
        stem(t2,y,'o')
        legend('Raw Azimuth','Resampling 250 [Hz]')
        xlabel('Time [ms]')
        ylabel('Angle (\circ)')
        hold off
        
        ax3=subplot(3,1,3);
        plot(result(:,1,1),'+')
        hold on 
        y = resample(result(:,1,1),100,1000);
        size(y)
        t2=(0:(size(y,1)-1))*1000/(100*fs);
        stem(t2,y,'o')
        legend('Raw Azimuth','Resampling 100 [Hz]')
        xlabel('Time [ms]')
        ylabel('Angle (\circ)')
        hold off
        
    end

linkaxes([ax1,ax2,ax3],'xy')

end
