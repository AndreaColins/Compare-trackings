
%%%%%%%%%%%smooth pole
 touchfile=strcat('./tr4s/tr4s/touchdetectionGUI/','TT3_20160812_133608','_touch.mat');
 %filename=strcat('./tr4s/processed/FrameCorrected/resultsmooth pole.mat');
 filename=strcat('./tr4s/processed/bandpass/upto5/resultsmooth pole.mat');
%%%%%%%%%%closed coil
%touchfile=strcat('./tr4s/tr4s/touchdetectionGUI/','TT3_20160819_143726','_touch.mat');
%filename=strcat('./tr4s/processed/RawData/resultclosed coil.mat');
%%%%%%%%%open coil
%touchfile=strcat('./tr4s/processed/touchdetection/Bamboo/','TT3_20160908_124030','_touch.mat');
%filename=strcat('./tr4s/processed/bandpass/upto5/resultsmooth pole.mat');
%filename=strcat('./tr4s/processed/FrameCorrected/resultBamboo.mat');
vtouch=load(touchfile,'-mat');
v=load(filename,'-mat');
%vraw=load(filename2,'-mat');
idx=[1:3488];
idx2=circshift(idx,size(idx,1)-vtouch.start_frame-1,2);
result=v.result(:,:,1);
%resultraw=vraw.result(idx2,:,1);
figure
plot(vtouch.touches(idx2))
I=find(vtouch.touches(idx2));
firstt=I(1);
lastt=I(end)
i=1;
touchescopy=vtouch.touches(idx2);
touchcount=0;
while i<lastt
    ini=find(touchescopy,1,'first')
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
touchperiods
I2=find(~vtouch.touches(idx2));
kappacoronal=result(:,3,1);
kappahorizontal=result(:,4,1);
kappacoronal(I2)=0;
kappahorizontal(I2)=0;
for i=1:touchcount
        kappacoronal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),3,1)-result(touchperiods(i,1)-1,3,1);
    kappahorizontal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),4,1)-result(touchperiods(i,1)-1,4,1);
end

figure
h1(1)=subplot(2,2,1);
plot(result(:,3,1))
hold on 
%plot(resultraw(:,3,1),'g')
plot(I,result(I,3,1),'r.')
title('Raw data')
ylabel('Kappa Coronal')
xlabel('Time[ms]')
legend('Raw Data','Touch periods (non centred k)')
hold off

h1(2)=subplot(2,2,2);
plot(result(:,3,1))
hold on 
plot(kappacoronal,'r')
title('Centred Curvature')
ylabel('Kappa Coronal')
xlabel('Time[ms]')
legend('Raw Data','Touch periods (centred k)')
hold off

h2(1)=subplot(2,2,3);
plot(result(:,4,1))
hold on 
%plot(resultraw(:,4,1),'g')
plot(I,result(I,4,1),'r.')
ylabel('Kappa Horizontal')
xlabel('Time[ms]')
legend('Raw Data','Touch periods (non centred k)')
hold off

h2(2)=subplot(2,2,4);
plot(result(:,4,1))
hold on
plot(kappahorizontal,'r')
ylabel('Kappa Horizontal')
xlabel('Time[ms]')
legend('Raw Data','Touch periods (centred k)')
hold off
linkaxes(h1,'xy');
linkaxes(h2,'xy');
clear touchperiods