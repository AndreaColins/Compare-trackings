function comparingPCAoverlapFreq
%%%%separar loss vectores en pedazos de 100 ms y 1 hz(?)

%%%%%%get different conditions 
pole=char('Toothpick','open coil','closed coil','black sandpaper','carbon pole','smooth pole','Bamboo');
%pole=char('air','carbon pole');
pole2=pole;
%pole2=char('Air','Pole');
color={'g','r','c','k','p'};
npole=size(pole,1);
nfiles=zeros(npole,1);
ms=5;%number of miliseconds
for i=1:npole
%%%%get the variables of every file
filename=strcat('./tr4s/tr4s/',pole(i,:),'/power',pole(i,:),'.mat');
v=load(filename,'-mat');
result(:,1,:)=v.azimuth;
result(:,2,:)=v.elevation;
result(:,3,:)=v.kcoronal;
result(:,4,:)=v.khorizontal;
nresult(i)=floor(size(result,1)/ms)
nfiles(i)=size(result,3)

result((nresult(i)*ms+1):end,:,:)=[];
%result(3451:end,:,:)=[];
%result(1:50,:,:)=[];
overlap=0;%miliseconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%azimuth
poleaz=reshape(squeeze(result(:,1,:)),nfiles(i).*nresult(i).*ms,1);

poleaz=poleaz./std(poleaz);
poleaz=reshape(poleaz,ms,nresult(i).*nfiles(i));

if overlap==0
    poleaz2=poleaz;
else
for j=1:nresult(i)*nfiles(i)-1
    poleaz2(:,j)=[poleaz(:,j);poleaz(1:overlap,j+1)];
end
%polekhor2(:,34*nfiles(i))=[polekhor(:,34*nfiles(i)); v.result(3401:3400+overlap,4,i)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%elevation
poleel=reshape(squeeze(result(:,2,:)),nfiles(i)*nresult(i)*ms,1);
poleel=poleel./std(poleel);
poleel=reshape(poleel,ms,nresult(i)*nfiles(i));

if overlap==0
    poleel2=poleel;
else
for j=1:nresult(i)*nfiles(i)-1
    poleel2(:,j)=[poleel(:,j);poleel(1:overlap,j+1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%kcoronal
polekcor=reshape(squeeze(result(:,3,:)),nfiles(i)*nresult(i)*ms,1);
polekcor=polekcor./std(polekcor);
polekcor=reshape(polekcor,ms,nresult(i)*nfiles(i));

if overlap==0
    polekcor2=polekcor;
else
for j=1:nresult(i)*nfiles(i)-1
    polekcor2(:,j)=[polekcor(:,j);polekcor(1:overlap,j+1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%khorizontal
polekhor=reshape(squeeze(result(:,4,:)),nfiles(i)*nresult(i)*ms,1);
polekhor=polekhor./std(polekhor);
polekhor=reshape(polekhor,ms,nresult(i)*nfiles(i));

if overlap==0
    polekhor2=polekhor;
else
for j=1:nresult(i)*nfiles(i)-1
    polekhor2(:,j)=[polekhor(:,j);polekhor(1:overlap,j+1)];
end
end

if i==1
tazimuth=poleaz2;
televation=poleel2;
tkcoronal=polekcor2;
tkhorizontal=polekhor2;
else
tazimuth=[tazimuth,poleaz2];
televation=[televation,poleel2];
tkcoronal=[tkcoronal,polekcor2];
tkhorizontal=[tkhorizontal,polekhor2];
end

clear result
end
%%%%apply PCA 
%%
%%%%%%%%%%%%%%%%%%%Fig1
%%%%%%%%%%%%%%%%%%azimuth,elevation
azel=[tazimuth',televation'];
 [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(azel,'Algorithm','eig');
 

figure
%%%graficar en los nuevos vectores 
subplot(2,2,1)
scatter(SCORE(1:nresult(1)*nfiles(1),1),SCORE(1:nresult(1)*nfiles(1),2))
%scatter3(SCORE(1:nresult(1)*nfiles(1),1),SCORE(1:nresult(1)*nfiles(1),2),SCORE(1:nresult(1)*nfiles(1),3))
totalfiles=nfiles(1)*nresult(1);

for i=2:npole
size(SCORE)
hold on
scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))
%scatter3(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),3))

title('[Azimuth,elevation]')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nfiles(i)*nresult(i);  
end
legend(pole2,'Location','best','Orientation','horizontal')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot Firsts components
subplot(2,2,2)
 plot(COEFF(:,1),'r')
 hold on
  plot(COEFF(:,2),'b')
   % plot(COEFF(:,3),'g')
  legend('First Component','Second Component')
 %plot(COEFF(:,1),'b*')
 %azel(1,:)
 %plot(0:1:199,azel(1,:),'b')
 xlabel('Time [milisecond]')
 vexplained=EXPLAINED(1)+EXPLAINED(2);
 title(strcat('Variance explained=',num2str(vexplained)))

 %plot(azel(101:200,1),'r*')
 %%%%%%%%%%%%%%%%%%%%%%%%explained variance
 subplot(2,2,3)
 plot(cumsum(EXPLAINED),'*-')
 xlabel('Number of principal component')
 ylabel('Variance explained [%]')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components
 
 subplot(2,2,4)
 scatter(COEFF(:,1),COEFF(:,2))
 xlabel('First Component')
 ylabel('Second Component')
 hold off
%%
%%%%%%%%%%%%%%%%%%%%%%%%%azimuth,khorizontal
%%%%%%%%%%%%%%%%%%%Fig2
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca([tazimuth',tkhorizontal'],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
subplot(2,2,1)
scatter(SCORE(1:nresult(1)*nfiles(1),1),SCORE(1:nresult(1)*nfiles(1),2))

totalfiles=nfiles(1)*nresult(1);

for i=2:npole
  
hold on
scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))
title('[Azimuth, Kappa Horizontal]')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nfiles(i)*nresult(i)
end
legend(pole2,'Location','best','Orientation','horizontal')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot firsts components
subplot(2,2,2)
 plot(COEFF(:,1),'r')
 hold on
  plot(COEFF(:,2),'b')
  legend('First Component','Second Component')
 %plot(COEFF(:,1),'b*')
 %azel(1,:)
 %plot(televation(:,1)','b')
 xlabel('Time [miliseconds]')
 vexplained=EXPLAINED(1)+EXPLAINED(2);
 title(strcat('Variance explained=',num2str(vexplained)))

 %plot(azel(101:200,1),'r*')
  %%%%%%%%%%%%%%%%%%%%%%%%explained variance
 subplot(2,2,3)
 plot(cumsum(EXPLAINED))
 xlabel('Number of principal component')
 ylabel('Variance explained [%]')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components
 
 subplot(2,2,4)
 scatter(COEFF(:,1),COEFF(:,2))
 xlabel('First Component')
 ylabel('Second Component')
 hold off
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%elevation, kappa
%%%%%%%%%%%%%%%%%%%Fig3
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca([televation',tkcoronal'],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
subplot(2,2,1)
scatter(SCORE(1:nresult(1)*nfiles(1),1),SCORE(1:nresult(1)*nfiles(1),2))

totalfiles=nfiles(1)*nresult(1); 

for i=2:npole
  
hold on
scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))

title('Elevation')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nfiles(i)*nresult(i);  
end
legend(pole2,'Location','best','Orientation','horizontal')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot firsts components
subplot(2,2,2)
 plot(COEFF(:,1),'r')
 hold on
  plot(COEFF(:,2),'b')
  legend('First Component','Second Component')
 %plot(COEFF(:,1),'b*')
 %azel(1,:)
 %plot(televation(:,1)','b')
 xlabel('Time [miliseconds]')
 vexplained=EXPLAINED(1)+EXPLAINED(2);
 title(strcat('Variance explained=',num2str(vexplained)))

 %plot(azel(101:200,1),'r*')
  %%%%%%%%%%%%%%%%%%%%%%%%explained variance
 subplot(2,2,3)
 plot(cumsum(EXPLAINED))
 xlabel('Number of principal component')
 ylabel('Variance explained [%]')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components
 
 subplot(2,2,4)
 scatter(COEFF(:,1),COEFF(:,2))
 xlabel('First Component')
 ylabel('Second Component')

 hold off
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%kappa coronal,kappa horizontal
%%%%%%%%%%%%%%%%%%%Fig4
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca([tkhorizontal', tkcoronal'],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
subplot(2,2,1)
scatter(SCORE(1:nresult(1)*nfiles(1),1),SCORE(1:nresult(1)*nfiles(1),2))

totalfiles=nfiles(1)*nresult(1);

for i=2:npole
  
hold on
scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))
title('[Kappa horizontal, Kappa coronal]')

totalfiles=totalfiles+nfiles(i)*nresult(i);
end
legend(pole2,'Location','best','Orientation','horizontal')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot firsts components
subplot(2,2,2)
 plot(COEFF(:,1),'r')
 hold on
  plot(COEFF(:,2),'b')
  legend('First Component','Second Component')
 %plot(COEFF(:,1),'b*')
 %azel(1,:)
 %plot(televation(:,1)','b')
 xlabel('Time [miliseconds]')
 vexplained=EXPLAINED(1)+EXPLAINED(2);
 title(strcat('Variance explained=',num2str(vexplained)))

 %plot(azel(101:200,1),'r*')
  %%%%%%%%%%%%%%%%%%%%%%%%explained variance
 subplot(2,2,3)
 plot(cumsum(EXPLAINED))
 xlabel('Number of principal component')
 ylabel('Variance explained [%]')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components
 
 subplot(2,2,4)
 scatter(COEFF(:,1),COEFF(:,2))
 xlabel('First Component')
 ylabel('Second Component')

 hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%azimuth, elevation, kappa horizontal, kappa coronal
%%%%%%%%%%%%%%%%%%%Fig5
[COEFF, SCORE]=pca([tazimuth',televation',tkhorizontal', tkcoronal'],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
subplot(2,2,1)
scatter(SCORE(1:nresult(1)*nfiles(1),1),SCORE(1:nresult(1)*nfiles(1),2))

totalfiles=nfiles(1)*nresult(1); 

for i=2:npole
  
hold on
scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))
title('[azimuth,elevation,k horizontal, k coronal]')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nfiles(i)*nresult(i);
end
legend(pole2,'Location','best','Orientation','horizontal')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot firsts components
subplot(2,2,2)
 plot(COEFF(:,1),'r')
 hold on
  plot(COEFF(:,2),'b')
  legend('First Component','Second Component')
 %plot(COEFF(:,1),'b*')
 %azel(1,:)
 %plot(televation(:,1)','b')
 xlabel('Time [miliseconds]')
 vexplained=EXPLAINED(1)+EXPLAINED(2);
 title(strcat('Variance explained=',num2str(vexplained)))

 %plot(azel(101:200,1),'r*')
  %%%%%%%%%%%%%%%%%%%%%%%%explained variance
 subplot(2,2,3)
 plot(cumsum(EXPLAINED))
 xlabel('Number of principal component')
 ylabel('Variance explained [%]')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components
 
 subplot(2,2,4)
 scatter(COEFF(:,1),COEFF(:,2))
 xlabel('First Component')
 ylabel('Second Component')
 
 hold off
end