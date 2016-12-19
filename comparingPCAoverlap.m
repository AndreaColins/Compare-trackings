function comparingPCAoverlap
%%%%separar loss vectores en pedazos de 100 ms y 1 hz(?)

%%%%%%get different conditions 
pole=char('air','smooth pole','carbon pole','black sandpaper','closed coil','open coil');
color={'g','r','c','k','p'};
npole=size(pole,1);
nfiles=zeros(npole,1);
ms=100;%number of miliseconds
for i=1:npole
%%%%get the variables of every file
filename=strcat('./tr4s/tr4s/',pole(i,:),'/result',pole(i,:),'.mat');
v=load(filename,'-mat');
result=v.result;
nfiles(i)=size(result,3);

result(3401:end,:,:)=[];
%result(3451:end,:,:)=[];
%result(1:50,:,:)=[];
overlap=0;%miliseconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%azimuth
poleaz=reshape(squeeze(result(:,1,:)),nfiles(i)*3400,1);

poleaz=poleaz./std(poleaz);
poleaz=reshape(poleaz,ms,34*nfiles(i));

if overlap==0
    poleaz2=poleaz;
else
for j=1:34*nfiles(i)-1
    poleaz2(:,j)=[poleaz(:,j);poleaz(1:overlap,j+1)];
end
%polekhor2(:,34*nfiles(i))=[polekhor(:,34*nfiles(i)); v.result(3401:3400+overlap,4,i)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%elevation
poleel=reshape(squeeze(result(:,2,:)),nfiles(i)*3400,1);
poleel=poleel./std(poleel);
poleel=reshape(poleel,ms,34*nfiles(i));

if overlap==0
    poleel2=poleel;
else
for j=1:34*nfiles(i)-1
    poleel2(:,j)=[poleel(:,j);poleel(1:overlap,j+1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%kcoronal
polekcor=reshape(squeeze(result(:,3,:)),nfiles(i)*3400,1);
polekcor=polekcor./std(polekcor);
polekcor=reshape(polekcor,ms,34*nfiles(i));

if overlap==0
    polekcor2=polekcor;
else
for j=1:34*nfiles(i)-1
    polekcor2(:,j)=[polekcor(:,j);polekcor(1:overlap,j+1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%khorizontal
polekhor=reshape(squeeze(result(:,4,:)),nfiles(i)*3400,1);
polekhor=polekhor./std(polekhor);
polekhor=reshape(polekhor,ms,34*nfiles(i));

if overlap==0
    polekhor2=polekhor;
else
for j=1:34*nfiles(i)-1
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


end
%%%%apply PCA 
%%
%%%%%%%%%%%%%%%%%%azimuth
azel=[televation',tazimuth'];
 [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(azel,'Algorithm','eig');
 

 figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
subplot(2,2,1)
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);

for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})
title('[Azimuth,elevation]')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
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
%%%%%%%%%%%%%%%%%%%%%%%%%elevation
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(tazimuth','Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
subplot(2,2,1)
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);
%color={'g','r','c','k','p'};
for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})

title('Azimuth')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
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
%%%%%%%%%%%%%%%%%%%%%%%%%elevation
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(televation','Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
subplot(2,2,1)
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);

for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})

title('Elevation')
xlabel('First component')
ylabel('Second component')
totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
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
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca([tkhorizontal', tkcoronal'],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
subplot(2,2,1)
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);

for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})

title('[Kappa horizontal, Kappa coronal]')

totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%kappa coronal
[COEFF, SCORE]=pca([tazimuth',televation',tkhorizontal', tkcoronal'],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
subplot(2,2,1)
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);

for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})

title('[azimuth,elevation,k horizontal, k coronal]')

totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
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