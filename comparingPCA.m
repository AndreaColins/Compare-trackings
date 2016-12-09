function comparingPCA
%%%%separar loss vectores en pedazos de 100 ms y 1 hz(?)

%%%%%%get different conditions 
pole=char('air','smooth pole','closed coil','black sandpaper');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%azimuth
poleaz=reshape(squeeze(result(:,1,:)),nfiles(i)*3400,1);
poleaz=reshape(poleaz,100,34*nfiles(i));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%elevation
poleel=reshape(squeeze(result(:,2,:)),nfiles(i)*3400,1);
poleel=reshape(poleel,100,34*nfiles(i));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%kcoronal
polekcor=reshape(squeeze(result(:,3,:)),nfiles(i)*3400,1);
polekcor=reshape(polekcor,100,34*nfiles(i));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%khorizontal
polekhor=reshape(squeeze(result(:,4,:)),nfiles(i)*3400,1);
polekhor=reshape(polekhor,100,34*nfiles(i));

if i==1
tazimuth=poleaz;
televation=poleel;
tkcoronal=polekcor;
tkhorizontal=polekhor;
else
tazimuth=[tazimuth,poleaz];
televation=[televation,poleel];
tkcoronal=[tkcoronal,polekcor];
tkhorizontal=[tkhorizontal,polekhor];
end


end
%%%%apply PCA 

%%%%%%%%%%%%%%%%%%azimuth
 [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(tazimuth','Algorithm','eig');

figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);
color={'g','r','c','k','p'};
for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})
title('Azimuth')

totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%elevation
[COEFF, SCORE]=pca(televation','Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);
color={'g','r','c','k','p'};
for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})

title('Elevation')

totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%kappa coronal
[COEFF, SCORE]=pca(tkcoronal','Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);
color={'g','r','c','k','p'};
for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})

title('Kappa Coronal')

totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%kappa coronal
[COEFF, SCORE]=pca(tkhorizontal','Algorithm','eig');
figure
%%%graficar en los nuevos vectores 
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
scatter(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),'b')
totalfiles=nfiles(1);
color={'g','r','c','k','p'};
for i=2:npole
  
hold on
%scatter3(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,3),color{i-1})
scatter(SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,1),SCORE(totalfiles*34+1:(totalfiles+nfiles(i))*34,2),color{i-1})

title('Kappa Horizontal')

totalfiles=totalfiles+nfiles(i);  
end
legend(pole)
hold off
end