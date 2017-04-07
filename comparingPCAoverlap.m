function comparingPCAoverlap
%%%%separar los vectores en pedazos de 100 ms y 1 hz(?)

%%%%%%get different conditions
pole=char('air','Bamboo','smooth pole', 'carbon pole', 'black sandpaper', 'closed coil','open coil');
%pole=char('smooth pole','air');
pole2=pole;
%pole2=char('Air','Pole');
npole=size(pole,1);
nfiles=zeros(npole,1);
ms=100;%number of miliseconds
overlap=0;%miliseconds
for i=1:npole
    %%%%get the variables of every file
    %%%%%%%%%%%%%%%%filtered data
    %filename=strcat('./tr4s/processed/bandpass/from20/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%raw data
    filename=strcat('./tr4s/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%touch periods1
    %filename=strcat('./tr4s/tr4s/touchdetection1/result',pole(i,:),'.mat');
    v=load(filename,'-mat');
    result=v.result;
    nresult(i)=floor(size(result,1)/ms);
    nfiles(i)=size(result,3);
    
    result((nresult(i)*ms+1):end,:,:)=[];
    %result(3451:end,:,:)=[];
    %result(1:50,:,:)=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%azimuth
    poleaz=reshape(squeeze(result(:,1,:)),nfiles(i).*nresult(i).*ms,1);
    poleaz=reshape(poleaz,ms,nresult(i).*nfiles(i));
     meanaz(i)=mean(mean(poleaz));
    if overlap==0
        poleaz2=poleaz;
    else
        for j=1:nresult(i)*nfiles(i)-1
            poleaz2(:,j)=[poleaz(:,j);poleaz(1:overlap,j+1)];
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%elevation
    poleel=reshape(squeeze(result(:,2,:)),nfiles(i)*nresult(i)*ms,1);
    poleel=reshape(poleel,ms,nresult(i)*nfiles(i));
    meanel(i)=mean(mean(poleel));
    if overlap==0
        poleel2=poleel;
    else
        for j=1:nresult(i)*nfiles(i)-1
            poleel2(:,j)=[poleel(:,j);poleel(1:overlap,j+1)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%kcoronal
    meank1=mean(mean(squeeze(result(1:10,3,:)),1));
    polekcor=reshape(squeeze(result(:,3,:)),nfiles(i)*nresult(i)*ms,1);
    polekcor=reshape(polekcor,ms,nresult(i)*nfiles(i));
    meankcor(i)=mean(mean(polekcor));
    if overlap==0
        polekcor2=polekcor;
    else
        for j=1:nresult(i)*nfiles(i)-1
            polekcor2(:,j)=[polekcor(:,j);polekcor(1:overlap,j+1)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%khorizontal
     meank2=mean(mean(squeeze(result(1:10,3,:)),1));
    polekhor=reshape(squeeze(result(:,4,:)),nfiles(i)*nresult(i)*ms,1);
    polekhor=reshape(polekhor,ms,nresult(i)*nfiles(i));
    meankhor(i)=mean(mean(polekhor));
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
    
 clear poleaz2 poleel2 polekcor2 polekhor2 
end
% tazimuth=tazimuth-mean(tazimuth(:));
% televation=televation-mean(televation(:));
% tkcoronal= tkcoronal-mean(tkcoronal(:));
% tkhorizontal=tkhorizontal-mean(tkhorizontal(:));
%%%%apply PCA
%%
%%%%%%%%%%%%%%%%%%%Fig1
%%%%%%%%%%%%%%%%%%azimuth,elevation
azel=[tazimuth',televation'];
size(azel)
 figure
 
% offsetangle=(meanel+meanaz)./2;
% bar(offsetangle)
% hold on
% set(gca,'XTickLabel',pole2,'FontSize',14)
% title('mean(Elevation,Azimuth)','FontSize',16)
% hold off


plot(azel(1,:)','b')
hold on 
plot(azel(nfiles(1)*nresult(1)+1,:)','r')
plot(azel(nfiles(1)*nresult(1)+nfiles(2)*nresult(2)+1,:)','g')
plot(azel(1:nfiles(1)*nresult(1),:)','b')
plot(azel(nfiles(1)*nresult(1)+1:nfiles(1)*nresult(1)+nfiles(2)*nresult(2),:)','r')
plot(azel(nfiles(1)*nresult(1)+nfiles(2)*nresult(2)+1:end,:)','g')
hold off
legend(pole2,'Location','best','Orientation','horizontal')
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(azel,'Algorithm','eig');
%[COEFF,EXPLAINED,SCORE,varcorr]=pcandrea(azel);
figure
%%%graficar en los nuevos vectores
subplot(2,2,1)
if overlap==0
    totalfiles=nfiles(1)*nresult(1);
else
    totalfiles=nfiles(1)*nresult(1)-1;
end
scatter(SCORE(1:totalfiles,1),SCORE(1:totalfiles,2),'b')

for i=2:npole
    
    hold on
    scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))   
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
legend('First Component','Second Component')
xlabel('Time [milisecond]')
vexplained=EXPLAINED(1)+EXPLAINED(2);
title(strcat('Variance explained=',num2str(vexplained)))


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
%%%%%%%%%%%%%%%%%%azimuth,khorizontal
%%%%%%%%%%%%%%%%%%Fig2
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca([(tazimuth'-mean(tazimuth(:)))./std(tazimuth(:)),(tkhorizontal'-mean(tkhorizontal(:)))./std(tkhorizontal(:))],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores
subplot(2,2,1)
if overlap==0
    totalfiles=nfiles(1)*nresult(1);
else
    totalfiles=nfiles(1)*nresult(1)-1;
end
scatter(SCORE(1:totalfiles,1),SCORE(1:totalfiles,2),'b')



for i=2:npole
    
    hold on
    scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))
    title('[Azimuth, Kappa Horizontal]')
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
%%%%%%%%%%%%%%%%%%%%%%%%%elevation, kappa
%%%%%%%%%%%%%%%%%%%Fig3
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca([tkhorizontal'],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores
subplot(2,2,1)
if overlap==0
    totalfiles=nfiles(1)*nresult(1);
else
    totalfiles=nfiles(1)*nresult(1)-1;
end
scatter(SCORE(1:totalfiles,1),SCORE(1:totalfiles,2),'b')

for i=2:npole
    
    hold on
    scatter(SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nfiles(i)*nresult(i)),2))
    
    title('[kappa horizontal')
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%kappa coronal,kappa horizontal
%%%%%%%%%%%%%%%%%%%Fig4
curvatures=[tkhorizontal',tkcoronal'];
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca([(curvatures-mean(curvatures(:)))./std(curvatures(:))],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
subplot(2,2,1)
if overlap==0
    totalfiles=nfiles(1)*nresult(1);
else
    totalfiles=nfiles(1)*nresult(1)-1;
end
scatter(SCORE(1:totalfiles,1),SCORE(1:totalfiles,2),'b')

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
xlabel('Time [miliseconds]')
vexplained=EXPLAINED(1)+EXPLAINED(2);
title(strcat('Variance explained=',num2str(vexplained)))

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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%azimuth, elevation, kappa horizontal, kappa coronal
%%%%%%%%%%%%%%%%%%%Fig5
[COEFF, SCORE]=pca([(azel-mean(azel(:)))./std(azel(:)),(curvatures-mean(curvatures(:)))./std(curvatures(:))],'Algorithm','eig');
figure
%%%graficar en los nuevos vectores
%scatter3(SCORE(1:34*nfiles(1),1),SCORE(1:34*nfiles(1),2),SCORE(1:34*nfiles(1),3),'b')
subplot(2,2,1)
if overlap==0
    totalfiles=nfiles(1)*nresult(1);
else
    totalfiles=nfiles(1)*nresult(1)-1;
end
scatter(SCORE(1:totalfiles,1),SCORE(1:totalfiles,2),'b')

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
plot(cumsum(EXPLAINED),'*-')
xlabel('Number of principal component')
ylabel('Variance explained [%]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components

subplot(2,2,4)
scatter(COEFF(:,1),COEFF(:,2))
xlabel('First Component')
ylabel('Second Component')

hold off
end