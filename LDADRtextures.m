function [classp,azel]=LDADRtextures
%%%%this version of LDA analyse data in exploratory/touch periods, i.e, and
%%%%is used to dimensionality reduction in times series may have different lenghts

pole=char( 'Smooth pole', 'Carbon Pole', 'Black Sandpaper','Closed coil','Open coil','Cardboard','Bamboo','Toothpick','Wood');
%pole=char('black sandpaper','smooth pole');
pole2=pole;
npole=size(pole,1);

nfiles=zeros(npole,1);
nresult=zeros(npole,1);
ms=20;%number of miliseconds
overlap=10;%miliseconds
poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
resultp=[];
touchlength=[];
for i=1:npole
    %%%%%%%%%%%%%%%raw data
    filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    %filename=strcat('./videosselectedtextures/processed/bandpass/from20/result',pole(i,:),'.mat');
        d=strcat('./videosselectedtextures/processed/touchdetection/',pole(i,:),'/');
        ff = dir([d '*.mat']);
        nfiles(i)=size(ff,1);
        string_list=cell(nfiles(i),1);
        for j=1:nfiles(i)
            string_list{j}=strcat(d,ff(j).name);
            touchfile=strcat(d,ff(j).name)
            v=load(filename,'-mat');
            result=v.result;
            vtouch=load(touchfile,'-mat');
            idx=[1:3488];
            idx2=circshift(idx,size(idx,1)-vtouch.start_frame-1,2);
            result=result(idx,:,j);
            vtouch.touches=vtouch.touches(idx2);
            %%%%%%%select touch periods
            startp=find(vtouch.touches,1,'first');
            endp=find(vtouch.touches,1,'last');

            %%%%%%%%%%%%%%%for exploratory periods
            [result2,touchp]=centredk(endp,vtouch.touches,result);
            if (3488-endp)>=50
                endp=endp+50;
            else
                endp=3488;
            end
            touchp=[startp endp];
            %result=result2(startp:endp,:);
%             
            %%%%%%%%%%%%for touch periods
             %[result2,touchp]=centredk(endp,vtouch.touches,result);
             %discard small touch periods
             counter=1;
             touchperiods=[];
             for t=1:size(touchp,1)
                 if (touchp(t,2)-touchp(t,1))>=ms
                     touchperiods(counter,:)=touchp(t,:);
                     counter=counter+1;
                 end
             end
           
             clear touchp 
             if ~isempty(touchperiods)
            for t=1:size(touchperiods,1)
                result=result2(touchperiods(t,1):touchperiods(t,2),:);

            nsamples=floor(size(result,1)/(ms-overlap));
            if mod(size(result,1),(ms-overlap))<overlap
                nsamples=nsamples-1;
            end
            nresult(i)=nsamples+nresult(i);
            result((nsamples*(ms-overlap)+overlap+1):end,:)=[];
            
             if overlap==0

                poleaz2=reshape(result(:,1),ms,nsamples);
                poleel2=reshape(result(:,2),ms,nsamples);
                polekcor2=reshape(result(:,3),ms,nsamples);
                polekhor2=reshape(result(:,4),ms,nsamples);
             else
                
                for j=0:nsamples-1
                    poleaz2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,1);
%                     plot(j*(ms-overlap)+1:1:(j+1)*ms-j*overlap,poleaz2(:,j+1))
%                     pause(2)
                    poleel2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,2);
                    polekcor2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,3);
                    polekhor2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,4);
                end
                
            end
    
            poleaz=[poleaz, poleaz2];

            poleel=[poleel,poleel2];
            
            polekcor=[polekcor,polekcor2];
            
            polekhor=[polekhor, polekhor2];
            
            clear poleaz2 poleel2 polekcor2 polekhor2 result 
            end
            end
            clear touchperiods
        end 
        
    
        
    
    %%%%%%%%%%%%%%%%%%Add new condition
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
    
    poleaz=[];
    poleel=[];
    polekcor=[];
    polekhor=[];
    resultp=[];
    clear result
end
%%
%%%[Azimuth,Elevation]
azel=zscore([tazimuth',televation']);

%azel=[tazimuth',televation'];
%%%%%%%%%Raw Data
size(azel)
nresult
plot(azel(1,:)','b')
hold on 
plot(azel(nresult(1)+1,:)','r')
%plot(azel(nresult(1)+nresult(2)+1,:)','g')
plot(azel(1:nresult(1),:)','b')
plot(azel(nresult(1)+1:end,:)','r')
%plot(azel(nresult(1)+nresult(2)+1:end,:)','g')
hold off
%%
%     totalfiles=nresult(1)
%     poleclass(1:totalfiles)={pole(1,:)};
% 
%     for i=2:npole
%         poleclass(totalfiles+1:totalfiles+nresult(i))={pole(i,:)};
%         totalfiles=totalfiles+nresult(i);
%     end
%%%%Apply PCA angles
%MdlLinear = fitcdiscr(azel,poleclass);
classp=LDA(azel,'[Azimuth,Elevation]',pole,nresult);

%%curvatures
%LDA(zscore(tkhorizontal'),'[Kappa Horizontal]',pole,nresult)

%LDA(tkcoronal','[Kappa Coronal]',pole,nresult)

LDA([tkhorizontal',tkcoronal'],'[Kappa Horizontal, Kappa Coronal]',pole,nresult);

LDA([zscore(azel),zscore([tkhorizontal',tkcoronal'])],'[Azimuth,Elevation,Kappa Horizontal, Kappa Coronal]',pole,nresult);

% figure
% edges=[0  200:1:2888];
% h1=histogram(touchlength,edges);
% hold on 
% 
% xlabel('Length of exploratory periods [ms]')
% ylabel('Number of exploratory periods')
% hold off
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
        elseif touchperiods(1,1)==1
            k0coronal=median(result(touchperiods(i+1,1)-10:touchperiods(i+1,1)-1,3));
            k0horizontal=median(result(touchperiods(i+1,1)-10:touchperiods(i+1,1)-1,4));
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
function [result]=notouchselection(touches,result)
% plot(result(:,4))
% hold on
I=find(touches,1,'first');
I2=find(touches,1,'last');
if I2==3488
result=result(1:I,:);
else
result=[result(1:I,:);result(I2:end,:)];
end
% plot(result(:,4),'.r')
% pause(5)
% hold off
clear I I2
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

function classp=LDA(data,titlevar,classes,nresult)
totalfiles=nresult;
means=ones(size(data,2),size(nresult,1));
S_classes=means;
means(:,1)=mean(data(1:nresult(1),:))';

S_classes=cov(data(1:nresult(1),:));
for i=2:size(nresult,1)
    means(:,i)=mean(data(totalfiles+1:nresult(i)+totalfiles,:))';
    S_classes=S_classes+cov(data(totalfiles+1:nresult(i)+totalfiles,:));
    totalfiles=totalfiles+nresult(i);
end   

mo=mean(data)';
S_w=S_classes;
%S_w1=inv(S_w);

%S_b=((m1-mo)*(m1-mo)'+(m2-mo)*(m2-mo)'+(m3-mo)*(m3-mo)')/3;
S_b=zeros(size(S_w));
for i=1:size(nresult,1)
    S_b=S_b+nresult(i).*(means(:,i)-mo)*(means(:,i)-mo)';
    
end
S_b=S_b./size(nresult,1);
A=S_w\S_b;
[V,D]=eig(A);
sample=A*V-V*D;
sample(1,1)

D=diag(D);
[Dsort,idx]=sort(D,'descend');

%sort the corresponding eigenvectors
Vsort=zeros(size(V));

for i=1:size(V,2)
    Vsort(:,i)=V(:,idx(i));
end

V=real(Vsort);

classp=V(:,1)'*data';
classp2=V(:,2)'*data';

colors=distinguishable_colors(size(nresult,1));
figure
subplot(1,3,1)
scatter(classp(1:nresult(1)),classp2(1:nresult(1)),[],'MarkerEdgeColor',colors(1,:))
totalfiles=nresult(1);
for i=2:size(nresult,1)
    hold on
scatter(classp(totalfiles+1:(totalfiles+nresult(i))),classp2(totalfiles+1:(totalfiles+nresult(i))),[],'MarkerEdgeColor',colors(i,:))   
%title('[Azimuth, Kappa Horizontal]')
xlabel('First LDA')
ylabel('Second LDA')
title(titlevar)
totalfiles=totalfiles+nresult(i);

end
legend(classes,'Location','best')%,'Orientation','horizontal') 
hold off

subplot(1,3,2)
Dsort=real(Dsort);
EXPLAINED=Dsort/sum(Dsort)*100;
plot(cumsum(EXPLAINED),'*-')
xlabel('Number of LDA Component')
ylabel('Eigenvalue [%]')

subplot(1,3,3)
plot(V(:,1),'-*')
hold on 
xlabel('Time [ms]')
ylabel('First LDA component value')
hold off
%subplot(2,2,4)
% [~,edges]=histcounts(data,200);
% hold on
% histogram(data(1:nresult(1),:),edges,'Normalization','pdf')
% totalfiles=nresult(1);
% hold on
% for i=2:size(nresult,1)
%     histogram(data(totalfiles+1:nresult(i)+totalfiles,:),edges,'Normalization','pdf');
%     totalfiles=totalfiles+nresult(i);
% end   
% title(titlevar)
% %xlabel('Elevation')
% legend(classes,'Location','best')
% hold off
% [~,edges]=histcounts([classp],200);
% hold on
% histogram(classp(1:nresult(1)),edges,'Normalization','pdf')
% totalfiles=nresult(1);
% hold on
% for i=2:size(nresult,1)
%     histogram(classp(totalfiles+1:nresult(i)+totalfiles),edges,'Normalization','pdf');
%     totalfiles=totalfiles+nresult(i);
% end   
% title(titlevar)
% xlabel('First LDA Component')
% legend(classes,'Location','best')
% hold off
end