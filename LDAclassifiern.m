function Mdlel=LDAclassifiern
%%%%this version of LDA analyse data in exploratory/touch periods, i.e, and
%%%%is used to dimensionality reduction in times series may have different lenghts

%pole=char('smooth pole', 'carbon pole', 'black sandpaper', 'closed coil','open coil','Bamboo','Toothpick');
pole=char('smooth pole', 'carbon pole', 'black sandpaper', 'closed coil','open coil','Bamboo','Toothpick');
pole2=pole;
colors={'b','r','g','k','c','m','y'};
npole=size(pole,1);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
ms=50;%number of miliseconds
overlap=0;%miliseconds
poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
resultp=[];
for i=1:npole
    %%%%%%%%%%%%%%%raw data
    filename=strcat('./tr4s/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    %filename=strcat('./tr4s/processed/bandpass/upto5/result',pole(i,:),'.mat');
  
        d=strcat('./tr4s/processed/touchdetection/',pole(i,:),'/');
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
%             result=centredk(endp,vtouch.touches,result);
%             if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
% 
%             result=result(startp:endp,:);
            
            %%%%%%%%%%%%for touch periods
             result=centredk(endp,vtouch.touches,result);
             result=touchselection(vtouch.touches,result);
             
            im=sum(imag(result(:)));
            if im>0
                display('Complex numbers found')
            end
            %%%%%%%%%%%for non touching periods
            %result=notouchselection(vtouch.touches,result);
            
            %%%%%%%%%for debugging
           % plot(result(:,4))
            %result=centredk(endp,vtouch.touches,result);
%              if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
            %result=notouchselection(vtouch.touches,result);
            %result=touchselection(vtouch.touches,result);
                         %hold on
%                           result=result(startp:endp,:);
                         %plot(result(:,1),colors{i})
%                          hold off
%                         pause(5)
           
             resultp=[resultp; result];
            plot(resultp(:,1))
            pause(1)
            clear result
        end
            result=resultp;
 
            
            nsamples=floor(size(result,1)/ms);
            nresult(i)=nresult(i)+nsamples;
            result((nsamples*ms+1):end,:,:)=[];
    
            poleaz2=reshape(squeeze(result(:,1,:)),nsamples.*ms,1);
            poleaz2=reshape(poleaz2,ms,nsamples);
          
            
            poleel2=reshape(squeeze(result(:,2,:)),nsamples.*ms,1);
            poleel2=reshape(poleel2,ms,nsamples);
            
            polekcor2=reshape(squeeze(result(:,3,:)),nsamples.*ms,1);
            polekcor2=reshape(polekcor2,ms,nsamples);
            
            polekhor2=reshape(squeeze(result(:,4,:)),nsamples.*ms,1);
            polekhor2=reshape(polekhor2,ms,nsamples);
            
            poleaz=[poleaz,poleaz2];
            poleel=[poleel,poleel2];
            polekcor=[polekcor,polekcor2];
            polekhor=[polekhor,polekhor2];
     
       
        

   
   
    
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
azel=[zscore(tazimuth'),zscore(televation')];

    totalfiles=nresult(1);
    poleclass(1:totalfiles)={pole(1,:)};

    for i=2:npole
        poleclass(totalfiles+1:totalfiles+nresult(i))={pole(i,:)};
        totalfiles=totalfiles+nresult(i);
    end
%%plot
%%%%%%%%%Raw Data
% size(azel)
% plot(azel(1,:)','b')
% hold on
% plot(azel(nresult(1)+1,:)','r')
% plot(azel(1:nresult(1),:)','b')
% plot(azel(nresult(1)+1:nresult(1)+nresult(2),:)','r')
% hold off
%%%%Apply LDA angles
%%%%%%%linear
Mdlazel = fitcdiscr(azel,poleclass);
resuberror = resubLoss(Mdlazel)
Rel = confusionmat(Mdlazel.Y,resubPredict(Mdlazel))
%pval=Barletttest(azel,poleclass,totalfiles)

%%%%%quadratic
% Mdlazelq = fitcdiscr(azel,poleclass,'DiscrimType','quadratic');
% resuberror = resubLoss(Mdlazelq)
% Rel = confusionmat(Mdlazelq.Y,resubPredict(Mdlazelq))

%%%%curvatures
%%%%%%Kappa horizontal
Mdlkhor = fitcdiscr(tkhorizontal',poleclass);
resuberror = resubLoss(Mdlkhor)
Rel = confusionmat(Mdlkhor.Y,resubPredict(Mdlkhor))

%%%%kappa coronal
Mdlkcor = fitcdiscr(tkcoronal',poleclass);
resuberror = resubLoss(Mdlkcor)
Rel = confusionmat(Mdlkcor.Y,resubPredict(Mdlkcor))

%%%%%%curvatures
Mdlkcor = fitcdiscr([tkhorizontal',tkcoronal'],poleclass);
resuberror = resubLoss(Mdlkcor)
Rel = confusionmat(Mdlkcor.Y,resubPredict(Mdlkcor))

%%%%%%all variables
Mdlall = fitcdiscr([azel,zscore(tkhorizontal'),zscore(tkcoronal')],poleclass);
resuberror = resubLoss(Mdlall)
Rel = confusionmat(Mdlall.Y,resubPredict(Mdlall))


end
function [result]=centredk(endp,touches,result)
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
function pval=Barletttest(data,poleclass,totalfiles)
L = fitcdiscr(data,poleclass);
Q = fitcdiscr(data,poleclass,'DiscrimType','quadratic');
D = size(data,2); % Number of dimensions of X
Nclass = [totalfiles size(data,1)-totalfiles];
N = L.NumObservations;
K = numel(L.ClassNames);
SigmaQ = Q.Sigma;
SigmaL = L.Sigma;
logV = (N-K)*log(det(SigmaL));
for k=1:K
    logV = logV - (Nclass(k)-1)*log(det(SigmaQ(:,:,k)));
end
nu = (K-1)*D*(D+1)/2;
pval = 1 - chi2cdf(logV,nu);
end
