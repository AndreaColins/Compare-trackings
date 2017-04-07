function timelag=trackcorr(files)

[~,result]=comparetr4(files,0);
nfiles=size(files,1);
acor=zeros(size(xcorr(result(:,:,1)),1),16,nfiles);


for i=1:nfiles
[acor(:,:,i),lag]=xcorr(result(:,:,i),'coeff');
[~,idlag(i,:)]=max(abs(acor(:,:,i)));
end

timelags=idlag-ones(size(idlag))*3488;


macor=zeros(size(acor(:,:,1)));
size(acor,2)

for i=1:size(acor,2)
macor(:,i)=mean(squeeze(acor(:,i,:)),2);
stdcor(:,i)=std(squeeze(acor(:,i,:)),0,2);
end

%timelag de el promedio de las correlaciones
[~,I] = max(abs(macor));
timelag=lag(I);
%seria bueno calcular la desviacion estandar de el time lag
%calcular el timelag de cada correlacion y para cada uno de esos vectores
%calcular la desviacion estandar
%%%%%%%%%%%%%%%%%%%%%%%%correlation Azimuth


        figure
        %patch([lag fliplr(lag)],[y+err fliplr(y-err)],[0.7 0.7 0.7]);
        hold on 
         
        % patch([lag fliplr(lag)],[macor(:,1)'+stdcor(:,1)' fliplr(macor(:,1)-stdcor(:,1))'],[ 0 0 0.1],'LineStyle' ,'none');
         %patch([lag fliplr(lag)],[macor(:,2)'+stdcor(:,2)' fliplr(macor(:,2)-stdcor(:,2))'],[0 0.1 0 ],'LineStyle' ,'none');
        plot(lag,acor(:,1:4))
        %errorbar(lag,macor(:,1:4),stdcor(:,1:4))
         xlabel('lag')
        ylabel('correlation')
        title('Correlation Azimuth')
        l1=legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
       
%%%%%%%%%%%%%%%%%%%%%%%%correlation Elevation
      
        figure
        plot(lag,macor(:,5:8))
        hold on 
        xlabel('lag')
        ylabel('correlation')
        title('Correlation Elevation')
        legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
        
        
%%%%%%%%%%%%%%%%%%%%%%%%correlation Kappa coronal
        figure
        plot(lag,macor(:,9:12))
        hold on 
        xlabel('lag')
        ylabel('correlation')
        title('Correlation Kappa coronal')
        legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
       
%%%%%%%%%%%%%%%%%%%%%%%%%correlation Kappa horizontal        
        figure
        plot(lag,macor(:,13:16))
        hold on 
        xlabel('lag')
        ylabel('correlation')
        title('Correlation Kappa horizontal')
       legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
       
        
end