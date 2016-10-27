function trackcorr(files)

[idx,result]=comparetr4(files,0);

nfiles=size(files,1);
nresult=size(result,3);
[acor,lag]=xcorr(result(:,:,1),'coeff');
%%%%%%%%%%%%%%%%%%%%%%%%correlation Azimuth
        figure
        plot(lag,acor(:,1:4))
        hold on 
        xlabel('lag')
        ylabel('correlation')
        title('Correlation Azimuth')
        l1=legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
       
%%%%%%%%%%%%%%%%%%%%%%%%correlation Elevation
      
        figure
        plot(lag,acor(:,5:8))
        hold on 
        xlabel('lag')
        ylabel('correlation')
        title('Correlation Elevation')
        legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
        
        
%%%%%%%%%%%%%%%%%%%%%%%%correlation Kappa coronal
        figure
        plot(lag,acor(:,9:12))
        hold on 
        xlabel('lag')
        ylabel('correlation')
        title('Correlation Kappa coronal')
        legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
       
%%%%%%%%%%%%%%%%%%%%%%%%%correlation Kappa horizontal        
        figure
        plot(lag,acor(:,13:16))
        hold on 
        xlabel('lag')
        ylabel('correlation')
        title('Correlation Kappa horizontal')
       legend('Azimuth','Elevation','Kappa coronal','Kappa horizontal');
        hold off
       
        
end