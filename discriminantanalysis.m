function MdlLinear=discriminantanalysis
%%%%try example with 1 result i.e. 1 file, it should show clusters at least
%%%%for touching and no touching periods
 touchfile=strcat('.\tr4s\processed\touchdetection\smooth pole\','TT3_20160812_133608','_touch.mat');
 filename=strcat('./tr4s/processed/FrameCorrected/resultsmooth pole.mat');
 %filename=strcat('./tr4s/processed/bandpass/upto5/resultsmooth pole.mat');
vtouch=load(touchfile,'-mat');
v=load(filename,'-mat');

idx=[1:3488];
idx2=circshift(idx,size(idx,1)-vtouch.start_frame-1,2);
I=vtouch.touches(idx2);
result=v.result(:,:,1);
azimuth=result(:,1);
elevation=result(:,2);

figure(1)
plot(azimuth,'b')
hold on 
plot(elevation,'r')
hold off

for i=1:3488
    if I(i)==1
    Iname(i,1)={'Touch'};
    else
    Iname(i,1)={'no Touch'};   
    end
end

%%%%graph original data
h1 = gscatter(azimuth,elevation,I,'krb','ov^',[],'off');
h1(1).LineWidth = 2;
h1(2).LineWidth = 2;
legend('No touch','Touch','Location','best')
hold on

%%%%%apply discriminant analysis
X = [azimuth, elevation];
I2=find(vtouch.touches(idx2));
I3=find(~vtouch.touches(idx2));
class1=X(I2,:);
class2=X(I3,:);
MdlLinear = fitcdiscr(X,Iname)
%%%%%%%%%%%%%%plot linear boundary between first and second class
%MdlLinear.ClassNames
%MdlLinear.Coeffs
K = MdlLinear.Coeffs(1,2).Const
L = MdlLinear.Coeffs(1,2).Linear
f = @(x1,x2) K+ L(1)*x1 + L(2)*x2;
h2 = ezplot(f,[0 110 0 0.001]);
h2.Color = 'g';
h2.LineWidth = 2;
Mdlazel=fitcsvm(X,I,'KernelFunction','rbf',...
    'BoxConstraint',Inf,'ClassNames',[0,1]);
x=76.88885:0.00001:76.88889;
y=L(2)*x/L(1)+7*1e6;
%plot(x,y)

V=[-L(2),L(1)]'
class1p=V(:,1)'*class1';
class2p=V(:,1)'*class2'

figure
plot(class1p,'r')
hold on
plot(class2p,'k')
hold off

figure
hi1=histogram(class1p)
hold on
hi2=histogram(class2p)
legend('Touch','No Touch','Location','best')
hold off

% f2 = @(x1,x2)-L(2)*x1 + L(1)*x2;
% h2 = ezplot(f2,[75 110 0 100]);
% h2.Color = 'r';
% h2.LineWidth = 2;
% f2(100,0.001)
% %%%%%%% plot linear boundary between first and second class
% f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% h3 = ezplot(f,[.9 7.1 0 2.5]);
% h3.Color = 'k';
% h3.LineWidth = 2;
% axis([.9 7.1 0 2.5])
% xlabel('Petal Length')
% ylabel('Petal Width')
% title('{\bf Linear Classification with Fisher Training Data}')
end


%%%%discriminant analysis example
% load fisheriris
% PL = meas(:,3);
% PW = meas(:,4);
% 
% figure
% plot(PL,'b')
% hold on 
% plot(PW,'r')
% hold off
% %%%%graph original data
% h1 = gscatter(PL,PW,species,'krb','ov^',[],'off');
% h1(1).LineWidth = 2;
% h1(2).LineWidth = 2;
% h1(3).LineWidth = 2;
% legend('Setosa','Versicolor','Virginica','Location','best')
% hold on
% 
% %%%%%apply discriminant analysis
% X = [PL,PW];
% MdlLinear = fitcdiscr(X,species);
% %%%%%%%%%%%%%%plot linear boundary between first and second class
% MdlLinear.ClassNames([2 3])
% K = MdlLinear.Coeffs(2,3).Const;
% L = MdlLinear.Coeffs(2,3).Linear;
% f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% h2 = ezplot(f,[.9 7.1 0 2.5]);
% h2.Color = 'r';
% h2.LineWidth = 2;
% 
% %%%%%%% plot linear boundary between first and second class
% f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% h3 = ezplot(f,[.9 7.1 0 2.5]);
% h3.Color = 'k';
% h3.LineWidth = 2;
% axis([.9 7.1 0 2.5])
% xlabel('Petal Length')
% ylabel('Petal Width')
% title('{\bf Linear Classification with Fisher Training Data}')
% end