%%%%%%%%%%%%%%%%%inputs
I=[ones(5,10);zeros(5,10)];
phi=rand(10,10,100);
a_i=linspace(-1,1,100)';

%%%parameter
lambda=1;

%%%%%%%%%%%%model of image
Imodel=zeros(size(phi(:,:,1)));
for i=1:size(a_i)
    Imodel=Imodel+a_i(i)*phi(:,:,i);
end
Isquare=sum(sum((I-Imodel).^2));

%%%%%%%%%%%obj:minimize E with respect to a_i and phi
E=Isquare+lambda*sum(log(1+a_i.^2))

% %%%%%%%%%%%%%%%%parameters
% z_beta=1;
% beta=10;
% nbasis=size(a_i,1);
% z_n=1;
% varN=0.2;
% 
% %%%%%%%%%Probability distributions activity coefficients
% P_a_i=exp(-beta*log(1+a_i.^2));
% P_a=P_a_i.^nbasis;
% figure
% plot(a_i,P_a_i)

% 
% %%%%%%%%%%%%%%%%%probability distirbution of error on the image
% Isquare2=a_i;
% P_a_phi=exp(-(Isquare2).^2/(2*varN));
% figure
% plot(Isquare2,P_a_phi)
% 
% %%%%%%%Kullback_Leibler divergence 
% PI_phi=P_a.*P_a_phi;
% 
% KL=mean(log(PI_phi))
%%%%%%%%%%%%%%problem: find the basis phi that maximize KL