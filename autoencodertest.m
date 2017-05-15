function autoencodertest(X,T)
%[X,T] = wine_dataset;
limit=sum(T);
hiddenSize =2;
autoenc1 = trainAutoencoder(X',hiddenSize,...
    'L2WeightRegularization',0.1,...
    'SparsityRegularization',4,...
    'SparsityProportion',0.01,...
    'DecoderTransferFunction','logsig',...
    'EncoderTransferFunction','logsig',...
    'MaxEpochs',6000)
XReconstructed = predict(autoenc1,X');
mseError = mse(X'-XReconstructed)
figure
subplot(4,4,1)
plot(X(2,:))
hold on 
plot(XReconstructed(:,2))
legend('Input','Reconstruction','Location','bestoutside')
hold off
subplot(4,4,2)
plot(X(10,:))
hold on 
plot(XReconstructed(:,10))
hold off
subplot(4,4,3)
plot(X(100,:))
hold on 
plot(XReconstructed(:,100))
hold off
subplot(4,4,4)
plot(X(500,:))
hold on 
plot(XReconstructed(:,500))
hold off
subplot(4,4,5)
plot(X(1000,:))
hold on 
plot(XReconstructed(:,1000))
hold off
features1 = encode(autoenc1,X');
subplot(4,4,6)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,7)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,8)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,9)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,10)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,11)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,12)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,13)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,14)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,15)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
subplot(4,4,16)
i=floor(rand(1,1)*size(X,1));
plot(X(i,:))
hold on 
plot(XReconstructed(:,i))
hold off
features1 = encode(autoenc1,X');
 figure
 subplot(1,2,1)
 scatter(features1(1,1:limit),features1(2,1:limit),'b')
 hold on 
 %scatter3(features1(1,60:130),features1(2,60:130),features1(3,60:130))
 scatter(features1(1,limit+1:end),features1(2,limit+1:end),'r')
 legend('Air','Pole')
 title('Autoencoder')
% hold off

%figure
%scatter3(features1(4,1:limit),features1(5,1:limit),features1(6,1:limit))
%hold on 
%scatter3(features1(1,60:130),features1(2,60:130),features1(3,60:130))
%scatter3(features1(4,limit+1:end),features1(5,limit+1:end),features1(6,limit+1:end))
%hold off


hiddenSize = 2;
 autoenc2 = trainAutoencoder(features1,hiddenSize,...
     'L2WeightRegularization',0.1,...
    'SparsityRegularization',4,...
    'SparsityProportion',0.01,...
    'DecoderTransferFunction','purelin',...
    'EncoderTransferFunction','logsig',...
    'MaxEpochs',3000)

features2 = encode(autoenc2,features1);
% figure
% subplot(1,2,1)
% scatter(features2(1,1:limit),features2(2,1:limit))
% hold on 
% %scatter3(features2(1,60:limit),features2(2,60:limit),features2(3,60:130))
% scatter(features2(1,limit+1:end),features2(2,limit+1:end))
% hold off
%softnet = trainSoftmaxLayer(features2,T,'LossFunction','crossentropy');
%deepnet = stack(autoenc1,autoenc2,softnet);
% deepnet = train(deepnet,X',T');
% wine_type = deepnet(X);
% figure
% plotconfusion(T,wine_type);

%%%%%%%%%%%%%%%%PCA
subplot(1,2,2)
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(X);

scatter(SCORE(1:limit,1),SCORE(1:limit,2),[],'MarkerEdgeColor','b')
hold on 
%scatter3(SCORE(60:130,1),SCORE(60:130,2),SCORE(60:130,3),[],'MarkerEdgeColor','b')
scatter(SCORE(limit+1:end,1),SCORE(limit+1:end,2),[],'MarkerEdgeColor','r')
legend('Air','Pole')
title('PCA')
end