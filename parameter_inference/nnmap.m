function [ net, tr, testP, testD ] = nnmap(theta,ppi,n,TransferFcn1,TransferFcn2,PSinputs)
%mapping using neural networks
%INPUTS: Theta parameters, P parameters(ppi), number of neurons (n)
%TransferFcn1-2: basis function for the first-second layer (hidden), possible
%inputs: 'logsig','radbas','poslin'
%PSinputs: struct for the normalization from mapminmax
out = theta;
in = ppi;
% out=theta; %network output
% in=ppi; %network input
%initialize
net = fitnet(n); %n neurons
net.layers{1}.transferFcn = TransferFcn1;
net.layers{2}.transferFcn = TransferFcn2;
% net.divideParam.valRatio=0;
net.trainParam.goal = 0.005;
net.trainParam.max_fail=100;
net.trainParam.mu_max=1e+15;
net.trainParam.min_grad=0.0001;
net.performParam.regularization=0;
[net,tr] = train(net,in,out);

testI = in(:,tr.testInd); %get test inputs
testD = out(:,tr.testInd); %get test outputs (desired) 
testP = net(testI); %simulate the network and get predicted outputs
testP=mapminmax('reverse',testP,PSinputs);
testD=mapminmax('reverse',testD,PSinputs);

tit={'Deff','Fimm','Timm'};
for j=1:size(theta,1)
subplot(3,1,j)
hold on
plot(testD(j,:),'b')
plot(testP(j,:),'r')
title(tit(j))
end 




