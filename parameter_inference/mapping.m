function [net,pind,PSinputs,PSppi]=mapping(inputs,pv,pc,nclusters,ncurves)
%SCRIPT FOR THE NN TESTS
%%
%prepare inputs for nnmap
%you have to load ThetaPvPc.mat 

%%
%first, normalize all

[inputs2,PSinputs]=mapminmax(inputs2);
pAc=cell2mat(pc2);
pAc=pAc(2:5,:);
[pAc, PSppi]=mapminmax(pAc);
pcNorm=mat2cell(pAc,4,ncurves*ones(1,nclusters));


%%
%train with an external independent test set
pind=randperm(nclusters);
%train set
inputstrain=inputs(:,pind(1:nclusters*0.8));
for i=1:nclusters*0.8
inputstrain2{i}=repmat(inputstrain(:,i),1,ncurves);
end
thetatrain=cell2mat(inputstrain2);
ppitrain=cell2mat(pcNorm(pind(1:nclusters*0.8)));

%test set
inputstest=inputs(:,pind((nclusters*0.8)+1:nclusters));
for i=1:nclusters*0.2 
inputstest2{i}=repmat(inputstest(:,i),1,ncurves);
end
thetatest=cell2mat(inputstest2);
ppitest=cell2mat(pcNorm(pind((nclusters*0.8)+1:nclusters)));
%run network
[ net, tr, testP, testD ] = nnmap(thetatrain,ppitrain,25,'tansig','tansig',PSinputs);
%%
%simulate the network and get test outputs
predicted = net(ppitest); 
%de-normalize
predicted=mapminmax('reverse',predicted,PSinputs);
thetatest=mapminmax('reverse',thetatest,PSinputs);
figure;
%plot
tit={'D_{eff}','F_{imm}','T_{res}'};
for j=1:3
subplot(3,1,j)
hold on
plot(thetatest(j,:),'b')
plot(predicted(j,:),'r')
xlabel(tit(j))
end 
%%
% inputs=mapminmax('reverse',inputs,PSinputs);
% PredAll = net(cell2mat(pcNorm)); %simulate for all clusters
% PredAll=mapminmax('reverse',PredAll,PSinputs);
% PredAll=mat2cell(PredAll,3,ncurves*ones(1,nclusters));

%%
% %confidence intervals of predictions
% alpha=0.95;
% beta=(1-alpha)/2;
% minindex=round(beta*ncurves);
% % maxindex=round((1-beta)*ncurves);
% % for j=1:ncurves
% % pred_sort{1,j}=sort(PredAll{j},2);
% % pred_ci{1,j}(:,2)=pred_sort{1,j}(:,maxindex);
% % pred_ci{1,j}(:,1)=pred_sort{1,j}(:,minindex);
% % end
% %%
% %plots of prediction clusters
% pvNorm=mapminmax('apply',pv(2:5,:),PSppi);
% Pred = net(pvNorm); %simulate the network and get predicted outputs
% PredInitial=mapminmax('reverse',Pred,PSinputs);
% figure;
% tit={'D_{eff}','F_{imm}','T_{res}'};
% for j=1:3
% subplot(3,1,j)
% hold on
% bar([PredInitial(j,pind); inputs(j,pind)]')
% ylabel(tit(j))
% %set(gca,'xlim',[0 112])
% colormap gray
% end 
% 
% for i=1:nclusters,
%     for j=1:3,
%         reler(j,i)=(PredInitial(j,i)-inputs(j,i))./inputs(j,i);
%     end,
% end
% 
% figure;
% for i=1:3,subplot(3,1,i),hold on,bar(reler(i,pind)),set(gca,'xlim',[0 112]),end
% % set(gca,'XTick',[0:1:100])
% % set(gca,'XTickLabel',num2cell(pind))
% 
% 
% 
% 



