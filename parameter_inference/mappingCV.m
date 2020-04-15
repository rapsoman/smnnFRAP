function [net,PSinputs,PSppi,indTest, indTrain, perfNN, pind, thetatrain,ppitrain]=mappingCV(inputs,pv,pc,nclusters,ncurves,k)
%SCRIPT FOR THE NN TESTS
%prepare inputs for nnmap
%you have to load ThetaPvPc.mat 
%%
%first, normalize all
[inputs,PSinputs]=mapminmax(inputs);
pAc=cell2mat(pc);
pAc=pAc(2:5,:);
[pAc, PSppi]=mapminmax(pAc);
pcNorm=mat2cell(pAc,4,ncurves*ones(1,nclusters));
pvNorm=mapminmax('apply',pv(2:5,:),PSppi);
%%
%train with k-fold Cross Validation 
sizek=nclusters/k;
%permute all first
pind=randperm(nclusters);
inputsP=inputs(:,pind);
pcNormP=pcNorm(pind);
pvP=pvNorm(:,pind);
%%
m=0:sizek:nclusters;
for n=1:k

    ind = circshift(1:nclusters,[0 m(n)]);  
    %test set
    indTest(:,n)=ind(1:sizek);
    inputstest=inputsP(:,indTest(:,n));
        for i=1:sizek
            inputstest2{i}=repmat(inputstest(:,i),1,ncurves);
        end
    thetatest=cell2mat(inputstest2);
    ppitest=cell2mat(pcNormP(:,1:sizek));

    %train set
    indTrain(:,n)=ind(sizek+1:end);
    inputstrain=inputsP(:,indTrain(:,n));
        for i=1:nclusters-sizek
        	inputstrain2{i}=repmat(inputstrain(:,i),1,ncurves);
        end
    thetatrain{n}=cell2mat(inputstrain2);
    ppitrain{n}=cell2mat(pcNormP(:,indTrain(:,n)));

    %run network
    [ net{n}, tr] = nnmap(thetatrain{n},ppitrain{n},25,'tansig','tansig',PSinputs);
    perfNN{n}=tr.perf;


end


   


