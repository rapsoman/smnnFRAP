function [pGSnew,ac,pAc]=artificialCurves2(trec,RFI,startpoints,ncurves,nclusters)
%computes the artificial curves, their respective parameter values and the
%confidence intervals with the bootstrapping way - initial conditions
%corrected

%INPUTS:
%alpha: confidence level (usually 0.95)
%RFI, trec: matrices of simulated curves - time vector sould be set so as
%t0=t(end of frap)
%ncurves: number of artificial curves
%nclusters: number of clusters: 100 for the scattered and 125 for the
%gridded data

%OUTPUTS:
%pGS: parameter vectors for the initial fit
%ac: artificial curves
%pAc: matrix of parameter vectors of each ac
%ci2: confidence intervals


%confidence level
% beta=(1-alpha)/2;
% minindex=round(beta*ncurves);
% maxindex=round((1-beta)*ncurves);
% 
% %preallocation
% pGSnew=zeros(5,nclusters);
% ac=cell(1,nclusters);
% pAc=cell(1,nclusters);
% ci2=cell(1,nclusters);
% pAc_sort=cell(1,nclusters);

for j=1:nclusters
j
%find the parameters for the curve of interest
[pGSnew(:,j)] = FitNewParam( trec(:,j), RFI(:,j), startpoints(:,j), 1 );
% %create artificial curves
[ac{1,j}]=boottest(trec(:,j),RFI(:,j),pGSnew(:,j),ncurves);
%find the parameters for the artificial curves
for i=1:ncurves
    
[pAc{1,j}(:,i)] = FitNewParam(trec(:,j), ac{j}(:,i),pGSnew(:,j), 0 );
end
% pAc_sort{1,j}=sort(pAc{1,j},2);
% ci2{1,j}(:,2)=pAc_sort{1,j}(:,maxindex);
% ci2{1,j}(:,1)=pAc_sort{1,j}(:,minindex);
end




