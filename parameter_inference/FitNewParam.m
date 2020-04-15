function [p,fval] = FitNewParam( trec, RFI, startpoints, plotmode )

if size(RFI,1)==302
    RFI(1:52,:)=[];
    trec(1:52,:)=[];
    trec=trec-3.5020;
end
problem = createOptimProblem('fmincon',...
    'objective',@(x) fun1(x,RFI,trec),...
    'x0',startpoints,'Aineq',-[1 -1 0 -1 0;0 0 1 0 -1],'bineq',-[0 eps]',...
    'lb',[1 1e-4 1e-4 1e-4 1e-4], 'ub',[1 10 100 10 10],...
    'options',optimset('Algorithm','sqp'));
gs = GlobalSearch;
[p fval] = run(gs,problem);
if plotmode==1
figure;
hold on
plot(trec,RFI)
plot(trec,fun2(p,trec),'r')
set(gca,'ylim', [0 1.1])
end
end

