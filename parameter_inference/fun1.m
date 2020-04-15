function [ error1 ] = fun1(x, RFI,trec)
%objective function
error1=norm(RFI-(fun2(x,trec)),2);
%disp(error1);
end

