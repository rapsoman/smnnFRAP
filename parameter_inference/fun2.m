function [ y ] = fun2(x,trec)
%returnes the curve for a given set of parameters
y=x(1)-x(2)*exp(-x(3)*trec)-x(4)*exp(-x(5)*trec);