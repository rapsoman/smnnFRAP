

function out = randsrc(m,n,distrib) 

%out = randsrc(m,n,[alphabet; prob]) generates an m-by-n matrix, each of whose entries is independently chosen 
%from the entries in the row vector alphabet. Duplicate values in alphabet are ignored. 
%The row vector prob lists corresponding probabilities, so that the symbol alphabet(k) occurs with probability prob(k), 
%where k is any integer between one and the number of columns of alphabet. The elements of prob must add up to 1.
%clone of randsrc from communications toolbox

%[alphabet; prob]
alphabet=distrib(1,:);
prob=distrib(2,:);

NotAssigned=ones(m,n);
out=zeros(m,n);
F=cumsum(prob);
N=length(alphabet);

x=rand(m,n);
for j=1:N
    Candidates=(x<=F(j));
    NewEntries=Candidates.*NotAssigned;
    out=out+alphabet(j)*NewEntries;
    NotAssigned=NotAssigned-NewEntries;
end