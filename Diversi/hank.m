% Generation of Hankel matrices
% function H = hank(X,col,row)
% H = [  X(col)      X(col-1)    ...  X(1)
%       X(col+1)      X(col)     ...  X(2)
%            .           .             .
%            .           .             .
%            .           .             .
%      X(col+row-1) X(col+row-2) ... X(row)]    
% X --> data vector
% col --> number of columns
% row --> number of rows

function H = hank(X,col,row)

N=length(X);

if row > N-col+1
    error('Too many rows. row must be <= length(X)-col+1');
end

H=[X(1:row)];

for i=1:col-1

H=[H X(i+1:row+i)];

end

H=fliplr(H);