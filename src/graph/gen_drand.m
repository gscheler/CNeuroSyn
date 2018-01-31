%
% generate directed random graph with (N,K)
%
function A=gen_drand(N,K)

A=zeros(N,N);

p=0;
while p < K,
        i = floor(N*rand)+1;
        j = floor(N*rand)+1;
%printf("%d %d\n",i,j);
        if ((A(i,j) == 0) && (i ~= j)),
%printf(" set: %d %d\n",i,j);
                A(i,j) = 1;
                p = p+1;
                end;
        end;

