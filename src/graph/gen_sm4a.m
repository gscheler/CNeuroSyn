function [A,A_rdir] = gen_sm4(N,d),
% preferential attachment
% from [BatageljBrandes2004]
% generates undirected graph
%
%N=160;
%d=5;

m=zeros(2*N*d,1);

for v=0:N-1,
	for i=0:d-1,
		m(2*(v*d+i)+1) = v;
		r = floor(rand *2*(v*d+i));
		m(2*(v*d+i)+1+1) = m(r+1);
	end;
end;

A=zeros(N,N);
for i=0:N*d-1,
	A(m(2*i+1)+1,m(2*i+1+1)+1) = 1;
	end;

%
% make this graph random-directed
%
A_rdir=zeros(N,N);
for i=1:N,
	for j=1:N,
		if (A(i,j)),
			if (sum(A(i,:)) < sum(A(j,:))), % (rand > 0.5),
				A_rdir(j,i) = 1;
				A_rdir(i,j) = 0;
			else
				A_rdir(j,i) = 0;
				A_rdir(i,j) = 1;
			end;
		end;
	end;
end;




