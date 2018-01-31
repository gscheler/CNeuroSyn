%
% preferential attachment
% from [BatageljBrandes2004]
% generates undirected graph
%
%N=210;
d=5;

N=1000;
d=12;

m=zeros(2*N*d,0);

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
for i=1:N,
	for j=1:N,
		if (A(i,j)),
			if (rand > 0.5),
				A(j,i) = 1;
				A(i,j) = 0;
			end;
		end;
	end;
end;




