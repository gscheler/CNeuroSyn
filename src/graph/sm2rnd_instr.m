%
% successively create random graph from A_sm toward A_r and
% measure sync as well as other graph properties
%
% N= number nodes
% std_target = desired sigma of target graph
%
% maximum of 5000 edge moves (hardcoded)
%
%


function [A,nn,ss] =sm2rnd(N,A_sm,A_r,std_target);
%std_target=7;
%std_target=13;

st_r=std(sum(A_r));
m_r=mean(sum(A_r));


A=A_sm;

std(sum(A))

ss=zeros(5000,1);
nn=zeros(5000,1);

touched = zeros(size(A));

kkk=0;
for i=1:5000,
	sa=sum(A);
	%thr = round(m_r + st_r*randn);
	thr = round(m_r);
	if (thr ==0),
		continue;
	end;

	
	id=find(sa>=thr);
	if (length(id)<10),
		continue;
	end;
	r=randperm(length(id));
	sel=id(r(1));

	id2 = find(sa < thr);
	sel2=id2;
	sel2=1:N;


	
	c=find(A(:,sel));
	I=randperm(length(sel2));
	J=randperm(length(sel2));
	kk=1;
	for k=sa(sel):-1:thr,
%		A(sel,c(kk)) = 0;
if (touched(c(kk),sel)),
	continue;
else
	touched(c(kk),sel);
end;

		A(c(kk),sel) = 0;
		if (kk>=length(sel2)),
			break;
		end;
		A(sel2(I(kk)),sel2(J(kk))) = 1;
		kk=kk+1;
		kkk=kkk+1;

		nn(kkk) = 1;
		od=sum(A);
		ss(kkk) = exp(std(log(od(od>0))));
	end;

%	if (kkk > 600),
%		break;
%	end;
	if (std(sa) < std_target),
		kkk
		break;
	end;
	if (1 == 0 & mod(i,10)==1),
		fprintf('%g  %g\n',mean(sa),std(sa));
		hist(sa,30);
		drawnow;
		input('xxx');
	end;
end;
%hist(sa,30);
		
