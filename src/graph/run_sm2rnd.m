
%
% successively create random graph from LG1 and
% measure sync as well as other graph properties
% N_nodes = 210 (hardcoded)
%
% plots dispersion for the graphs
%
% uses sm2rnd_instr(...)
%

% MUST PRELOAD SOME A_* stuff
% thr=70;
% [A_LN, A_R] = gen1(350,210,thr, 1);
%
% load RG_hist.mat
% load LN_hist.mat

load sm2rnd_exp2.mat

	% number of intermediate steps
N = 10;

	% sigma^2 for the target graph
s_target = 3;
s_target = 4;

ss=zeros(5000,N);
nnn=zeros(N,1);

for i=1:N,
	if (i>1),
	[A,nn1,ss1]= sm2rnd_instr(210,A_LN,A_R,s_target);

	l = max(find(nn1));
	else
	l=1200;
	ss1(:,1)=2.89+zeros(5000,1);
	end;

	nnn(i)=l;
%	ss(1:l,i)=ss1(l:-1:1);
	ss(1:l,i)=ss1(1:l);
end;

figure;
for i=2:N,
%	plot(ss(1:nnn(i),i));
	plot(1200:-1:1,ss(1:1200,i),'k--');
	hold on
end;

	
ms=mean(ss(1:1200,2:N)');
%figure;
plot(1200:-1:1,ms,'r','Linewidth',3);

xlabel('percentage of edges changed')
ylabel('dispersion \sigma^*')
set(gca,'Fontsize',14);


set(gca,'YTick',[1.44,2,2.5,2.89,3]);
set(gca,'YTickLabel',{'RG',2,2.5,'LG1',3});

K=sum(sum(A_LN));
set(gca,'XTick',[0,0.2*K,0.4*K,0.6*K]);
set(gca,'XTickLabel',{0,20,40,60});

axis([0,1200,1.3,3])
