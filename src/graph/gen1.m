
function [G_LN,G_R] = gen1(N,N_target,thr,k)

%thr=70;

for i=k:k+10,
	[A,B] = gen_sm4a(N,i);
	A0 = chop(B,thr);
	if (size(A0,1) >= N_target),
		break;
	end;
end;


G_LN = A0(1:N_target,1:N_target);
[d,NN,K]=density(G_LN);
G_R = gen_drand(N_target,K);

[d,NN,K]=density(G_LN);
od=sum(G_LN);
sigma_star = exp(std(log(od(od>0))));
fprintf('LN Graph with N=%d K=%d\n',NN,K);
fprintf('LN density	= %f\n',d);
fprintf('LN sigma* 	= %f\n',sigma_star);

[d,NN,K]=density(G_R);
od=sum(G_R);
sigma_star = exp(std(log(od(od>0))));
fprintf('RND Graph with N=%d K=%d\n',NN,K);
fprintf('RND density	= %f\n',d);
fprintf('RND sigma* 	= %f\n',sigma_star);


	
