%
% run simulations for different graphs and different densities
% generates "res" data structure
%
N=210;

FN_GRAPHS='Y_210_15_sm4a';

load(sprintf('%s.dat',FN_GRAPHS),'-mat');

verbose=0;


W0=[0.06,0.05,0.04];

kS=10:2:30;
sS=1:5:45;

% for selection and table:
kS=15
kS=14
%sS=1:50:45;
%sS=[10,30];
sS=10;

dodo=1;


%-------------------------------


run=1;

for run_s=sS,
	i_k =1;
for run_k=kS,
	thr = 5* run_k;
%	thr = 4* run_k;
%	thr = 2* run_k;
%	thr = 10* run_k;
%thr=999;
%	thr = 4* run_k -5;
%	[A_LN, A_R] = gen1(220,210,thr, run_k);
	[A_LN, A_R] = gen1(350,210,thr, run_k);
%	[A_LN, A_R] = gen1(220,210,thr, run_k);

	A = A_LN;

%	A = A_R;

%	A = sm2rnd(N,A_LN,A_R,run_s);

	FN_GRAPHS=sprintf('foo_G_foo3_MIX_%d',run);
	fn_graphs=sprintf('%s.dat',FN_GRAPHS);
	save(fn_graphs, 'A', 'E', 'E_inh', 'E_bg');

	W=0.08;
	W=0.04;
	W=0.06;

	[kden,N,K] = density(A);

%	W= W0(run);
%	W= 0.04 + (run-1)*0.02/7;
%	W= 0.06 - (run-1)*0.005/7;
	W= 0.06 - (i_k-1)*0.02/length(kS);
	i_k = i_k+1;

	res.W(run) = W;

fprintf('WEIGHT w=%f\n',W);
	
if (dodo ==1),
% $Id:$
% gr7_001
%
% input for 1st 1000 ms; then no input afterwards
%

do_sim =         1;
do_gen_inputs =  0;
do_save_inputs = 0;
do_gen_graphs =  0;
do_save_graphs = 0;

sel_input=1;


% random graph
%FN='gr7_001_5000_rnd_002';
%FN_GAUSS='gr7_gauss_001_rnd';
%FN_GRAPHS='gr7_graph_002_rnd';

%FN='../0709/gr7_001_5000_sw';
%ORIG FN='rnd_g16_001';
%ORIG FN_GAUSS='gaussian_rnd_g16_001';
%FN_GRAPHS='../0709/gr7_graph_001_sw';
%FN_GRAPHS='../0709/gr7_graph_001_rnd';

%FN_GRAPHS='gr7_graph_002_rnd';
%FN_GRAPHS='gr7_001_5000_rnd';
%
%FN_GRAPHS='gr7_graph_002_rnd';

%-------------------------------------------
% take that graph and run several mods on that
%
%FN_GRAPHS='gr7_graph_001_sw';
%FN_GRAPHS='X_2016_swg1_12';
%FN_GRAPHS='X_2016_swg2_10';
%FN_GRAPHS='X_2016_swg3_7';
%FN_GRAPHS='X_2016_swg4_5';
%FN_GRAPHS='gr7_001_5000_rnd';

%FN_GRAPHS='Y_210_15_sm4a';
%FN_GRAPHS='Y_210_15_drand';

% reuse OLD input
%
FN_INP='gr_inp_104d';
FN_INP='gr_inp_104u';

FN_PLOT='foofoo';
FN_GAUSS='foo';

FN='foobar';


sim.exp = 'Description';

N_nn	= 210;	% number of neurons
N_nn_inh= 210;	% number of neurons
T_upd   = 5000;  % length of each update cycle [ms]

N_upd	= 1;	% number of update cycles

	%
	% definition of graph parameters
	%
%sim.reg.k = 11;
sim.reg.k = 10;
sim.reg.k = 9;
sim.reg.k = 8;
%BUILD_GRAPH='build_graph_r(sim,mksmall_pa_u(sim.N_nn, sim.reg.k))';

%BUILD_GRAPH='build_graph_r(sim,mkrand(sim.N_nn, 800))';
%BUILD_GRAPH='build_graph_r(sim,mkrand(sim.N_nn, 1600))';
%BUILD_GRAPH='build_graph_r(sim,mkrand(sim.N_nn, 1800))';

	% definition of neuron mu parameters
nn_mu_params=zeros(N_nn,4);

nn_mu_params(:,1) = 0.02;     % a
nn_mu_params(:,2) = 0.2;      % b

r=0.0;
%r=0.08*rand(N_nn,1);
nn_mu_params(:,3) = - 65 + 15*r.^2;  %c
nn_mu_params(:,4) = + 8  - 6*r.^2 ;  % d

%-----------------------------------------

nn_mu_params(1:30,1) = 0.022;
nn_mu_params(1:30,2) = 0.3;
nn_mu_params(1:30,4) = 9.5;

nn_mu_params(31:60,1) = 0.022;
nn_mu_params(31:60,2) = 0.3;
nn_mu_params(31:60,4) = 14;

nn_mu_params(61:90,1) = 0.015;
nn_mu_params(61:90,2) = 0.15;
nn_mu_params(61:90,4) = 14;

nn_mu_params(91:120,1) = 0.015;
nn_mu_params(91:120,2) = 0.2;
nn_mu_params(91:120,4) = 12;

nn_mu_params(121:150,1) = 0.02;
nn_mu_params(121:150,2) = 0.2;
nn_mu_params(121:150,4) = 9;

nn_mu_params(151:180,1) = 0.025;
nn_mu_params(151:180,2) = 0.2;
nn_mu_params(151:180,4) = 6;

nn_mu_params(181:210,1) = 0.02;
nn_mu_params(181:210,2) = 0.2;
nn_mu_params(181:210,4) = 8;



%-----------------------------------------

w_ext = zeros(N_nn,1) +0.08;
w_ext_inh = zeros(N_nn,1) - 0.002;
w_bg = zeros(N_nn,1) + 0.002;
w_net = zeros(N_nn,1) + W;
w_ext_inh = zeros(N_nn,1) - 0.002;
%w_bg = zeros(N_nn,1) + 0.02;


sim.delay = 7;


	%
	% definition of external inputs
	% for excitatory 
	%
	% medium correlated Poisson
	%

if (1==0),
input_params.Mp = 50;   % Mp
input_params.Mn = 10;   % Mn
input_params.lambdan = 400;     % lambda_n
input_params.lambdap = 350;     % lambda_p
input_params.corrp = 0.8;       % rel. correlation for Mp
input_params.corrn = 0.8;       % rel. correlation for Mn
input_params.g0 = -2.0; % g_0
input_params.eta        = 0.3;          % sigma^2 of randn noise
input_params.start      = 200;  % start offset
input_params.dc = -0.39;
input_params.sin_freq =  25; %Hz
input_params.sin_ampl = 0.0;
input_params.sin_width=20;  %ms
end;

%============= driven ===================

if (1==1),
input_params.Mp = 50;   % Mp
input_params.Mn = 10;   % Mn
input_params.lambdan = 400;     % lambda_n
input_params.lambdap = 350;     % lambda_p
%input_params.corrp = 0.6;      % rel. correlation for Mp
%input_params.corrn = 0.6;      % rel. correlation for Mn
input_params.corrp = 0.8;       % rel. correlation for Mp
input_params.corrn = 0.8;       % rel. correlation for Mn
%input_params.g0 = -4.5;        % g_0
%input_params.g0 = 0;   % g_0
input_params.g0 = -2.5; % g_0
input_params.g0 = -2.0; % g_0
input_params.g0 = -2.5; % g_0
input_params.eta        = 0.3;          % sigma^2 of randn noise
input_params.eta        = 0.1;          % sigma^2 of randn noise
input_params.start      = 200;  % start offset
input_params.dc = -0.30;
input_params.dc = -0.39;
input_params.dc = -0.35;
input_params.dc = -0.30;
input_params.sin_freq =  25; %Hz
input_params.sin_ampl = -2.2;
input_params.sin_ampl = -2.0;
input_params.sin_ampl = -2.2;
input_params.sin_ampl = -2.4;
input_params.sin_ampl = -2.3;
input_params.sin_ampl = -2.2;
input_params.sin_width=20;  %ms
input_params.sin_width=10;  %ms
%input_params.sin_width=8;  %ms
end;

offset = 1;

%============= /driven ===================
	%
	% definition of external inputs
	% for inhib.
	%
	% medium correlated Poisson
	%

input_params_inh.Mp = 50;   % Mp
input_params_inh.Mn = 10;   % Mn
input_params_inh.lambdan = 400;     % lambda_n
input_params_inh.lambdap = 350;     % lambda_p
input_params_inh.corrp = 0.8;       % rel. correlation for Mp
input_params_inh.corrn = 0.8;       % rel. correlation for Mn
input_params_inh.g0 = -2.0; % g_0
input_params_inh.eta        = 0.3;          % sigma^2 of randn noise
input_params_inh.start      = 200;  % start offset
input_params_inh.dc = -0.39;
input_params_inh.sin_freq =  25; %Hz
input_params_inh.sin_ampl = 0.0;
input_params_inh.sin_width=20;  %ms

	%
	% definition of background
	% for inhib.
	%
	% medium correlated Poisson
	%

input_params_bg.Mp = 50;   % Mp
input_params_bg.Mn = 10;   % Mn
input_params_bg.lambdan = 400;     % lambda_n
input_params_bg.lambdap = 350;     % lambda_p
input_params_bg.corrp = 0.8;       % rel. correlation for Mp
input_params_bg.corrn = 0.8;       % rel. correlation for Mn
input_params_bg.g0 = -2.0; % g_0
input_params_bg.eta        = 0.3;          % sigma^2 of randn noise
input_params_bg.start      = 200;  % start offset
input_params_bg.dc = -0.39;
input_params_bg.sin_freq =  25; %Hz
input_params_bg.sin_ampl = 0.0;
input_params_bg.sin_width=20;  %ms






offset = 1;


%-------------------------------------------------------------


run_izh_graph7_5000_g16;

sim.isi_min = 5;
sim.isi_max=300;
%
if (verbose),
[mn,st] = plot_nfig11(FN_PLOT, sim, all_nn_inputs', 1, 1, 1);
axis([1101,2100,0,210]);
%%%%%%%plot_gaussian_ax(FN_GAUSS,mn,st,0,300);
end;

sim.isi_max=300;
sim.isi_min=3;

[m,s] = calc_isi_plot(sim,1000,0);

mm(1) = mean(m(1:30));
ss(1) = std(m(1:30));

mm(2) = mean(m(31:60));
ss(2) = std(m(31:60));

mm(3) = mean(m(61:90));
ss(3) = std(m(61:90));

mm(4) = mean(m(91:120));
ss(4) = std(m(91:120));

mm(5) = mean(m(121:150));
ss(5) = std(m(121:150));

mm(6) = mean(m(151:180));
ss(6) = std(m(151:180));

mm(7) = mean(m(181:200));
ss(7) = std(m(181:200));

if (verbose), plot_gaussian_ax_2018(FN_GAUSS, mm,ss,10,270); end;

fprintf('%d (generic), %d (type1), %d (type2), %d (type3), %d (type4), %d (type5), %d (type6)\n', ...
	round(mm(7)), round(mm(6)), round(mm(5)), round(mm(4)), ...
	round(mm(3)), round(mm(2)), round(mm(1)));



c_max = 5;
W = 5;

sim_rnd = sim;

[m_rnd, w_act_rnd] = calc_spikecorr2(sim_rnd, 1000, c_max,W);

sp_rnd = sum(w_act_rnd,2);
mm_rnd=m_rnd;
for j=1:sim.N_nn,
	if (sp_rnd(j) > 0),
		mm_rnd(j,:) = mm_rnd(j,:) ./sp_rnd(j);
	end;
end;
gs_rnd =sum(sum(mm_rnd))/((sim.N_nn^2-sim.N_nn));

fprintf('G   %.2f\n', gs_rnd);

sync_G=gs_rnd;

%----------------------------------
ds=52;
gr=4;
for i=1:gr,
	sel=(i-1)*ds+1:i*ds;

	sp_rnd = sum(w_act_rnd(sel,:),2);
	mm_rnd=m_rnd(sel,sel);

	for j=1:ds,
		if (sp_rnd(j) > 0),
			mm_rnd(j,:) = mm_rnd(j,:) ./sp_rnd(j);
		end;
	end;
	gs_rnd=sum(sum(mm_rnd))/((ds^2-ds));

	fprintf('Gr%d %.2f\n',i, gs_rnd);
end;

else

sync_G= -1;
end;


[kden,N,K] = density(A);

fprintf('Graph with N=%d K=%d\n',N,K);
fprintf('density	= %f\n',kden);

od=sum(A);
sigma_star = exp(std(log(od(od>0))));
fprintf('sigma* 	= %f\n',sigma_star);

% find indegree, outdegree, and degree for each vertex
[id,od,deg,J,J_od,J_id,J_bl] = degrees(A);
fprintf('mean indegree	= %f [%d %d]\n',mean(id), min(id), max(id));
fprintf('mean outdegree	= %f [%d %d]\n',mean(od), min(od), max(od));

% +++3.9+++ Cluster Index
disp(['clustind...']);
[gamma,gammaG] = clustind(A);
fprintf('cluster index	= %f\n',gammaG);

[R,D] = breadthdist(A);

[lambda,ecc,radius,diameter] = charpath(D);
fprintf('char. path length= %f\n',lambda);
fprintf('radius 	= %f\n',radius);
fprintf('diameter 	= %f\n',diameter);

[strconn,mulcomp,ncomps,ind_comps,siz_comps] = components(A,R);
fprintf('# components	= %f\n',ncomps);

%-----------------------
res.N(run) = N;
res.K(run) = K;
res.in_mean(run) = mean(id);
res.in_min(run) = min(id);
res.in_max(run) = max(id);
res.out_mean(run) = mean(od);
res.out_min(run) = min(od);
res.out_max(run) = max(od);
res.kden(run) = kden;
res.sigma_star(run) = sigma_star;
res.ncomps(run) = ncomps;
res.radius(run) = radius;
res.lambda(run) = lambda;
res.clust_index(run) = gammaG;
res.sync_G(run) = sync_G;
%res.W(run) = W;
	

	run=run+1;

end;
end;
