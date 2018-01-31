%-------------------------------------------------------------
% run script for running IZH in a graph structure
% for gr7... scripts
% **** stop input after 1000
%-------------------------------------------------------------
% NOTE: background to all
%

fprintf('JSC: FN=%s\n',FN)

fn=sprintf('%s.dat', FN);
fprintf('JSC: fn=%s\n',fn)
fn_inp=sprintf('%s.dat', FN_INP);
fn_graphs=sprintf('%s.dat', FN_GRAPHS);

if (do_sim == 0),
	jsc
	load(fn, '-mat');
else 

	ts = 1;

	sim.N_nn = N_nn;
	sim.N_nn_inh = N_nn_inh;
	sim.T_upd = T_upd;
	sim.ts = ts;

	sim.activity_thr =0;
	sim.activity_win = 3;
	sim.isi_max=100;
	sim.isi_min=3;
	sim.isi_min=11;

	N_upd = 1;
%---------------------------------------------------------------
	%
	% initialize instrumentation
	%
	sim.instrument.saved_act = zeros(N_upd,N_nn);
	
	%
	% generate the inputs
	%
	if (do_gen_inputs == 0),
%OK JSC
		load(fn_inp, '-mat');
	else 
   		%----------------------------------------
                % create noisy poisson_distributed input
		% 
                %
                inp= ...
                  input_params.g0 * ...
                   inp_poisson(sim.T_upd+input_params.start, ...
                        input_params.Mp, ...
                        input_params.Mn, ...
                        input_params.lambdap, ...
                        input_params.lambdan, ...
                        1);

                inp = inp + ...
                        input_sin(sim.T_upd+input_params.start,input_params.sin_freq, input_params.sin_width, input_params.sin_ampl) + ...
                        input_params.dc + ...
                        input_params.eta*randn(sim.T_upd+input_params.start,1);


                inp_inh= ...
                  input_params_inh.g0 * ...
                   inp_poisson(sim.T_upd+input_params.start, ...
                        input_params_inh.Mp, ...
                        input_params_inh.Mn, ...
                        input_params_inh.lambdap, ...
                        input_params_inh.lambdan, ...
                        1);

                inp_inh = inp_inh + ...
                        input_sin(sim.T_upd+input_params.start,input_params_inh.sin_freq, input_params_inh.sin_width, input_params_inh.sin_ampl) + ...
                        input_params_inh.dc + ...
                        input_params_inh.eta*randn(sim.T_upd+input_params.start,1);



                inp_bg= ...
                  input_params_inh.g0 * ...
                   inp_poisson(sim.T_upd+input_params.start, ...
                        input_params_inh.Mp, ...
                        input_params_inh.Mn, ...
                        input_params_inh.lambdap, ...
                        input_params_inh.lambdan, ...
                        1);

                inp_bg = inp_bg + ...
                        input_sin(sim.T_upd+input_params.start,input_params_inh.sin_freq, input_params_inh.sin_width, input_params_inh.sin_ampl) + ...
                        input_params_inh.dc + ...
                        input_params_inh.eta*randn(sim.T_upd+input_params.start,1);


		all_nn_inputs = inp(input_params.start+1:end);
		all_nn_inputs_inh = inp_inh(input_params.start+1:end);
		all_nn_inputs_bg = inp_bg(input_params.start+1:end);

		if (do_save_inputs > 0),
		 save(fn_inp, 'all_nn_inputs','all_nn_inputs_inh','all_nn_inputs_bg');
		end;

	end;

	%
	% generate the graph(s)
	%
	if (do_gen_graphs == 0),
%OK JSC
		load(fn_graphs, '-mat');

	% graph statistics

	



	else
		[A, E, E_inh, E_bg] = eval(BUILD_GRAPH);
		if (do_save_graphs),
			save(fn_graphs, 'A', 'E', 'E_inh', 'E_bg');
			end;
			
		end;

% background to all neurons
%	E_bg = E_bg + 1;
	E_bg(1:10) = 1;
%	E(1:52) = 1;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% simulation loop
	%
	i_upd = 1;
	dt=1;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	I_syn = zeros(N_nn,T_upd);
	vv=zeros(N_nn,T_upd);
	uu=zeros(N_nn,T_upd);
	gg=zeros(N_nn,T_upd);
	fev = zeros(N_nn,T_upd);

	b = nn_mu_params(:,2);
	v0 = (1/0.08)*(-(5-b.^2) - sqrt(abs((5-b).^2 - 4*0.04*140)));
	v = zeros(N_nn,1) + v0;
	u = nn_mu_params(:,2).*v;

%	g = zeros(N_nn,1) +0.1;
	g = zeros(N_nn,1) +0.01;
	tau = zeros(N_nn,1) +5; % AMPA

	%	e_sig = external input

	all_nn_inputs(1000:sim.T_upd) = 0;

	esig = -all_nn_inputs';


%	all_nn_inputs_inh(1000:sim.T_upd) = 0;
%	all_nn_inputs_inh(1001:2000) = all_nn_inputs_inh;
	all_nn_inputs_inh(2001:3000) = all_nn_inputs_inh(1001:2000);
	all_nn_inputs_inh(3001:4000) = all_nn_inputs_inh(1001:2000);
	all_nn_inputs_inh(4001:5000) = all_nn_inputs_inh(1001:2000);
	esig_inh = -all_nn_inputs_inh';

%	all_nn_inputs_bg(1000:sim.T_upd) = 0;
%	all_nn_inputs_bg(1001:2000) = all_nn_inputs_bg(1001:2000);
	all_nn_inputs_bg(2001:3000) = all_nn_inputs_bg(1001:2000);
	all_nn_inputs_bg(3001:4000) = all_nn_inputs_bg(1001:2000);
	all_nn_inputs_bg(4001:5000) = all_nn_inputs_bg(1001:2000);
	esig_bg = -all_nn_inputs_bg';


	%===========================================
	% inner loop
	%===========================================
	for i=1:T_upd,
		%
		% calc I_syn
		%
		I_s = g .* v;
	
		dg = -g./tau;
	
		g = g + dg*dt;
		gg(:,i)  = g;
	
		I_syn(:,i) = I_s;
		[v,u,ind] = dneuron(v,u, I_s, nn_mu_params(:,1), nn_mu_params(:,2),nn_mu_params(:,3),nn_mu_params(:,4));
		vv(:,i) = v;
		vv(ind,i) = 30;
		uu(:,i) = u;
		
		fev(ind,i) = 1;
	
	%
	% currently: fixed delay
	%
	if (i > sim.delay),
		sf = fev(:,i-sim.delay);
	else
		sf = zeros(sim.N_nn,1);
	end;

	%
	% get pre-synaptic info
	% plus external signals
	pre_fire = w_net .* (A * sf) + ...
	 w_ext .* (E *esig(:,i)) + ...
	 w_ext_inh .* (E_inh * esig_inh(:,i)) + ...
	 w_bg .* (E_bg * esig_bg(:,i));

	g = g + pre_fire;

	end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% evaluate activity
	%
	for nn=1:N_nn,
%	[spi, spt, act] = calc_spiketrain(vm, sim);
	sim.instrument.allvm(i_upd,nn,:) = vv(nn,:);
	sim.instrument.spiketrain(i_upd,nn,:) = fev(nn,:);
			
	end;

fprintf('JSC: fn=%s\n',fn)
%JSC	save(fn);

end;

